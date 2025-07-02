"""
Requires that the target_setting cna_clip_and_analyze.. and
_people_within_travel_time ... are set.
"""

from collections import defaultdict
from datetime import datetime
import logging
import os
import sys

from ecoshard import geoprocessing
from ecoshard import taskgraph
from osgeo import gdal
import pandas as pd
import numpy
import psutil

gdal.SetCacheMax(2**27)

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        "%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s"
        " [%(funcName)s:%(lineno)d] %(message)s"
    ),
    stream=sys.stdout,
)
LOGGER = logging.getLogger(__name__)
logging.getLogger("ecoshard.taskgraph").setLevel(logging.INFO)
logging.getLogger("rasterio").setLevel(logging.WARNING)
logging.getLogger("fiona").setLevel(logging.WARNING)


MASK_RASTER_DICT = {
    "A_25_60m_traveltime": "./workspace_target_setting_dist_to_hab_with_friction/A_25/A_25_max_reach_60min.tif",
    "A_50_60m_traveltime": "./workspace_target_setting_dist_to_hab_with_friction/A_50/A_50_max_reach_60min.tif",
    "A_90_60m_traveltime": "./workspace_target_setting_dist_to_hab_with_friction/A_90/A_90_max_reach_60min.tif",
    "A_25_ds_coverage": "./workspace_clip_and_analyze_CNA/A_25_aoi_ds_coverage_1000m_1000m.tif",
    "A_50_ds_coverage": "./workspace_clip_and_analyze_CNA/A_50_aoi_ds_coverage_1000m_1000m.tif",
    "A_90_ds_coverage": "./workspace_clip_and_analyze_CNA/A_90_aoi_ds_coverage_1000m_1000m.tif",
}

COMBINE_PAIRS = [
    ("A_25_ds_coverage", "A_25_60m_traveltime"),
    ("A_50_ds_coverage", "A_50_60m_traveltime"),
    ("A_90_ds_coverage", "A_90_60m_traveltime"),
]

POPULATION_RASTER_PATH = "./data/pop_rasters/landscan-global-2023.tif"
WORKSPACE_DIR = "./workspace_target_setting_ds_travel_time_combined_pop_analysis"
os.makedirs(WORKSPACE_DIR, exist_ok=True)

COUNTRY_VECTOR_PATH = (
    "data/countries/countries_iso3_md5_6fb2431e911401992e6e56ddf0a9bcda.gpkg"
)
COUNTRY_NAME_FIELD_ID = "iso3"


def _mask_op(value_array, mask_array):
    return numpy.where(mask_array > 0, value_array, 0)


def _merge_mask_op(mask_a, mask_b):
    return (mask_a > 0) | (mask_b > 0)


def main():
    task_graph = taskgraph.TaskGraph(
        WORKSPACE_DIR,
        min(
            psutil.cpu_count(logical=False),
            len(MASK_RASTER_DICT) * len(COMBINE_PAIRS),
        ),
    )
    population_raster_info = geoprocessing.get_raster_info(POPULATION_RASTER_PATH)
    LOGGER.info(population_raster_info["projection_wkt"])
    aligned_raster_task_dict = {}
    for prefix, mask_raster_path in MASK_RASTER_DICT.items():
        aligned_mask_raster_path = os.path.join(WORKSPACE_DIR, f"{prefix}_aligned.tif")
        align_task = task_graph.add_task(
            func=geoprocessing.warp_raster,
            args=(
                mask_raster_path,
                population_raster_info["pixel_size"],
                aligned_mask_raster_path,
                "nearest",
            ),
            kwargs={
                "target_bb": population_raster_info["bounding_box"],
                "target_projection_wkt": population_raster_info["projection_wkt"],
            },
            target_path_list=[aligned_mask_raster_path],
            task_name=f"align {prefix}",
        )
        aligned_raster_task_dict[prefix] = (
            align_task,
            aligned_mask_raster_path,
        )

    for prefix_a, prefix_b in COMBINE_PAIRS:
        task_a, path_a = aligned_raster_task_dict[prefix_a]
        task_b, path_b = aligned_raster_task_dict[prefix_b]
        joined_mask_raster_path = os.path.join(
            WORKSPACE_DIR, f"{prefix_a}_{prefix_b}.tif"
        )
        join_mask_task = task_graph.add_task(
            func=geoprocessing.raster_calculator,
            args=(
                [(path_a, 1), (path_b, 1)],
                _merge_mask_op,
                joined_mask_raster_path,
                gdal.GDT_Byte,
                2,
            ),
            dependent_task_list=[task_a, task_b],
            target_path_list=[joined_mask_raster_path],
            task_name=f"mask pop for {prefix}",
        )
        aligned_raster_task_dict[f"{prefix_a}_{prefix_b}"] = (
            join_mask_task,
            joined_mask_raster_path,
        )

    zonal_task_dict = {}
    for prefix, (
        mask_creation_task,
        mask_raster_path,
    ) in aligned_raster_task_dict.items():
        masked_pop_raster_path = os.path.join(
            WORKSPACE_DIR,
            f"{prefix}_{os.path.basename(POPULATION_RASTER_PATH)}",
        )

        mask_pop_task = task_graph.add_task(
            func=geoprocessing.raster_calculator,
            args=(
                [(POPULATION_RASTER_PATH, 1), (mask_raster_path, 1)],
                _mask_op,
                masked_pop_raster_path,
                gdal.GDT_Float32,
                population_raster_info["nodata"][0],
            ),
            kwargs={
                "allow_different_blocksize": True,
            },
            dependent_task_list=[mask_creation_task],
            target_path_list=[masked_pop_raster_path],
            task_name=f"mask pop for {prefix}",
        )

        zonal_stats_task = task_graph.add_task(
            func=geoprocessing.zonal_statistics,
            args=(
                (masked_pop_raster_path, 1),
                COUNTRY_VECTOR_PATH,
            ),
            kwargs={
                "working_dir": WORKSPACE_DIR,
                "polygons_might_overlap": False,
            },
            store_result=True,
            dependent_task_list=[mask_pop_task],
            task_name=f"calc zonal for {prefix}",
        )
        zonal_task_dict[prefix] = zonal_stats_task

    country_vector = gdal.OpenEx(COUNTRY_VECTOR_PATH, gdal.OF_VECTOR)
    fid_to_country = {
        feature.GetFID(): feature.GetField(COUNTRY_NAME_FIELD_ID)
        for feature in country_vector.GetLayer()
    }
    pop_sum_per_country = defaultdict(dict)
    for prefix, zonal_stats_task in zonal_task_dict.items():
        stats = zonal_stats_task.get()
        col_id = f"pop_sum_{prefix}"
        for fid, value in stats.items():
            pop_sum_per_country[fid_to_country[fid]][col_id] = value["sum"]

    totals = defaultdict(float)
    for cols in pop_sum_per_country.values():
        for col_id, val in cols.items():
            totals[col_id] += val

    pop_sum_per_country["ALL"] = dict(totals)

    df = (
        pd.DataFrame.from_dict(pop_sum_per_country, orient="index")
        .rename_axis("country name")
        .reset_index()
    )

    # put the 'ALL' row first
    df = pd.concat(
        [
            df[df["country name"] == "ALL"],
            df[df["country name"] != "ALL"].sort_values("country name"),
        ],
        ignore_index=True,
    )

    timestamp = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    table_path = f"combined_ds_area_travel_time_mask_pop_count_{timestamp}.csv"
    df.to_csv(table_path, index=False)
    LOGGER.info(f"result in {table_path}")

    task_graph.join()
    task_graph.close()


def combine_masks(raster_path_list, target_raster_path):
    def combine_op(*array_list):
        result = array_list[0] > 0
        for array in array_list[1:]:
            result |= array > 0
        return result

    geoprocessing.raster_calculator(
        [(path, 1) for path in raster_path_list],
        combine_op,
        target_raster_path,
        gdal.GDT_Byte,
        2,
    )


if __name__ == "__main__":
    main()
