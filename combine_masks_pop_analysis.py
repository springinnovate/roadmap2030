"""Combine masks to get population counts."""

from collections import defaultdict
from datetime import datetime
import glob
import logging
import os
import sys

from pilot_area_downstream_pop_and_es_summary import (
    align_and_resize_raster_stack,
)
from pilot_area_downstream_pop_and_es_summary import clean_it
from pilot_area_downstream_pop_and_es_summary import mask_by_nonzero_and_sum
from ecoshard import geoprocessing
from ecoshard import taskgraph
from osgeo import gdal
import pandas as pd

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


MASK_DIRS = [
    "workspace_downstream_es_analysis",
    "workspace_dist_to_hab_with_friction",
]
SUFFIXES = {
    "_max_reach_60min.tif",
    "_aoi_ds_coverage_1000m.tif",
}

POPULATION_RASTER_PATH = "./data/pop_rasters/landscan-global-2023.tif"
WORKSPACE_DIR = "./workspace_combined_masks_pop_analysis"
os.makedirs(WORKSPACE_DIR, exist_ok=True)

matched_rasters = defaultdict(dict)

for raster_dir in MASK_DIRS:
    pattern = os.path.join(raster_dir, "**", "*.tif")
    for file_path in glob.glob(pattern, recursive=True):
        basename = os.path.basename(file_path)
        for suffix in SUFFIXES:
            if basename.endswith(suffix):
                prefix = basename[: -len(suffix)]
                matched_rasters[prefix][suffix] = file_path

# Filter out incomplete sets if necessary
COMPLETE_SETS = {
    prefix: suffix_path_dict.values()
    for prefix, suffix_path_dict in matched_rasters.items()
    if SUFFIXES.issubset(suffix_path_dict)
}

INCOMPLETE_SETS = {
    prefix: suffix_path_dict.values()
    for prefix, suffix_path_dict in matched_rasters.items()
    if not SUFFIXES.issubset(suffix_path_dict)
}


def combine_mask_and_sum_pop(
    prefix,
    path_list,
    workspace_dir,
    population_raster_info,
    population_raster_path,
):
    local_dir = os.path.join(workspace_dir, prefix)
    os.makedirs(local_dir, exist_ok=True)
    LOGGER.info(f"processing {prefix}")
    LOGGER.info(path_list)

    aligned_raster_path_list = [
        clean_it(
            os.path.join(
                local_dir,
                f"%s_{prefix}_aligned%s" % os.path.splitext(os.path.basename(path)),
            )
        )
        for path in path_list
    ]

    align_and_resize_raster_stack(
        path_list,
        aligned_raster_path_list,
        ["nearest"] * len(path_list),
        population_raster_info["pixel_size"],
        "union",
        target_projection_wkt=population_raster_info["projection_wkt"],
    )

    aligned_raster_info = geoprocessing.get_raster_info(aligned_raster_path_list[0])

    aligned_pop_raster_path = os.path.join(
        local_dir, f"aligned_{os.path.basename(population_raster_path)}"
    )

    geoprocessing.warp_raster(
        population_raster_path,
        aligned_raster_info["pixel_size"],
        aligned_pop_raster_path,
        "nearest",
        target_bb=aligned_raster_info["bounding_box"],
        target_projection_wkt=aligned_raster_info["projection_wkt"],
    )

    LOGGER.info(f"look in {local_dir}")

    combined_mask_path = os.path.join(local_dir, f"{prefix}_combined_mask.tif")
    masked_pop_raster_path = os.path.join(local_dir, f"{prefix}_masked_pop.tif")

    combine_masks(aligned_raster_path_list, combined_mask_path)

    results = {}

    # pop count for combined mask
    pop_count_all = mask_by_nonzero_and_sum(
        prefix,
        aligned_pop_raster_path,
        combined_mask_path,
        masked_pop_raster_path,
    )
    results["all"] = pop_count_all

    # pop count for each individual mask
    for mask_path in aligned_raster_path_list:
        mask_name = os.path.splitext(os.path.basename(mask_path))[0]
        indiv_masked_pop_path = os.path.join(local_dir, f"{mask_name}_masked_pop.tif")
        pop_count_indiv = mask_by_nonzero_and_sum(
            mask_name,
            aligned_pop_raster_path,
            mask_path,
            indiv_masked_pop_path,
        )
        results[mask_name] = pop_count_indiv

    return results


def main():
    task_graph = taskgraph.TaskGraph(WORKSPACE_DIR, -1)
    population_raster_info = geoprocessing.get_raster_info(POPULATION_RASTER_PATH)

    results_by_area_task_list = []
    for prefix, path_list in COMPLETE_SETS.items():
        pop_count_task = task_graph.add_task(
            func=combine_mask_and_sum_pop,
            args=(
                prefix,
                path_list,
                WORKSPACE_DIR,
                population_raster_info,
                POPULATION_RASTER_PATH,
            ),
            store_result=True,
            transient_run=True,
            task_name=f"calculate {prefix}",
        )
        results_by_area_task_list.append((prefix, pop_count_task))

    rows = []
    for prefix, pop_dict_task in results_by_area_task_list:
        pop_dict = pop_dict_task.get()
        row = {"area": prefix}
        row.update(pop_dict)
        rows.append(row)

    df = pd.DataFrame(rows)

    timestamp = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    table_path = f"combined_ds_area_travel_time_mask_pop_count_{timestamp}.csv"
    df.to_csv(table_path, index=False)
    LOGGER.info(f"result in {table_path}")


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
