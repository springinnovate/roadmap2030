"""
https://springinnovate.slack.com/archives/C093HTWH4KS/p1750888349804319

* clip the CNA raster to the operational landscapes
* repeat the downstream and travel time analyses for that resulting masks to
  get the total number of people potentially supported by critical natural
  assets within all of our operational landscapes
"""

from datetime import datetime
import logging
import os
import sys

from ecoshard import geoprocessing
from ecoshard import taskgraph
from ecoshard.geoprocessing import routing
from osgeo import gdal
from rasterio.warp import reproject, Resampling
import numpy as np
import pandas as pd
import rasterio

from pilot_indicators_area_downstream_pop_and_es_summary import (
    create_circular_kernel,
)
from pilot_indicators_area_downstream_pop_and_es_summary import (
    subset_subwatersheds,
)
from pilot_indicators_area_downstream_pop_and_es_summary import calc_flow_dir
from pilot_indicators_area_downstream_pop_and_es_summary import (
    mask_by_nonzero_and_sum,
)

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

CNA_RASTER_PATH_DICT = {
    "A_25": "./data/A_25_md5_737a2625eda959bc0e9c024919e5e7b9.tif",
    "A_50": "./data/A_50_md5_e26abccbc56f1188e269458fe427073f.tif",
    "A_90": "./data/A_90_md5_79f5e0d5d5029d90e8f10d5932da93ff.tif",
}
# Prod_scapes_EE_wflags is the "operational landscapes"
PRODSCAPE_VECTOR_PATH = "./data/Prod_scapes_EE_wflags/Prod_scapes_EE_wflags_fixed.gpkg"
GLOBAL_SUBWATERSHEDS_VECTOR_PATH = "./dem_precondition/data/merged_lev06.shp"
DEM_RASTER_PATH = "./dem_precondition/data/astgtm_compressed.tif"

POP_PIXEL_SIZE_IN_DEG = 0.008333333333333
BUFFER_SIZE_IN_M = 1000
BUFFER_SIZE_IN_PX = int(np.ceil(BUFFER_SIZE_IN_M / 111000 * POP_PIXEL_SIZE_IN_DEG))

POPULATION_RASTER_PATH = "./data/pop_rasters/landscan-global-2023.tif"

COUNTRY_VECTOR_PATH = (
    "data/countries/countries_iso3_md5_6fb2431e911401992e6e56ddf0a9bcda.gpkg"
)
COUNTRY_NAME_FIELD_ID = "iso3"

GLOBAL_TRAVEL_TIME_MASK_1HR_PATH = "."

WORKSPACE_DIR = "./workspace_clip_and_analyze_CNA"
os.makedirs(WORKSPACE_DIR, exist_ok=True)


def warp_to_match(
    base_raster_path,
    warp_to_this_raster_path,
    target_raster_path,
    resampling=Resampling.nearest,
):
    with rasterio.open(warp_to_this_raster_path) as target:
        target_meta = target.meta.copy()

    target_meta.update(
        {
            "tiled": True,
            "blockxsize": 256,
            "blockysize": 256,
        }
    )

    # Initialize destination array based on target metadata
    destination = np.full(
        (target_meta["height"], target_meta["width"]),
        target_meta["nodata"],
        dtype=target_meta["dtype"],
    )

    with rasterio.open(base_raster_path) as src:
        reproject(
            source=rasterio.band(src, 1),
            destination=destination,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=target_meta["transform"],
            dst_crs=target_meta["crs"],
            resampling=resampling,
            dst_nodata=target_meta["nodata"],
        )

    # Write the reprojected data to the new raster file
    with rasterio.open(target_raster_path, "w", **target_meta) as dst:
        dst.write(destination, 1)


def main():
    LOGGER.info("starting")
    task_graph = taskgraph.TaskGraph(WORKSPACE_DIR, 4, 10.0)

    zonal_stats_task_dict = {}
    for cna_label, cna_raster_path in CNA_RASTER_PATH_DICT.items():
        zonal_stats_task, mask_by_nonzero_task = process_cna_raster(
            cna_label, task_graph, cna_raster_path
        )
        zonal_stats_task_dict[cna_label] = (
            zonal_stats_task,
            mask_by_nonzero_task,
        )

    df_final = None
    for cna_label, (
        zonal_stats_task,
        mask_by_nonzero_task,
    ) in zonal_stats_task_dict.items():
        results_by_country = zonal_stats_task.get()
        results_by_country.insert(
            0, {"country": "ALL", "pop sum": mask_by_nonzero_task.get()}
        )

        temp_df = pd.DataFrame(results_by_country)

        # Rename pop sum column with cna_label
        temp_df = temp_df.rename(columns={"pop sum": f"pop sum_{cna_label}"})

        if df_final is None:
            df_final = temp_df
        else:
            # Merge on 'country' column to ensure correct alignment
            df_final = pd.merge(df_final, temp_df, on="country", how="outer")
    timestamp = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    table_path = f"CNA_pop_by_country_{timestamp}.csv"
    df_final.to_csv(table_path, index=False)
    LOGGER.info(f"result at {table_path}")

    task_graph.join()
    task_graph.close()


def process_cna_raster(job_id, task_graph, cna_raster_path):
    # reproject Prod_scapes_EE to cna_raster_path
    # mask out CNA by prodscape vector

    cna_raster_info = geoprocessing.get_raster_info(cna_raster_path)
    target_projected_prodscape_vector_path = os.path.join(
        WORKSPACE_DIR, f"{job_id}_projected_prodscape_vector.gpkg"
    )
    reproject_prodscape_task = task_graph.add_task(
        func=geoprocessing.reproject_vector,
        args=(
            PRODSCAPE_VECTOR_PATH,
            cna_raster_info["projection_wkt"],
            target_projected_prodscape_vector_path,
        ),
        ignore_path_list=[
            PRODSCAPE_VECTOR_PATH,
            target_projected_prodscape_vector_path,
        ],
        target_path_list=[target_projected_prodscape_vector_path],
        task_name=f"reproject {PRODSCAPE_VECTOR_PATH}",
    )

    target_projected_subwatershed_vector_path = os.path.join(
        WORKSPACE_DIR, f"{job_id}_projected_subwatersheds_vector.gpkg"
    )
    reproject_subwatershed_task = task_graph.add_task(
        func=geoprocessing.reproject_vector,
        args=(
            GLOBAL_SUBWATERSHEDS_VECTOR_PATH,
            cna_raster_info["projection_wkt"],
            target_projected_subwatershed_vector_path,
        ),
        ignore_path_list=[
            GLOBAL_SUBWATERSHEDS_VECTOR_PATH,
            target_projected_subwatershed_vector_path,
        ],
        target_path_list=[target_projected_subwatershed_vector_path],
        task_name=f"reproject {GLOBAL_SUBWATERSHEDS_VECTOR_PATH}",
    )

    target_masked_cna_raster_path = os.path.join(
        WORKSPACE_DIR, f"{job_id}_masked_{os.path.basename(cna_raster_path)}"
    )
    mask_cna_raster_task = task_graph.add_task(
        func=geoprocessing.mask_raster,
        args=(
            (cna_raster_path, 1),
            target_projected_prodscape_vector_path,
            target_masked_cna_raster_path,
        ),
        dependent_task_list=[reproject_prodscape_task],
        target_path_list=[target_masked_cna_raster_path],
        task_name=f"mask out {cna_raster_path}",
    )

    downstream_subwatershed_vector_path = os.path.join(
        WORKSPACE_DIR, f"{job_id}_downstream_subwatersheds.gpkg"
    )
    downstream_subwatershed_task = task_graph.add_task(
        func=subset_subwatersheds,
        args=(
            target_projected_prodscape_vector_path,
            target_projected_subwatershed_vector_path,
            downstream_subwatershed_vector_path,
        ),
        dependent_task_list=[
            reproject_subwatershed_task,
            reproject_prodscape_task,
        ],
        ignore_path_list=[
            target_projected_prodscape_vector_path,
            target_projected_subwatershed_vector_path,
            downstream_subwatershed_vector_path,
        ],
        target_path_list=[downstream_subwatershed_vector_path],
        task_name=f"subset {downstream_subwatershed_vector_path}",
    )
    target_clipped_dem_path = os.path.join(
        WORKSPACE_DIR, f"{job_id}_clipped_{os.path.basename(DEM_RASTER_PATH)}"
    )
    target_flow_dir_path = os.path.join(WORKSPACE_DIR, f"{job_id}_flow_dir.tif")

    dem_projected_downstream_subwatershed_vector_path = os.path.join(
        WORKSPACE_DIR, f"{job_id}_dem_projected_downstream_subwatershed.gpkg"
    )

    dem_raster_info = geoprocessing.get_raster_info(DEM_RASTER_PATH)
    reproject_ds_watershed_task = task_graph.add_task(
        func=geoprocessing.reproject_vector,
        args=(
            downstream_subwatershed_vector_path,
            dem_raster_info["projection_wkt"],
            dem_projected_downstream_subwatershed_vector_path,
        ),
        ignore_path_list=[dem_projected_downstream_subwatershed_vector_path],
        dependent_task_list=[downstream_subwatershed_task],
        target_path_list=[dem_projected_downstream_subwatershed_vector_path],
        task_name=f"reproject {PRODSCAPE_VECTOR_PATH}",
    )

    flow_dir_task = task_graph.add_task(
        func=calc_flow_dir,
        args=(
            f"prodscape_cna_{job_id}",
            DEM_RASTER_PATH,
            dem_projected_downstream_subwatershed_vector_path,
            target_clipped_dem_path,
            target_flow_dir_path,
        ),
        ignore_path_list=[dem_projected_downstream_subwatershed_vector_path],
        dependent_task_list=[reproject_ds_watershed_task],
        target_path_list=[
            target_clipped_dem_path,
            target_flow_dir_path,
        ],
        task_name="calc flow dir",
    )

    warped_masked_cna_raster_path = os.path.join(
        WORKSPACE_DIR,
        f"{job_id}_warped_{os.path.basename(target_masked_cna_raster_path)}",
    )
    warp_cna_task = task_graph.add_task(
        func=warp_to_match,
        args=(
            target_masked_cna_raster_path,
            target_flow_dir_path,
            warped_masked_cna_raster_path,
        ),
        dependent_task_list=[flow_dir_task, mask_cna_raster_task],
        target_path_list=[warped_masked_cna_raster_path],
        task_name=f"warp {warped_masked_cna_raster_path}",
    )

    aoi_downstream_flow_mask_path = os.path.join(
        WORKSPACE_DIR, f"{job_id}_aoi_ds_coverage_{BUFFER_SIZE_IN_M}m.tif"
    )
    flow_accum_task = task_graph.add_task(
        func=routing.flow_accumulation_mfd,
        args=((target_flow_dir_path, 1), aoi_downstream_flow_mask_path),
        kwargs={"weight_raster_path_band": (warped_masked_cna_raster_path, 1)},
        dependent_task_list=[
            warp_cna_task,
            flow_dir_task,
        ],
        target_path_list=[aoi_downstream_flow_mask_path],
        task_name="downstream coverage",
    )

    # convolve it out 1km
    buffered_downstream_flow_mask_path = f"%s_{BUFFER_SIZE_IN_M}m%s" % os.path.splitext(
        aoi_downstream_flow_mask_path
    )
    # make a kernel raster that is a circle kernel that's all 1s
    # within buffer_size_in_px from the center
    # dimensions should be
    #    buffer_size_in_px*2+1 X buffer_size_in_px*2+1
    kernel_path = os.path.join(WORKSPACE_DIR, f"{job_id}_kernel_{BUFFER_SIZE_IN_M}.tif")
    kernel_task = task_graph.add_task(
        func=create_circular_kernel,
        args=(kernel_path, BUFFER_SIZE_IN_PX),
        target_path_list=[kernel_path],
        task_name=f"kernel for {kernel_path}",
    )

    buffer_task = task_graph.add_task(
        func=geoprocessing.convolve_2d,
        args=(
            (aoi_downstream_flow_mask_path, 1),
            (kernel_path, 1),
            buffered_downstream_flow_mask_path,
        ),
        kwargs={"n_workers": 1},
        target_path_list=[buffered_downstream_flow_mask_path],
        dependent_task_list=[flow_accum_task, kernel_task],
        task_name=f"buffer {buffered_downstream_flow_mask_path}",
    )

    masked_population_raster_path = os.path.join(
        WORKSPACE_DIR,
        f"{job_id}_cna_pop_downstream_of_prod_{BUFFER_SIZE_IN_M}m.tif",
    )
    mask_by_nonzero_task = task_graph.add_task(
        func=mask_by_nonzero_and_sum,
        args=(
            "cna_pop",
            POPULATION_RASTER_PATH,
            buffered_downstream_flow_mask_path,
            masked_population_raster_path,
        ),
        dependent_task_list=[buffer_task],
        target_path_list=[masked_population_raster_path],
        store_result=True,
        task_name=f"ds_cna_population mask {BUFFER_SIZE_IN_M}",
    )

    zonal_stats_task = task_graph.add_task(
        func=zonal_stats,
        args=(
            masked_population_raster_path,
            COUNTRY_VECTOR_PATH,
            COUNTRY_NAME_FIELD_ID,
        ),
        store_result=True,
        dependent_task_list=[mask_by_nonzero_task],
        task_name="calc raster stats by country",
    )

    LOGGER.info(
        f"\n{target_masked_cna_raster_path} -- CNA raster masked by prodscapes"
        f"\n{aoi_downstream_flow_mask_path} -- areas downstream of the CNA masked prodscapes (any positive value is downstream)"
        f"\n{masked_population_raster_path} -- population masked by the downstream areas of the CNA masked by prodcsapes (values are people counts per pixel)"
    )
    return zonal_stats_task, mask_by_nonzero_task


def zonal_stats(base_raster_path, vector_path, field_name):
    zonal_stats_by_fid = geoprocessing.zonal_statistics(
        (base_raster_path, 1),
        vector_path,
    )
    vector = gdal.OpenEx(vector_path, gdal.OF_VECTOR)
    layer = vector.GetLayer()

    zonal_results = []

    for feature in layer:
        fid = feature.GetFID()
        feature_id = feature.GetField(field_name)
        zonal_results.append(
            {"country": feature_id, "pop sum": zonal_stats_by_fid[fid]["sum"]}
        )

    return zonal_results


if __name__ == "__main__":
    main()
