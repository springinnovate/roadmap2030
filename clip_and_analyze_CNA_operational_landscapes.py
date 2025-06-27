"""
https://springinnovate.slack.com/archives/C093HTWH4KS/p1750888349804319

* clip the CNA raster to the operational landscapes
* repeat the downstream and travel time analyses for that resulting masks to
  get the total number of people potentially supported by critical natural
  assets within all of our operational landscapes
"""

import logging
import os
import sys

from ecoshard.geoprocessing import routing
from ecoshard import taskgraph
from ecoshard import geoprocessing
from osgeo import gdal
import numpy as np

from pilot_area_downstream_pop_and_es_summary import create_circular_kernel
from pilot_area_downstream_pop_and_es_summary import subset_subwatersheds
from pilot_area_downstream_pop_and_es_summary import calc_flow_dir
from pilot_area_downstream_pop_and_es_summary import mask_by_nonzero_and_sum

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

CNA_RASTER_PATH = "./data/A_90_md5_79f5e0d5d5029d90e8f10d5932da93ff.tif"
PRODSCAPE_VECTOR_PATH = "./data/Prod_scapes_EE_wflags"
GLOBAL_SUBWATERSHEDS_VECTOR_PATH = "./dem_precondition/data/merged_lev06.shp"
DEM_RASTER_PATH = "./dem_precondition/data/astgtm_compressed.tif"

POP_PIXEL_SIZE = [0.008333333333333, -0.008333333333333]
km_in_deg = POP_PIXEL_SIZE[0] * 1000 / 900  # because it's 900m so converting to 1km
BUFFER_SIZE_IN_M = 1000
BUFFER_SIZE_IN_PX = int(np.round(BUFFER_SIZE_IN_M * km_in_deg / POP_PIXEL_SIZE[0]))

POPULATION_RASTER_PATH = "./data/pop_rasters/landscan-global-2023.tif"


WORKSPACE_DIR = "./workspace_clip_and_analyze_CNA"
os.makedirs(WORKSPACE_DIR, exist_ok=True)


def main():
    # reproject Prod_scapes_EE to CNA_RASTER_PATH
    # mask out CNA by prodscape vector

    task_graph = taskgraph.TaskGraph(WORKSPACE_DIR, 3, 10.0)

    cna_raster_info = geoprocessing.get_raster_info(CNA_RASTER_PATH)
    target_projected_prodscape_vector_path = os.path.join(
        WORKSPACE_DIR, "projected_prodscape_vector.gpkg"
    )
    reproject_prodscape_task = task_graph.add_task(
        func=geoprocessing.reproject_vector,
        args=(
            PRODSCAPE_VECTOR_PATH,
            cna_raster_info["projection_wkt"],
            target_projected_prodscape_vector_path,
        ),
        target_path_list=[target_projected_prodscape_vector_path],
        task_name=f"reproject {PRODSCAPE_VECTOR_PATH}",
    )

    target_projected_subwatershed_vector_path = os.path.join(
        WORKSPACE_DIR, "projected_subwatersheds_vector.gpkg"
    )
    reproject_subwatershed_task = task_graph.add_task(
        func=geoprocessing.reproject_vector,
        args=(
            GLOBAL_SUBWATERSHEDS_VECTOR_PATH,
            cna_raster_info["projection_wkt"],
            target_projected_subwatershed_vector_path,
        ),
        target_path_list=[target_projected_subwatershed_vector_path],
        task_name=f"reproject {GLOBAL_SUBWATERSHEDS_VECTOR_PATH}",
    )

    target_masked_cna_raster_path = os.path.join(
        WORKSPACE_DIR, f"masked_{os.path.basename(CNA_RASTER_PATH)}"
    )
    mask_cna_raster_task = task_graph.add_task(
        func=geoprocessing.mask_raster,
        args=(
            (CNA_RASTER_PATH, 1),
            target_projected_prodscape_vector_path,
            target_masked_cna_raster_path,
        ),
        dependent_task_list=[reproject_prodscape_task],
        target_path_list=[target_masked_cna_raster_path],
        task_name=f"mask out {CNA_RASTER_PATH}",
    )

    downstream_subwatershed_vector_path = os.path.join(
        WORKSPACE_DIR, "downstream_subwatersheds.gpkg"
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
        target_path_list=[downstream_subwatershed_vector_path],
        task_name=f"subset {downstream_subwatershed_vector_path}",
    )
    target_clipped_dem_path = os.path.join(
        WORKSPACE_DIR, f"clipped_{os.path.basename(DEM_RASTER_PATH)}"
    )
    target_flow_dir_path = os.path.join(WORKSPACE_DIR, "flow_dir.tif")

    dem_projected_downstream_subwatershed_vector_path = os.path.join(
        WORKSPACE_DIR, "dem_projected_downstream_subwatershed.gpkg"
    )

    dem_raster_info = geoprocessing.get_raster_info(DEM_RASTER_PATH)
    reproject_ds_watershed_task = task_graph.add_task(
        func=geoprocessing.reproject_vector,
        args=(
            downstream_subwatershed_vector_path,
            dem_raster_info["projection_wkt"],
            dem_projected_downstream_subwatershed_vector_path,
        ),
        dependent_task_list=[downstream_subwatershed_task],
        target_path_list=[dem_projected_downstream_subwatershed_vector_path],
        task_name=f"reproject {PRODSCAPE_VECTOR_PATH}",
    )

    flow_dir_task = task_graph.add_task(
        func=calc_flow_dir,
        args=(
            "prodscape_cna",
            DEM_RASTER_PATH,
            dem_projected_downstream_subwatershed_vector_path,
            target_clipped_dem_path,
            target_flow_dir_path,
        ),
        dependent_task_list=[reproject_ds_watershed_task],
        target_path_list=[
            target_clipped_dem_path,
            target_flow_dir_path,
        ],
    )

    aoi_downstream_flow_mask_path = os.path.join(WORKSPACE_DIR, "aoi_ds_coverage.tif")
    flow_accum_task = task_graph.add_task(
        func=routing.flow_accumulation_mfd,
        args=((target_flow_dir_path, 1), aoi_downstream_flow_mask_path),
        kwargs={"weight_raster_path_band": (target_masked_cna_raster_path, 1)},
        dependent_task_list=[flow_dir_task, mask_cna_raster_task],
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
    kernel_path = os.path.join(WORKSPACE_DIR, f"kernel_{BUFFER_SIZE_IN_M}.tif")
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
        f"cna_pop_downstream_of_prod_{BUFFER_SIZE_IN_M}m.tif",
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

    task_graph.join()
    task_graph.close()
    LOGGER.info(f"total sum is: {mask_by_nonzero_task.get()}")


if __name__ == "__main__":
    main()
