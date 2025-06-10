"""
Run with docker container:

docker run --rm -it -v .:/usr/local/esos_c_models -v "D:\repositories\dem_precondition":/usr/local/esos_c_models/dem_precondition therealspring/roadmap2030_executor:latest
"""

import csv
import collections
import datetime
import logging
import os
import sys

from ecoshard import geoprocessing
from ecoshard import taskgraph
from ecoshard.geoprocessing import routing
from osgeo import gdal
from shapely.geometry import box
from shapely.ops import transform
import fiona
import geopandas as gpd
import numpy
import pandas as pd
import pyproj

logging.basicConfig(
    level=logging.DEBUG,
    stream=sys.stdout,
    format=(
        "%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s"
        " [%(pathname)s.%(funcName)s:%(lineno)d] %(message)s"
    ),
)
LOGGER = logging.getLogger(__name__)
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
logging.getLogger("PIL").setLevel(logging.ERROR)
logging.getLogger("ecoshard.taskgraph").setLevel(logging.INFO)
logging.getLogger("fiona").setLevel(logging.WARN)

POPULATION_RASTERS = [
    "./data/pop_rasters/landscan-global-2023.tif",
]

POP_PIXEL_SIZE = [0.008333333333333, -0.008333333333333]
km_in_deg = POP_PIXEL_SIZE[0] * 1000 / 900  # because it's 900m so converting to 1km
BUFFER_AMOUNTS_IN_PIXELS_M = [
    (int(numpy.round(x * km_in_deg / POP_PIXEL_SIZE[0])), x * 1000) for x in range(0, 6)
]  # try buffer zones of 0-5km

# This is relative because Docker will map a volume
GLOBAL_SUBWATERSHEDS_VECTOR_PATH = "./dem_precondition/data/merged_lev06.shp"
DEM_RASTER_PATH = "./dem_precondition/data/astgtm_compressed.tif"
ANALYSIS_AOIS = {
    "C217_Bengo_Territorios": "./data/WWF-Int_Pilot/Colombia/Colombia/C217_Bengo_Territorios/C217_Bengo_Territorios.shp",
}


OUTPUT_DIR = "./workspace_downstream_es_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)


def calc_flow_dir(
    analysis_id, base_dem_raster_path, aoi_vector_path, target_flow_dir_path
):
    local_workspace_dir = os.path.dirname(target_flow_dir_path)
    os.makedirs(local_workspace_dir, exist_ok=True)
    clipped_dem_path = os.path.join(local_workspace_dir, f"{analysis_id}_dem.tif")
    dem_info = geoprocessing.get_raster_info(base_dem_raster_path)
    aoi_ds = gdal.OpenEx(aoi_vector_path, gdal.OF_VECTOR | gdal.GA_ReadOnly)
    aoi_layer = aoi_ds.GetLayer()
    minx, maxx, miny, maxy = aoi_layer.GetExtent()
    aoi_bb = [minx, miny, maxx, maxy]
    bounding_box = geoprocessing.merge_bounding_box_list(
        [dem_info["bounding_box"], aoi_bb], "intersection"
    )
    geoprocessing.warp_raster(
        base_dem_raster_path,
        POP_PIXEL_SIZE,
        clipped_dem_path,
        "nearest",
        target_bb=bounding_box,
        vector_mask_options={
            "mask_vector_path": aoi_vector_path,
            "all_touched": True,
            "target_mask_value": 0,
        },
    )
    r = gdal.OpenEx(clipped_dem_path, gdal.OF_RASTER | gdal.OF_UPDATE)
    b = r.GetRasterBand(1)
    b.SetNoDataValue(0)
    b = None
    r = None

    filled_dem_path = os.path.join(local_workspace_dir, f"{analysis_id}_dem_filled.tif")
    routing.fill_pits(
        (clipped_dem_path, 1),
        filled_dem_path,
        working_dir=local_workspace_dir,
        max_pixel_fill_count=10000,
    )

    routing.flow_dir_mfd(
        (filled_dem_path, 1),
        target_flow_dir_path,
        working_dir=local_workspace_dir,
    )

    return clipped_dem_path


def rasterize(aoi_vector_path, dem_raster_path, aoi_raster_mask_path):
    geoprocessing.new_raster_from_base(
        dem_raster_path,
        aoi_raster_mask_path,
        datatype=gdal.GDT_Byte,
        band_nodata_list=[-1],
    )
    geoprocessing.rasterize(aoi_vector_path, aoi_raster_mask_path, burn_values=[1])


def mask_by_nonzero_and_sum(
    analysis_id, base_raster_path, mask_raster_path, target_masked_path
):
    base_raster_info = geoprocessing.get_raster_info(base_raster_path)
    mask_raster_info = geoprocessing.get_raster_info(mask_raster_path)
    nodata = base_raster_info["nodata"][0]

    def _mask_by_nonzero(base_array, mask_array):
        result = base_array.copy()
        result[mask_array <= 0] = nodata
        return result

    working_dir_path = os.path.dirname(target_masked_path)
    aligned_raster_path_list = [
        os.path.join(
            working_dir_path,
            analysis_id,
            f"%s_{analysis_id}_aligned%s" % os.path.splitext(os.path.basename(path)),
        )
        for path in [base_raster_path, mask_raster_path]
    ]
    geoprocessing.align_and_resize_raster_stack(
        [base_raster_path, mask_raster_path],
        aligned_raster_path_list,
        ["near"] * 2,
        base_raster_info["pixel_size"],
        mask_raster_info["bounding_box"],
    )

    geoprocessing.raster_calculator(
        [(path, 1) for path in aligned_raster_path_list],
        _mask_by_nonzero,
        target_masked_path,
        base_raster_info["datatype"],
        nodata,
        allow_different_blocksize=True,
        skip_sparse=True,
    )

    array = gdal.OpenEx(target_masked_path).ReadAsArray()
    array = array[array != nodata]
    return numpy.sum(array)


def create_circular_kernel(kernel_path, buffer_size_in_px):
    diameter = buffer_size_in_px * 2 + 1
    kernel_array = numpy.zeros((diameter, diameter), dtype=numpy.float32)
    cx, cy = buffer_size_in_px, buffer_size_in_px

    for i in range(diameter):
        for j in range(diameter):
            if (i - cx) ** 2 + (j - cy) ** 2 <= buffer_size_in_px**2:
                kernel_array[i, j] = 1.0

    driver = gdal.GetDriverByName("GTiff")
    out_raster = driver.Create(kernel_path, diameter, diameter, 1, gdal.GDT_Float32)
    out_raster.GetRasterBand(1).WriteArray(kernel_array)
    out_raster.FlushCache()
    out_raster = None


import geopandas as gpd
import fiona
import pyproj
from shapely.geometry import box
from shapely.ops import transform


def subset_subwatersheds(
    aoi_vector_path, subwatershed_vector_path, subset_subwatersheds_vector_path
):
    # Prepare AOI
    aoi_vector = gpd.read_file(aoi_vector_path)
    aoi_vector.geometry = aoi_vector.geometry.buffer(0)
    aoi_crs = aoi_vector.crs
    aoi_union = aoi_vector.geometry.union_all()
    aoi_bbox_geom = box(*aoi_vector.total_bounds)

    # Retrieve subwatershed CRS
    with fiona.open(subwatershed_vector_path, "r") as subwatershed_vector:
        subwatershed_crs = subwatershed_vector.crs

    if aoi_crs != subwatershed_crs:
        transformer = pyproj.Transformer.from_crs(
            aoi_crs, subwatershed_crs, always_xy=True
        ).transform
        aoi_bbox_geom = transform(transformer, aoi_bbox_geom)

    # Initial filter based on bbox for fast lookup, then slower geometry
    # itersection
    subwatershed_filtered = gpd.read_file(
        subwatershed_vector_path, bbox=aoi_bbox_geom.bounds
    )
    initial_subwatersheds = subwatershed_filtered[
        subwatershed_filtered.intersects(aoi_union)
    ]

    all_hybas_ids = set(initial_subwatersheds["HYBAS_ID"])
    downstream_ids = set(initial_subwatersheds["NEXT_DOWN"]) - {0}

    # Create lookup of ID -> NEXT_DOWN
    with fiona.open(subwatershed_vector_path, "r") as src:
        hybas_to_nextdown = {
            f["properties"]["HYBAS_ID"]: f["properties"]["NEXT_DOWN"] for f in src
        }

    # breadth first graph walk to pick up all the ids
    while downstream_ids:
        all_hybas_ids.update(downstream_ids)
        downstream_ids = (
            {
                hybas_to_nextdown[hybas_id]
                for hybas_id in downstream_ids
                if hybas_id in hybas_to_nextdown
            }
            - all_hybas_ids
            - {0}
        )

    # Load all relevant subwatersheds in one go
    with fiona.open(subwatershed_vector_path, "r") as subwatershed_vector:
        fetched_subwatersheds = [
            f
            for f in subwatershed_vector
            if f["properties"]["HYBAS_ID"] in all_hybas_ids
        ]

    all_subwatersheds_gdf = gpd.GeoDataFrame.from_features(
        fetched_subwatersheds, crs=subwatershed_crs
    )
    all_subwatersheds_gdf.geometry = all_subwatersheds_gdf.geometry.buffer(0)

    if subwatershed_crs != aoi_crs:
        all_subwatersheds_gdf = all_subwatersheds_gdf.to_crs(aoi_crs)

    all_subwatersheds_gdf.to_file(subset_subwatersheds_vector_path, driver="GPKG")


def main():
    """Entry point."""
    task_graph = taskgraph.TaskGraph(
        OUTPUT_DIR, os.cpu_count(), reporting_interval=10.0
    )
    kernel_task_map = {}
    result = collections.defaultdict(
        lambda: collections.defaultdict(lambda: collections.defaultdict(dict))
    )
    file = open("log.txt", "w")
    clipped_dem_work_list = []
    for analysis_id, aoi_vector_path in ANALYSIS_AOIS.items():
        subset_subwatersheds_vector_path = None

        local_workspace_dir = os.path.join(OUTPUT_DIR, analysis_id)
        os.makedirs(local_workspace_dir, exist_ok=True)

        aoi_raster_mask_path = os.path.join(
            local_workspace_dir, f"{analysis_id}_aoi_mask.tif"
        )
        target_projection_wkt = geoprocessing.get_raster_info(DEM_RASTER_PATH)[
            "projection_wkt"
        ]
        reprojected_aoi_vector_path = (
            "%s_projected.gpkg" % os.path.splitext(aoi_raster_mask_path)[0]
        )
        reproject_task = task_graph.add_task(
            func=geoprocessing.reproject_vector,
            args=(
                aoi_vector_path,
                target_projection_wkt,
                reprojected_aoi_vector_path,
            ),
            ignore_path_list=[aoi_vector_path, reprojected_aoi_vector_path],
            target_path_list=[reprojected_aoi_vector_path],
            task_name=f"reproject {analysis_id}",
        )

        subset_subwatersheds_vector_path = os.path.join(
            local_workspace_dir, f"subwatershed_{analysis_id}.gpkg"
        )
        subset_task = task_graph.add_task(
            func=subset_subwatersheds,
            args=(
                reprojected_aoi_vector_path,
                GLOBAL_SUBWATERSHEDS_VECTOR_PATH,
                subset_subwatersheds_vector_path,
            ),
            dependent_task_list=[reproject_task],
            target_path_list=[subset_subwatersheds_vector_path],
            task_name=f"subset subwatersheds for {analysis_id}",
        )
        task_graph.join()
        LOGGER.debug(
            f"all done with testing, look in {subset_subwatersheds_vector_path}"
        )
        return
        flow_dir_path = os.path.join(
            local_workspace_dir, f"{analysis_id}_mfd_flow_dir.tif"
        )
        flow_dir_task = task_graph.add_task(
            func=calc_flow_dir,
            args=(
                analysis_id,
                dem_raster_path,
                subset_subwatersheds_vector_path,
                flow_dir_path,
            ),
            dependent_task_list=[subset_task],
            target_path_list=[flow_dir_path],
            store_result=True,
            task_name=f"calculate flow dir for {analysis_id}",
        )
        clipped_dem_work_list.append(
            (
                analysis_id,
                (
                    local_workspace_dir,
                    flow_dir_path,
                    reprojected_aoi_vector_path,
                    aoi_raster_mask_path,
                    flow_dir_task,
                ),
            )
        )
    for analysis_id, (
        local_workspace_dir,
        flow_dir_path,
        reprojected_aoi_vector_path,
        aoi_raster_mask_path,
        flow_dir_task,
    ) in clipped_dem_work_list:
        clipped_dem_path = flow_dir_task.get()
        print(analysis_id)
        rasterize_task = task_graph.add_task(
            func=rasterize,
            args=(
                reprojected_aoi_vector_path,
                clipped_dem_path,
                aoi_raster_mask_path,
            ),
            dependent_task_list=[reproject_task],
            ignore_path_list=[reprojected_aoi_vector_path],
            target_path_list=[aoi_raster_mask_path],
            task_name=f"{analysis_id} raster mask",
        )

        aoi_downstream_flow_mask_path = os.path.join(
            local_workspace_dir, f"{analysis_id}_aoi_ds_coverage.tif"
        )
        flow_accum_task = task_graph.add_task(
            func=routing.flow_accumulation_mfd,
            args=((flow_dir_path, 1), aoi_downstream_flow_mask_path),
            kwargs={"weight_raster_path_band": (aoi_raster_mask_path, 1)},
            dependent_task_list=[flow_dir_task, rasterize_task],
            target_path_list=[aoi_downstream_flow_mask_path],
            task_name=f"flow accum for {analysis_id}",
        )

        for buffer_size_in_px, buffer_size_in_m in BUFFER_AMOUNTS_IN_PIXELS_M:
            if buffer_size_in_px > 0:
                buffered_downstream_flow_mask_path = (
                    f"%s_{buffer_size_in_m}m%s"
                    % os.path.splitext(aoi_downstream_flow_mask_path)
                )
                # make a kernel raster that is a circle kernel that's all 1s within buffer_size_in_px from the center
                # dimensions should be buffer_size_in_px*2+1 X buffer_size_in_px*2+1
                # no projection is necessary
                kernel_path = os.path.join(
                    local_workspace_dir, f"kernel_{buffer_size_in_m}.tif"
                )
                convolve_dependent_task_list = [flow_accum_task]
                if kernel_path not in kernel_task_map:
                    kernel_task = task_graph.add_task(
                        func=create_circular_kernel,
                        args=(kernel_path, buffer_size_in_px),
                        target_path_list=[kernel_path],
                        task_name=f"kernel for {kernel_path}",
                    )
                    kernel_task_map[kernel_path] = kernel_task
                    convolve_dependent_task_list.append(kernel_task)
                else:
                    convolve_dependent_task_list.append(kernel_task_map[kernel_path])
                buffer_task = task_graph.add_task(
                    func=geoprocessing.convolve_2d,
                    args=(
                        (aoi_downstream_flow_mask_path, 1),
                        (kernel_path, 1),
                        buffered_downstream_flow_mask_path,
                    ),
                    kwargs={"n_workers": 1},
                    target_path_list=[buffered_downstream_flow_mask_path],
                    dependent_task_list=convolve_dependent_task_list,
                    task_name=f"buffer {buffered_downstream_flow_mask_path}",
                )
                mask_by_nonzero_and_sum_dependent_task_list = [buffer_task]
            else:
                buffered_downstream_flow_mask_path = aoi_downstream_flow_mask_path
                mask_by_nonzero_and_sum_dependent_task_list = [flow_accum_task]

            for population_raster_path in POPULATION_RASTERS:
                pop_basename = os.path.splitext(
                    os.path.basename(population_raster_path)
                )[0]
                print(f"processing {pop_basename} for {buffer_size_in_m}m buffer")
                masked_population_raster_path = os.path.join(
                    local_workspace_dir,
                    f"{analysis_id}_{buffer_size_in_m}m_{os.path.basename(population_raster_path)}",
                )
                mask_by_nonzero_task = task_graph.add_task(
                    func=mask_by_nonzero_and_sum,
                    args=(
                        f"{analysis_id}_{pop_basename}_{buffer_size_in_m}",
                        population_raster_path,
                        buffered_downstream_flow_mask_path,
                        masked_population_raster_path,
                    ),
                    dependent_task_list=mask_by_nonzero_and_sum_dependent_task_list,
                    target_path_list=[masked_population_raster_path],
                    store_result=True,
                    task_name=f"{analysis_id}_population mask",
                )
                file.write(
                    f"{analysis_id}_{pop_basename}, {population_raster_path}, {buffered_downstream_flow_mask_path}, {masked_population_raster_path}\n"
                )

                result[analysis_id][pop_basename][
                    buffer_size_in_m
                ] = mask_by_nonzero_task

    for analysis_id in result:
        for pop_basename in result[analysis_id]:
            for buffer_size_in_m in result[analysis_id][pop_basename]:
                result[analysis_id][pop_basename][buffer_size_in_m] = result[
                    analysis_id
                ][pop_basename][buffer_size_in_m].get()

    timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    output_filename = os.path.join(OUTPUT_DIR, f"pop_results_{timestamp}.csv")
    with open(output_filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["analysis_id", "pop_basename", "buffer in m", "value"])
        for analysis_id in result:
            for pop_basename, buffer_val_dict in result[analysis_id].items():
                for buffer, val in buffer_val_dict.items():
                    writer.writerow([analysis_id, pop_basename, buffer, val])
    task_graph.join()
    task_graph.close()
    LOGGER.info(f"all done results at {output_filename}")
    file.close()


if __name__ == "__main__":
    main()
