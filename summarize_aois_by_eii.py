"""6/13/205 (3:51pm) (Friday before vacation) Becky needed this an emergency analysis, it's nearly copied and pasted from `doit.py` except it walks a directory path looking for vectors."""

import collections
import csv
import datetime
import glob
import logging
import os
import sys

from ecoshard import geoprocessing
from ecoshard import taskgraph
from osgeo import gdal
from osgeo import osr
from shapely.geometry import shape
import fiona
import numpy


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

PERCENTILES_LIST = [10, 25, 75, 90]
# report the area in the clipped raster that has values >= to these values
THRESHOLD_AREA_LIST = [0.9]

BASE_RASTER_LOOKUP = {
    "eii": r"./data_platform/Nature/eii_padj_v5140524_epsg_3395.tif",
}

ANALYSIS_AOIS = {}
BAD_AOIS = {}  # trash, but you can look if you want

AOI_DIRS = ["./data/WWF-Int_Pilot", "./data/aoi_by_country"]

for aoi_dir in AOI_DIRS:
    for ext in ["shp", "gpkg"]:
        pattern = os.path.join(aoi_dir, "**", f"*.{ext}")
        for file_path in glob.glob(pattern, recursive=True):
            basename = os.path.splitext(os.path.basename(file_path))[0]
            try:
                with fiona.open(file_path, "r") as src:
                    if len(src) == 0:
                        raise ValueError("No features found")

                    valid_geometry_found = False
                    for feature in src:
                        geom = feature["geometry"]
                        if geom and geom != {}:
                            shapely_geom = shape(geom)
                            if not shapely_geom.is_empty:
                                valid_geometry_found = True
                                break

                    if not valid_geometry_found:
                        raise ValueError("No valid geometry found in features")

                ANALYSIS_AOIS[basename] = file_path
            except Exception:
                BAD_AOIS[basename] = file_path


OUTPUT_DIR = "./eii_summary_results"
CLIPPED_DIR = os.path.join(OUTPUT_DIR, "clipped")
for dirpath in [OUTPUT_DIR, CLIPPED_DIR]:
    os.makedirs(dirpath, exist_ok=True)


def vector_area_in_ha(vector_path):
    dataset = gdal.OpenEx(vector_path, gdal.OF_VECTOR)
    layer = dataset.GetLayer()

    source_srs = layer.GetSpatialRef()
    target_srs = osr.SpatialReference()
    target_srs.ImportFromProj4(
        "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +R=6371007 +units=m +no_defs"
    )

    total_area_m2 = 0

    for feature in layer:
        geometry = feature.GetGeometryRef()
        if geometry is not None:
            geom_clone = geometry.Clone()
            geom_clone.AssignSpatialReference(source_srs)
            geom_clone.TransformTo(target_srs)
            total_area_m2 += geom_clone.GetArea()

    total_area_ha = total_area_m2 / 10000
    return total_area_ha


def create_subset(gdf, name, target_vector_path):
    LOGGER.info(f"creating subset of {name}")
    subset_gdf = gdf[gdf["Name"] == name]
    subset_gdf.to_file(target_vector_path, driver="GPKG")
    LOGGER.info(f"done with subset of {name}")


def clip_raster(base_raster_path, summary_vector_path, temp_clip_path):
    base_raster_info = geoprocessing.get_raster_info(base_raster_path)
    summary_vector_info = geoprocessing.get_vector_info(summary_vector_path)
    target_pixel_size = base_raster_info["pixel_size"]
    base_vector_bb = summary_vector_info["bounding_box"]

    target_bb = geoprocessing.transform_bounding_box(
        base_vector_bb,
        summary_vector_info["projection_wkt"],
        base_raster_info["projection_wkt"],
    )

    geoprocessing.warp_raster(
        base_raster_path,
        target_pixel_size,
        temp_clip_path,
        "near",
        target_bb=target_bb,
        vector_mask_options={
            "mask_vector_path": summary_vector_path,
            "all_touched": True,
        },
    )


def get_area_stats(raster_path, thresholds):
    r = gdal.Open(raster_path)
    proj_wkt = r.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(proj_wkt)
    if srs.IsGeographic():
        warp_srs = osr.SpatialReference()
        # Mollweide (World), epsg code does not work
        warp_srs.ImportFromProj4(
            "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +R=6371007 +units=m +no_defs"
        )
        warped_ds = gdal.Warp(
            "",
            r,
            format="MEM",
            dstSRS=warp_srs,
            resampleAlg=gdal.GRA_NearestNeighbour,
        )
    else:
        warped_ds = r

    band = warped_ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype(float)
    nodata = band.GetNoDataValue()

    valid_mask = numpy.ones_like(arr, dtype=bool)
    if nodata is not None:
        valid_mask &= arr != nodata
    valid_mask &= ~numpy.isnan(arr)
    arr = arr[valid_mask]

    gt = warped_ds.GetGeoTransform()
    pixel_area_m2 = abs(gt[1] * gt[5])
    results = {}

    for thr in thresholds:
        pix_count = numpy.count_nonzero(arr >= thr)
        area_ha = (pix_count * pixel_area_m2) / 10000.0  # 1 ha = 10,000 m^2
        results[f"area_ge_{thr}"] = area_ha

    area_ha = (arr.size * pixel_area_m2) / 10000.0  # 1 ha = 10,000 m^2
    results["area_ha"] = area_ha

    return results


def get_stats(raster_path):
    r = gdal.OpenEx(raster_path)
    b = r.GetRasterBand(1)
    array = b.ReadAsArray()
    nodata = b.GetNoDataValue()
    array = array[(array != nodata) & ~numpy.isnan(array)]
    if array.size == 0:
        # guard against a nodata array
        array = numpy.array([0.0])
    stats = {
        "min": numpy.min(array),
        "max": numpy.max(array),
        "sum": numpy.sum(array),
        "mean": numpy.mean(array),
    }

    percentile_dict = {
        f"p{percentile}": value
        for percentile, value in zip(
            PERCENTILES_LIST, numpy.percentile(array, PERCENTILES_LIST)
        )
    }
    stats.update(percentile_dict)

    value_thresholds = get_area_stats(raster_path, THRESHOLD_AREA_LIST)

    stats.update(value_thresholds)

    return stats


def dump_results_to_csv(results, vector_path_lookup, csv_path):
    with open(csv_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            ["vector_id", "raster_name", "area_ha", "min", "max", "mean", "sum"]
            + [f"p{percentile}" for percentile in PERCENTILES_LIST]
            + [f"area_ge_{threshold}" for threshold in THRESHOLD_AREA_LIST]
        )
        for vector_id, info_dict in results.items():
            area_ha = info_dict.get("area_ha", None)
            if area_ha is not None:
                writer.writerow(
                    [
                        vector_id,
                        "",  # raster_name is empty
                        area_ha,
                        "",
                        "",
                        "",
                        "",  # no min/max/mean/sum for area
                    ]
                )

            # Next rows: stats for each raster
            # (anything that's not "area_ha" in `results[vector_id]`)
            for raster_basename, stats_dict in info_dict.items():
                if raster_basename == "area_ha":
                    continue
                if not isinstance(stats_dict, dict):
                    continue

                # Extract stats if they exist
                r_min = stats_dict.get("min", "")
                r_max = stats_dict.get("max", "")
                r_mean = stats_dict.get("mean", "")
                r_sum = stats_dict.get("sum", "")
                area_ha = stats_dict.get("area_ha", "")
                writer.writerow(
                    [
                        vector_id,
                        raster_basename,  # raster_name
                        area_ha,  # area_ha is empty here
                        r_min,
                        r_max,
                        r_mean,
                        r_sum,
                    ]
                    + [
                        stats_dict.get(f"p{percentile}", "")
                        for percentile in PERCENTILES_LIST
                    ]
                    + [
                        stats_dict.get(f"area_ge_{thr}", "")
                        for thr in THRESHOLD_AREA_LIST
                    ]
                )


def main():
    """Entry point."""
    print(os.cpu_count())
    task_graph = taskgraph.TaskGraph(
        OUTPUT_DIR, os.cpu_count(), reporting_interval=10.0
    )
    results = collections.defaultdict(lambda: collections.defaultdict(dict))
    for vector_id, vector_path in ANALYSIS_AOIS.items():
        results[vector_id]["area_ha"] = vector_area_in_ha(vector_path)
        LOGGER.info(f"processing {vector_id}")
        for raster_basename, raster_path in BASE_RASTER_LOOKUP.items():
            if not os.path.exists(raster_path):
                raise RuntimeError(f"{raster_path} not found")
            LOGGER.info(f"clipping {raster_basename} to {vector_id}")
            clipped_raster_path = os.path.join(
                CLIPPED_DIR, f"{vector_id}_{raster_basename}.tif"
            )
            clipped_task = task_graph.add_task(
                func=clip_raster,
                args=(raster_path, vector_path, clipped_raster_path),
                ignore_path_list=[vector_path],
                target_path_list=[clipped_raster_path],
                task_name=f"clipping {raster_path} to {vector_path}",
            )

            stats_task = task_graph.add_task(
                func=get_stats,
                args=(clipped_raster_path,),
                dependent_task_list=[clipped_task],
                store_result=True,
                task_name=f"stats for {raster_path}",
            )
            results[vector_id][raster_basename]["stats"] = stats_task

    for vector_id in ANALYSIS_AOIS:
        for raster_basename in BASE_RASTER_LOOKUP:
            stats_task = results[vector_id][raster_basename]["stats"]
            results[vector_id][raster_basename] = {}
            for fn_str, value in stats_task.get().items():
                results[vector_id][raster_basename][fn_str] = value

    task_graph.join()
    print(results)
    timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    output_filename = f"eii_summary_results_{timestamp}.csv"
    dump_results_to_csv(results, ANALYSIS_AOIS, output_filename)
    print(f"all done -- results in {output_filename}!")


if __name__ == "__main__":
    main()
