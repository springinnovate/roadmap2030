import datetime
import re
import glob
import logging
import numpy as np
import os
import sys
import traceback

from ecoshard import taskgraph
from ecoshard import geoprocessing
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
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

AOI_PATH, AOI_NAME_KEY = (
    # r"D:\repositories\roadmap2030\data\aois\SOKNOT_landscapes.gpkg",
    # "Sublandsca",
    r"D:\repositories\data_platform\Conservation_Activities\SOKNOT\SOKNOT_ProjectAreas.gpkg",
    "Name",
)


# Make a list here if you like
RASTER_PATHS_TO_SUMMARIZE = glob.glob(
    r"D:\repositories\data_platform\Nature\eii_soknot"
) + glob.glob(r"Z:/data_platform/Nature/global_pasture_watch_rasters/*.tif")
# r"Z:\data_platform\Nature\global_pasture_watch_rasters\gpw_gpp.daily.grass_lue.model_m_30m_s_20000101_20000228_go_epsg.4326_v1.tif"


PERCENTILES_LIST = [1, 5, 10, 25, 75, 90, 95, 99]

"""
- mean
- standard deviation
- max
- min
- 99th, 95th, 90th, 75th, 25th, 10th, 5th, 1st percentiles (would be nice if this tool let me pick the percentiles, actually now that I'm saying this I'm almost sure you did this for me before... I'm going to look and will update!)
I would like each of those to be a column, and each polygon-raster combo to be a row. And then also have two columns reporting what those are:
- polygon id
- raster name
"""

OUTPUT_DIR = "./results"
CLIPPED_DIR = os.path.join(OUTPUT_DIR, "clipped")
for dirpath in [OUTPUT_DIR, CLIPPED_DIR]:
    os.makedirs(dirpath, exist_ok=True)


def extract_raster_array_by_feature(
    raster_path_band, vector_path, fid, output_tif=None
):
    """Extracts and masks a raster subset array using a specified vector feature.

    The function extracts the minimal bounding window from the raster dataset
    that fully contains the specified feature geometry from a vector dataset.
    It projects the vector geometry to match the raster coordinate system if
    needed, reads only the overlapping portion of the raster as a NumPy array,
    and applies a mask so that raster values outside the feature geometry
    are set to NoData. Optionally, the resulting array can be exported as a
    GeoTIFF file for debugging or visualization purposes. Note, this code is
    optimized so the minimal data access is done on raster_path to maximize
    disk throughput.

    Args:
        raster_path_band (tuple): (str, integer) path band tuple to analyze.
        vector_path (str): Path to the vector dataset (e.g., GeoPackage, Shapefile).
        fid (int): Feature ID of the target geometry in the first layer of the vector dataset.
        output_tif (str, optional): Path to output a GeoTIFF file of the clipped raster
            subset for verification. Defaults to None.

    Returns:
        numpy.ndarray: Masked array representing the raster data clipped by
            the specified vector feature geometry. Pixels outside the feature
            geometry are assigned the raster's NoData value.
    """
    try:
        raster_path, raster_band_id = raster_path_band
        raster = gdal.Open(raster_path, gdal.GA_ReadOnly)
        gt = raster.GetGeoTransform()
        inv_gt = gdal.InvGeoTransform(gt)

        raster_projection = osr.SpatialReference()
        raster_projection.ImportFromWkt(raster.GetProjection())
        raster_projection.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

        vector = ogr.Open(vector_path)
        layer = vector.GetLayer(0)
        feature = layer.GetFeature(fid)
        geometry = feature.GetGeometryRef().Clone()
        if geometry.GetSpatialReference() is None:
            geometry.AssignSpatialReference(layer.GetSpatialRef())
        if not geometry.GetSpatialReference().IsSame(raster_projection):
            geometry.TransformTo(raster_projection)

        minx, maxx, miny, maxy = geometry.GetEnvelope()
        LOGGER.debug(f"gt: {inv_gt}, minx: {minx}, maxy: {maxy}")
        payload = gdal.ApplyGeoTransform(inv_gt, minx, maxy)
        LOGGER.debug(f"payload: {payload}")
        px_min, py_max = gdal.ApplyGeoTransform(inv_gt, minx, maxy)
        px_max, py_min = gdal.ApplyGeoTransform(inv_gt, maxx, miny)

        xoff = int(np.floor(min(px_min, px_max)))
        yoff = int(np.floor(min(py_min, py_max)))
        xsize = int(np.ceil(abs(px_max - px_min)))
        ysize = int(np.ceil(abs(py_max - py_min)))
        xoff = max(0, xoff)
        yoff = max(0, yoff)
        xsize = min(xsize, raster.RasterXSize - xoff)
        ysize = min(ysize, raster.RasterYSize - yoff)

        LOGGER.debug(f"{xoff}, {yoff}, {xsize}, {ysize}")

        band = raster.GetRasterBand(raster_band_id)
        arr = band.ReadAsArray(xoff, yoff, xsize, ysize)
        nodata = (
            band.GetNoDataValue()
            if band.GetNoDataValue() is not None
            else -9999
        )

        subset_gt = (
            gt[0] + xoff * gt[1],
            gt[1],
            0.0,
            gt[3] + yoff * gt[5],
            0.0,
            gt[5],
        )

        mem_drv = ogr.GetDriverByName("Memory")
        mem_ds = mem_drv.CreateDataSource("")
        mem_layer = mem_ds.CreateLayer(
            "clip", raster_projection, geom_type=geometry.GetGeometryType()
        )
        mem_feature = ogr.Feature(mem_layer.GetLayerDefn())
        mem_feature.SetGeometry(geometry)
        mem_layer.CreateFeature(mem_feature)

        mask_ds = gdal.GetDriverByName("MEM").Create(
            "", xsize, ysize, 1, gdal.GDT_Byte
        )
        mask_ds.SetGeoTransform(subset_gt)
        mask_ds.SetProjection(raster_projection.ExportToWkt())
        gdal.RasterizeLayer(mask_ds, [1], mem_layer, burn_values=[1])
        mask = mask_ds.ReadAsArray().astype(bool)

        arr_masked = np.where(mask, arr, nodata)

        # get equal area hectars
        equal_area_sr = osr.SpatialReference()
        equal_area_sr.ImportFromEPSG(6933)  # World Equal Area CRS
        geom_equal_area = geometry.Clone()
        geom_equal_area.TransformTo(equal_area_sr)
        feature_area_m2 = geom_equal_area.GetArea()
        feature_area_ha = feature_area_m2 / 10000.0

        if raster_projection.IsGeographic():
            # raster is geographic (degrees), so approximate pixel area via reprojection
            # LOGGER.info("it is in degrees so approximating via reprojection")
            ul_x, ul_y = subset_gt[0], subset_gt[3]
            lr_x = subset_gt[0] + xsize * subset_gt[1]
            lr_y = subset_gt[3] + ysize * subset_gt[5]

            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(ul_x, ul_y)
            ring.AddPoint(lr_x, ul_y)
            ring.AddPoint(lr_x, lr_y)
            ring.AddPoint(ul_x, lr_y)
            ring.AddPoint(ul_x, ul_y)
            bbox_geom = ogr.Geometry(ogr.wkbPolygon)
            bbox_geom.AddGeometry(ring)
            bbox_geom.AssignSpatialReference(raster_projection)
            bbox_geom.TransformTo(equal_area_sr)
            bbox_area_m2 = bbox_geom.GetArea()
            pixel_area_m2 = bbox_area_m2 / (xsize * ysize)
        else:
            # raster projection in linear units (meters)
            pixel_area_m2 = abs(gt[1] * gt[5])

        # Valid Pixel Area (ha)
        valid_pixel_count = np.count_nonzero((arr != nodata) & mask)
        valid_pixel_area_ha = (valid_pixel_count * pixel_area_m2) / 10000.0

        if output_tif:
            drv = gdal.GetDriverByName("GTiff")
            out_ds = drv.Create(
                output_tif,
                xsize,
                ysize,
                1,
                band.DataType,
                options=["COMPRESS=LZW"],
            )
            out_ds.SetGeoTransform(subset_gt)
            out_ds.SetProjection(raster_projection.ExportToWkt())
            out_band = out_ds.GetRasterBand(1)
            out_band.WriteArray(arr_masked)
            out_band.SetNoDataValue(nodata)
            out_band.FlushCache()
            out_ds = None

        running_stats = {
            # "array": arr_masked,
            "nodata": nodata,
            "feature_area_ha": feature_area_ha,
            "valid_pixel_area_ha": valid_pixel_area_ha,
        }
        running_stats["stats"] = calculate_summary_stats(arr_masked, nodata)

        raster = None
        vector = None
        mask_ds = None
        mem_ds = None
        arr_masked = None
        return running_stats
    except Exception as e:
        return f"Failure on {raster_path}: {e}: \n{traceback.format_exc()}"


def calculate_summary_stats(array, nodata):
    """Calculate the min/max/std/mean/percentile over the non-nodata pixel in raster."""
    array = array[(array != nodata) & ~np.isnan(array)]

    if array.size == 0:
        # guard against a nodata array
        array = np.array([0.0])
    stats = {
        "min": np.min(array),
        "max": np.max(array),
        "std": np.std(array),
        "mean": np.mean(array),
    }

    percentile_dict = {
        f"p{percentile}": value
        for percentile, value in zip(
            PERCENTILES_LIST, np.percentile(array, PERCENTILES_LIST)
        )
    }
    stats.update(percentile_dict)
    return stats


def main():
    """Entry point."""
    print(os.cpu_count())
    task_graph = taskgraph.TaskGraph(
        OUTPUT_DIR, os.cpu_count(), reporting_interval=15.0
    )

    vector = gdal.OpenEx(AOI_PATH)
    layer = vector.GetLayer()
    fid_name_set = set(
        [
            (feature.GetFID(), feature.GetField(AOI_NAME_KEY))
            for feature in layer
        ]
    )
    layer.ResetReading()
    LOGGER.debug(fid_name_set)

    task_list = []
    for raster_path in RASTER_PATHS_TO_SUMMARIZE:
        raster = gdal.OpenEx(raster_path)
        raster_band_description_list = [
            raster.GetRasterBand(band_index).GetDescription()
            or f"band_{band_index}"
            for band_index in range(1, raster.RasterCount + 1)
        ]
        raster = None
        for band_index, band_description in enumerate(
            raster_band_description_list
        ):
            for fid, name in fid_name_set:
                LOGGER.debug(raster_path)
                task = task_graph.add_task(
                    func=extract_raster_array_by_feature,
                    args=((raster_path, band_index + 1), AOI_PATH, fid),
                    store_result=True,
                    task_name=f"stats for {raster_path}/{band_description}",
                )
                task_list.append((raster_path, band_description, name, task))

    stats_list = []
    error_list = []
    for raster_path, band_description, name, task in task_list:
        payload = task.get()
        if isinstance(payload, str):
            error_list.append(payload)
            continue
        stat_dict = payload
        raster_name = os.path.splitext(os.path.basename(raster_path))[0]
        match = re.search(
            r"_(\d{4})(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_", raster_name
        )
        year = "unknown"
        period = "unknown"
        if match:
            year_start, month_start, day_start, year_end, month_end, day_end = (
                match.groups()
            )
            year = year_start
            period = f"{month_start}-{day_start} -- {month_end}-{day_end}"
        summary_stats = {
            "feature name": name,
            "raster name": raster_name,
            "band description": band_description,
            "year": year,
            "period": period,
            "feature_area_ha": stat_dict["feature_area_ha"],
            "valid_pixel_area_ha": stat_dict["valid_pixel_area_ha"],
        }
        summary_stats.update(stat_dict["stats"])
        stats_list.append(summary_stats)
        LOGGER.info(f"done with {raster_path}")

    stats_df = pd.DataFrame(stats_list)
    timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    output_filename = f"results_{timestamp}.csv"
    stats_df.to_csv(output_filename, index=False)
    task_graph.close()
    task_graph.join()
    print(f"all done -- results in {output_filename}!")
    print("these are the errors:\n" + "\n".join(error_list))


if __name__ == "__main__":
    main()
