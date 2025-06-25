"""Distance to habitat with a friction layer."""

import glob
import logging
import os
import sys

from osgeo import gdal
from pyproj import CRS
from rasterio.warp import calculate_default_transform, reproject, Resampling
from shapely.geometry import box
from shapely.geometry import shape
import fiona
import geopandas as gpd
import rasterio
import rasterio.mask

import shortest_distances

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
logging.getLogger("taskgraph").setLevel(logging.DEBUG)
logging.getLogger("rasterio").setLevel(logging.WARNING)
logging.getLogger("fiona").setLevel(logging.WARNING)


AOI_DIRS = ["./data/WWF-Int_Pilot"]

ANALYSIS_AOIS = {}
BAD_AOIS = {}

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


FRICTION_SURFACE_RASTER_PATH = "./data/travel_time/friction_surface_2019_compressed_md5_1be7dd230178a5d395529be7a5e3fb0a.tif"
POPULATION_RASTER_PATH = "./data/pop_rasters/landscan-global-2023.tif"
WORKSPACE_DIR = "workspace_dist_to_hab_with_friction"

CHURN_DIR = os.path.join(WORKSPACE_DIR, "churn")
TARGET_NODATA = -1

# used to avoid computing paths where the population is too low
POPULATION_COUNT_CUTOFF = 0
MAX_TIME = 8 * 60  # 8 hours
BUFFER_DISTANCE_M = 104 * 8 * 1000  # driving 8 hours straight


def get_utm_crs(geometry):
    centroid = geometry.centroid
    utm_zone = int((centroid.x + 180) // 6) + 1
    hemisphere = "north" if centroid.y >= 0 else "south"
    return CRS.from_dict(
        {"proj": "utm", "zone": utm_zone, "south": hemisphere == "south"}
    )


def clip_and_reproject_raster(
    src_path, bbox_geom, dst_crs, dst_path, reference_meta=None
):
    with rasterio.open(src_path) as src:
        bbox_geom_src_crs = bbox_geom.to_crs(src.crs)
        geom = [bbox_geom_src_crs.geometry.iloc[0]]
        out_image, out_transform = rasterio.mask.mask(src, geom, crop=True)

        if reference_meta:
            dst_transform = reference_meta["transform"]
            width = reference_meta["width"]
            height = reference_meta["height"]
        else:
            dst_transform, width, height = calculate_default_transform(
                src.crs,
                dst_crs,
                out_image.shape[2],
                out_image.shape[1],
                *bbox_geom_src_crs.total_bounds,
            )

        out_meta = src.meta.copy()
        out_meta.update(
            {
                "driver": "GTiff",
                "crs": dst_crs,
                "transform": dst_transform,
                "width": width,
                "height": height,
            }
        )

        with rasterio.open(dst_path, "w", **out_meta) as dst:
            reproject(
                source=out_image,
                destination=rasterio.band(dst, 1),
                src_transform=out_transform,
                src_crs=src.crs,
                dst_transform=dst_transform,
                dst_crs=dst_crs,
                resampling=Resampling.nearest,
            )


def main():
    """Entry point."""
    for dir_path in [WORKSPACE_DIR, CHURN_DIR]:
        os.makedirs(dir_path, exist_ok=True)
    for aoi_basename, aoi_path in ANALYSIS_AOIS.items():
        LOGGER.debug(aoi_path)
        gdf = gpd.read_file(aoi_path)
        projected_gdf = gdf.to_crs(get_utm_crs(gdf.union_all()))
        bbox = projected_gdf.total_bounds
        buffered_bbox = box(
            bbox[0] - BUFFER_DISTANCE_M,
            bbox[1] - BUFFER_DISTANCE_M,
            bbox[2] + BUFFER_DISTANCE_M,
            bbox[3] + BUFFER_DISTANCE_M,
        )
        bbox_gdf = gpd.GeoDataFrame(
            {"geometry": [buffered_bbox]}, crs=projected_gdf.crs
        )
        LOGGER.debug(bbox_gdf)
        target_pop_clipped_raster_path = f"{aoi_basename}_population_clipped.tif"
        target_friction_clipped_raster_path = (
            f"{aoi_basename}_friction_surface_clipped.tif"
        )
        target_aoi_raster_path = f"{aoi_basename}_aoi_rasterized.tif"
        clip_and_reproject_raster(
            POPULATION_RASTER_PATH,
            bbox_gdf,
            projected_gdf.crs,
            target_pop_clipped_raster_path,
        )

        with rasterio.open(target_pop_clipped_raster_path) as pop_ref:
            ref_meta = pop_ref.meta.copy()

        clip_and_reproject_raster(
            FRICTION_SURFACE_RASTER_PATH,
            bbox_gdf,
            projected_gdf.crs,
            target_friction_clipped_raster_path,
            reference_meta=ref_meta,
        )

        rasterized_shape = rasterio.features.rasterize(
            ((geom, 1) for geom in projected_gdf.geometry),
            out_shape=(ref_meta["height"], ref_meta["width"]),
            transform=ref_meta["transform"],
            fill=0,
            dtype=rasterio.uint8,
        )

        # Write rasterized AOI
        aoi_meta = ref_meta.copy()
        aoi_meta.update(
            {
                "count": 1,
                "dtype": rasterio.uint8,
                "nodata": 0,
                "compress": "lzw",
            }
        )

        with rasterio.open(target_aoi_raster_path, "w", **aoi_meta) as dst:
            dst.write(rasterized_shape, 1)

        with rasterio.open(target_friction_clipped_raster_path) as friction_ds:
            friction_array = friction_ds.read(1)
            transform = friction_ds.transform
            cell_length_m = transform.a  # pixel size (assumes square pixels)
            n_rows, n_cols = friction_array.shape

        with rasterio.open(target_aoi_raster_path) as mask_ds:
            mask_array = mask_ds.read(1).astype("int8")

        LOGGER.info(f"about to calculate {n_cols} X {n_rows} raster")
        target_max_reach_raster_path = f"{aoi_basename}_max_reach.tif"
        array = shortest_distances.find_mask_reach(
            friction_array,
            mask_array,
            cell_length_m,
            n_cols,
            n_rows,
            MAX_TIME,
        )
        with rasterio.open(target_max_reach_raster_path, "w", **aoi_meta) as max_reach:
            max_reach.write(array, 1)

        LOGGER.info("all done")
        return


if __name__ == "__main__":
    main()
