"""Distance to habitat with a friction layer."""

from datetime import datetime
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
import numpy as np
import pandas as pd
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
BUFFER_DISTANCE_M = 104 * 8 * 1000  # driving 8 hours straight


def get_utm_crs(geometry, original_crs):
    # Ensure centroid is calculated in lat/lon
    if original_crs.is_geographic:
        centroid = geometry.centroid
    else:
        # temporarily project to lat/lon for correct centroid calculation
        geometry_latlon = (
            gpd.GeoSeries([geometry], crs=original_crs).to_crs("EPSG:4326").iloc[0]
        )
        centroid = geometry_latlon.centroid

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
    results_dict = {}
    for aoi_basename, aoi_path in ANALYSIS_AOIS.items():
        for max_hours in range(1, 9):
            max_time_mins = max_hours * 60
            local_workspace_dir = os.path.join(WORKSPACE_DIR, aoi_basename)
            os.makedirs(local_workspace_dir, exist_ok=True)
            LOGGER.info(f"processing {aoi_basename}")
            gdf = gpd.read_file(aoi_path)
            projected_gdf = gdf.to_crs(get_utm_crs(gdf.union_all(), gdf.crs))
            bbox = projected_gdf.total_bounds
            buffer_distance_m = (
                max_hours * 104 * 1000
            )  # drive 65mph for that many hours
            buffered_bbox = box(
                bbox[0] - buffer_distance_m,
                bbox[1] - buffer_distance_m,
                bbox[2] + buffer_distance_m,
                bbox[3] + buffer_distance_m,
            )
            bbox_gdf = gpd.GeoDataFrame(
                {"geometry": [buffered_bbox]}, crs=projected_gdf.crs
            )

            target_pop_clipped_raster_path = os.path.join(
                local_workspace_dir,
                f"{aoi_basename}_population_clipped_{max_time_mins}min.tif",
            )
            target_friction_clipped_raster_path = os.path.join(
                local_workspace_dir,
                f"{aoi_basename}_friction_surface_clipped.tif",
            )
            target_aoi_raster_path = os.path.join(
                local_workspace_dir,
                f"{aoi_basename}_aoi_rasterized_{max_time_mins}min.tif",
            )
            clip_and_reproject_raster(
                POPULATION_RASTER_PATH,
                bbox_gdf,
                projected_gdf.crs,
                target_pop_clipped_raster_path,
            )

            with rasterio.open(target_pop_clipped_raster_path) as pop_ref:
                ref_meta = pop_ref.meta.copy()
                pop_array = pop_ref.read(1).astype(np.int64)

            clip_and_reproject_raster(
                FRICTION_SURFACE_RASTER_PATH,
                bbox_gdf,
                projected_gdf.crs,
                target_friction_clipped_raster_path,
                reference_meta=ref_meta,
            )

            mask_array = rasterio.features.rasterize(
                ((geom, 1) for geom in projected_gdf.geometry),
                out_shape=(ref_meta["height"], ref_meta["width"]),
                transform=ref_meta["transform"],
                fill=0,
                dtype=rasterio.uint8,
            ).astype(np.int8)

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
                dst.write(mask_array, 1)

            with rasterio.open(target_friction_clipped_raster_path) as friction_ds:
                friction_array = friction_ds.read(1)
                transform = friction_ds.transform
                cell_length_m = (
                    transform.a
                )  # this is the 'a' in the affine transformation, so cell x length
                n_rows, n_cols = friction_array.shape

            LOGGER.info(f"about to calculate {n_cols} X {n_rows} raster")
            target_max_reach_raster_path = (
                f"{aoi_basename}_max_reach_{max_time_mins}min.tif"
            )
            travel_reach_array = shortest_distances.find_mask_reach(
                friction_array,
                mask_array,
                cell_length_m,
                n_cols,
                n_rows,
                max_time_mins,
            )
            with rasterio.open(
                target_max_reach_raster_path, "w", **aoi_meta
            ) as max_reach:
                max_reach.write(travel_reach_array, 1)
            in_range_pop_count = np.sum(
                pop_array[(travel_reach_array > 0) & (pop_array > 0)]
            )

            LOGGER.debug(
                f"{aoi_basename} {max_time_mins//60}hours {in_range_pop_count}"
            )
            if aoi_basename not in results_dict:
                results_dict[aoi_basename] = {}
            results_dict[aoi_basename][max_hours] = in_range_pop_count

    df = pd.DataFrame(results_dict).T.sort_index(axis=0)
    df = df.sort_index(axis=1)

    # Rename columns clearly
    df.columns = [f"{hour}_hours" for hour in df.columns]

    # Save DataFrame to CSV with current timestamp in filename
    timestamp = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    csv_filename = f"people_within_travel_time_{timestamp}.csv"
    df.to_csv(csv_filename)

    LOGGER.info(f"Saved results to {csv_filename}")


if __name__ == "__main__":
    main()
