"""Generates combined AOI vectors clipped by specified country boundaries.

This script crawls through a directory containing vector files
identifies valid vectors, then iteratively combines AOI vectors intersecting
predefined countries. The AOI vectors are clipped to match country
boundaries from a provided country vector dataset. Results are saved as
GeoPackage file for each country.
"""

import glob
import logging
import os
import sys

from shapely.geometry import shape
import fiona
import geopandas as gpd
import pandas as pd

COUNTRY_VECTOR_PATH = (
    "data/countries/countries_iso3_md5_6fb2431e911401992e6e56ddf0a9bcda.gpkg"
)
COUNTRY_NAME_FIELD_ID = "iso3"
AOI_DIR = "./data/WWF-Int_Pilot"
OUTPUT_DIR = "./data/aoi_by_country"


logging.basicConfig(
    level=logging.DEBUG,
    stream=sys.stdout,
    format=(
        "%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s"
        " [%(pathname)s.%(funcName)s:%(lineno)d] %(message)s"
    ),
)
LOGGER = logging.getLogger(__name__)
logging.getLogger("fiona").setLevel(logging.WARNING)


def crawl_for_valid_vectors(dir_path):
    valid_aois = {}
    for ext in ["shp", "gpkg"]:
        pattern = os.path.join(dir_path, "**", f"*.{ext}")
        for file_path in glob.glob(pattern, recursive=True):
            basename = os.path.splitext(os.path.basename(file_path))[0]
            try:
                with fiona.open(file_path, "r") as src:
                    if any(
                        shape(feat["geometry"]).is_valid
                        and not shape(feat["geometry"]).is_empty  # noqa: W503
                        for feat in src
                    ):
                        valid_aois[basename] = file_path
            except Exception:
                LOGGER.exception(f"issue with {file_path}")
                continue
    return valid_aois


def main():
    """Entry point."""
    base_vector_dict = crawl_for_valid_vectors(AOI_DIR)
    country_gdf = gpd.read_file(COUNTRY_VECTOR_PATH)
    country_crs = country_gdf.crs
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for country_name in country_gdf[COUNTRY_NAME_FIELD_ID]:
        country_geom = gpd.GeoDataFrame(
            country_gdf[country_gdf[COUNTRY_NAME_FIELD_ID] == country_name].geometry
        ).union_all()

        intersections = []

        for aoi_name, aoi_path in base_vector_dict.items():
            aoi_gdf = gpd.read_file(aoi_path)
            if aoi_gdf.crs != country_crs:
                aoi_gdf = aoi_gdf.to_crs(country_crs)

            intersection = aoi_gdf[aoi_gdf.geometry.intersects(country_geom)][
                ["geometry"]
            ].copy()
            if not intersection.empty:
                intersection.geometry = intersection.geometry.intersection(country_geom)
                intersections.append(intersection)

        if intersections:
            combined_gdf = gpd.GeoDataFrame(
                pd.concat(intersections, ignore_index=True), crs=country_crs
            )
            output_filename = os.path.join(
                OUTPUT_DIR, f"{country_name}_combined_aoi.gpkg"
            )
            combined_gdf.to_file(output_filename, driver="GPKG")
            print(f"Saved combined AOI for {country_name} to {output_filename}")


if __name__ == "__main__":
    main()
