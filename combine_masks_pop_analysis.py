"""Combine masks to get population counts."""

from collections import defaultdict
from datetime import datetime
import glob
import logging
import os
import sys

from ecoshard import taskgraph
from osgeo import gdal
from pyproj import CRS
from rasterio.warp import calculate_default_transform, reproject, Resampling
from shapely.geometry import box
from shapely.geometry import shape
import fiona
import geopandas as gpd
import numpy as np
import pandas as pd
import psutil
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
logging.getLogger("ecoshard.taskgraph").setLevel(logging.INFO)
logging.getLogger("rasterio").setLevel(logging.WARNING)
logging.getLogger("fiona").setLevel(logging.WARNING)


MASK_DIRS = [
    "workspace_downstream_es_analysis",
    "workspace_dist_to_hab_with_friction",
]
SUFFIXES = {
    "_aoi_rasterized_60min.tif",
    "_aoi_ds_coverage_1000m.tif",
}

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
complete_sets = {
    prefix: paths
    for prefix, paths in matched_rasters.items()
    if SUFFIXES.issubset(paths)
}

print(complete_sets)


def main():
    # iterate over complete sets and clip/resize them against the population raster so they all fit nicely
    pass


if __name__ == "__main__":
    main()
