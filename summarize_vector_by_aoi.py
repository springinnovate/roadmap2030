"""# noqa: D205, D210, D415
This script calculates summary statistics for vector data features within
AOIs overlaid with separate vector layers containing attributes. For each
pair of AOI x vector data layer the script computes intersections and
generates statistics (min, max, mean (weighted by area), percentiles
(weighted by area) the intersected features. The output provides a
cross-tabulated summary of how data attributes vary across AOIs.
Use py311 conda environment: mamba activate py311
"""

from datetime import datetime
from pathlib import Path
import logging
import os
import sys

from shapely.geometry import box
import geopandas as gdf
import numpy as np
import pandas as pd


logging.basicConfig(
    level=logging.INFO,
    stream=sys.stdout,
    format=(
        "%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s "
        "[%(filename)s.%(funcName)s:%(lineno)d] %(message)s"
    ),
)
LOGGER = logging.getLogger(__name__)
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
logging.getLogger("PIL").setLevel(logging.ERROR)
logging.getLogger("ecoshard.taskgraph").setLevel(logging.INFO)
logging.getLogger("fiona").setLevel(logging.WARN)


AOI_PATHS = [
    r"D:\repositories\roadmap2030\data\WWF-US_Pilot\aois\small_set\NGP.gpkg",
    r"D:\repositories\roadmap2030\data\WWF-US_Pilot\aois\small_set\RGRB.gpkg",
]

VECTOR_DATA_FIELD_PATHS = [  # this expects a tuple: (path, field)
    (
        r"D:\repositories\data_platform\Climate\Resilience\TNC_Freshwater_Resilience_Data_Nov2023\Scored_Units.gpkg",
        "Resil_Z",
    ),
]

# percentiles are in % units
PERCENTILES_LIST = [1, 5, 10, 25, 50, 75, 90, 95, 99]

VALUE_THRESHOLDS_FOR_AREA_SUM = [
    0.15,
    0.3,
    0.96,
]


def _validate_data():
    errors = []
    for path in AOI_PATHS:
        if not os.path.exists(path):
            errors.append(f"Missing AOI file: {path}")

    for path, field in VECTOR_DATA_FIELD_PATHS:
        if not os.path.exists(path):
            errors.append(f"Missing vector file: {path}")
            continue
        try:
            vector_gdf = gdf.read_file(path, rows=1)
        except Exception as e:
            errors.append(f"Could not read vector file {path}: {e}")
            continue
        if field not in vector_gdf.columns:
            errors.append(f'Missing field "{field}" in {path}')

    if errors:
        message = "Validation failed with the following issues:\n" + "\n".join(errors)
        raise RuntimeError(message)


def summarize_aoi_vector(
    aoi_gdf: gdf.GeoDataFrame, vec_gdf: gdf.GeoDataFrame, vector_field: str
):
    """Given the prepped geodataframes, do the summary."""
    if aoi_gdf.crs is None:
        raise ValueError("AOI has no CRS")
    if vec_gdf.crs is None:
        raise ValueError("Vector has no CRS")

    # if the vector is in a different crs, use the AOI's as the governing one
    if vec_gdf.crs != aoi_gdf.crs:
        LOGGER.debug(f"reprojecting vector to aoi's crs: {aoi_gdf.crs}")
        vec_gdf = vec_gdf.to_crs(aoi_gdf.crs)

    # compute AOI bbox for filtered read
    aoi_gdf = gdf.GeoDataFrame(geometry=[aoi_gdf.union_all()], crs=aoi_gdf.crs)
    minx, miny, maxx, maxy = aoi_gdf.total_bounds
    bbox = (minx, miny, maxx, maxy)

    # quick bbox prefilter to toss polys outside of the aoi
    idx = vec_gdf.sindex
    cand_ix = list(idx.query(box(*bbox)))
    vec_gdf = vec_gdf.iloc[cand_ix]

    # geometries are always invalid...
    if not aoi_gdf.geometry.is_valid.all():
        LOGGER.debug("fixing invalid aoi geometry")
        aoi_gdf["geometry"] = aoi_gdf.buffer(0)
    if not vec_gdf.geometry.is_valid.all():
        LOGGER.debug("fixing invalid vector geometry")
        vec_gdf["geometry"] = vec_gdf.buffer(0)

    vec_gdf = gdf.clip(vec_gdf, aoi_gdf)
    # drop empty/zero-area geometries if any snuck through
    vec_gdf = vec_gdf[~vec_gdf.geometry.is_empty & vec_gdf.geometry.notna()]

    # calculate area field
    vec_gdf["__area_m2"] = vec_gdf.geometry.area
    vec_gdf = vec_gdf[vec_gdf["__area_m2"] > 0]

    values = vec_gdf[vector_field].to_numpy(dtype="float64")
    areas = vec_gdf["__area_m2"].to_numpy(dtype="float64")

    area_threshold_sums = {}
    for t in VALUE_THRESHOLDS_FOR_AREA_SUM:
        mask = values >= t
        area_threshold_sums[f"area_sum_ge_{t}"] = float(areas[mask].sum())

    def _weighted_percentiles(x, w, p):
        order = np.argsort(x)
        x_sorted = x[order]
        w_sorted = w[order]
        cw = np.cumsum(w_sorted)
        total = cw[-1]
        # calculate cutoffs based on total WEIGHT
        targets = np.asarray(p, dtype="float64") / 100.0 * total
        # interpolate because we might not be exact
        return np.interp(targets, cw, x_sorted)

    unweighted_percentiles = dict(
        zip(
            [f"p{int(p)}" for p in PERCENTILES_LIST],
            np.percentile(values, PERCENTILES_LIST).astype(float),
        )
    )

    weighted_percentiles = dict(
        zip(
            [f"wp{int(p)}" for p in PERCENTILES_LIST],
            _weighted_percentiles(values, areas, PERCENTILES_LIST).astype(float),
        )
    )

    stats = (
        {
            "count": len(vec_gdf),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
            "sum": float(np.sum(values)),
            "mean": float(np.mean(values)),
            "std": (float(np.std(values, ddof=1)) if values.size > 1 else float("nan")),
            "area_sum_m2": float(areas.sum()),
            "area_weighted_mean": float(np.average(values, weights=areas)),
            "area_weighted_sum": float(np.sum(values * areas)),
        }
        | unweighted_percentiles  # noqa: W503
        | weighted_percentiles  # noqa: W503
        | area_threshold_sums  # noqa: W503
    )
    return stats


def _area_error_pct_one(vector, tol, idx=0):
    geom = vector.geometry.iloc[idx]
    simp = geom.simplify(tol, preserve_topology=True)
    if geom.is_empty or simp.is_empty:
        return 0.0
    area_err = abs(geom.area - simp.area) / max(geom.area, 1e-9) * 100
    return area_err


def choose_tolerance(vector_path, target_error_pct=0.1, guess=30):
    """Calcualte a valid tolerance that doesn't change the area too much."""
    vector = gdf.read_file(vector_path, rows=1)
    if vector.crs is None or not vector.crs.is_projected:
        raise ValueError("Use a projected CRS (meters) before simplifying.")

    # expand upward until error exceeds target
    low = 0
    tol = guess
    err = _area_error_pct_one(vector, tol)
    while err <= target_error_pct:
        low = tol
        tol *= 2
        err = _area_error_pct_one(vector, tol)

    high = tol

    # binary search between low and high
    while high - low > 1:  # stop when 1 meter apart
        mid = (low + high) // 2
        err = _area_error_pct_one(vector, mid)
        if err <= target_error_pct:
            low = mid
        else:
            high = mid

    return low


def load_and_simplify(vector_path):
    """Reads the vector path then simplifies it to a reasonable tolerance."""
    tolerance = choose_tolerance(vector_path)
    LOGGER.debug(f"reading full file at {vector_path}")
    vector = gdf.read_file(vector_path, engine="pyogrio")
    LOGGER.debug(f"simplifying to tolerance {tolerance}")
    vector["geometry"] = vector.geometry.simplify(tolerance, preserve_topology=False)
    return vector


def main():
    """Entry point."""
    LOGGER.info("validating data...")
    _validate_data()
    results = []
    for aoi_path in AOI_PATHS:
        for vector_path, vector_field in VECTOR_DATA_FIELD_PATHS:
            LOGGER.info(f"pre-processing {aoi_path}, {vector_path}, {vector_field}")
            LOGGER.debug(f"prepping {aoi_path}")
            aoi_gdf = load_and_simplify(aoi_path)
            LOGGER.debug(f"prepping {vector_path}")
            vec_gdf = load_and_simplify(vector_path)
            LOGGER.info("starting summary")
            stats = summarize_aoi_vector(aoi_gdf, vec_gdf, vector_field)
            stats["aoi"] = Path(aoi_path).stem
            stats["vector"] = Path(vector_path).stem
            stats["field"] = vector_field
            results.append(stats)
    df = pd.DataFrame(results)
    core_columns = ["aoi", "vector", "field"]
    stat_keys = [k for k in results[0].keys() if k not in core_columns]
    df = df[core_columns + stat_keys]

    timestamp = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    out_path = f"results_{timestamp}.csv"
    df.to_csv(out_path, index=False)
    LOGGER.info(f"wrote {out_path}")


if __name__ == "__main__":
    main()
