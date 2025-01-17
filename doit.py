import logging
import sys
from ecoshard import taskgraph
import os
import geopandas as gpd
from ecoshard import geoprocessing

logging.basicConfig(
    level=logging.DEBUG,
    stream=sys.stdout,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(pathname)s.%(funcName)s:%(lineno)d] %(message)s'))
LOGGER = logging.getLogger(__name__)
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
logging.getLogger('PIL').setLevel(logging.ERROR)
logging.getLogger('ecoshard.taskgraph').setLevel(logging.INFO)
logging.getLogger('fiona').setLevel(logging.WARN)


BASE_RASTER_LOOKUP = {
    'sed_export_marine': r"D:\repositories\roadmap2030\data\ndv_0.0_sed_export_marineESA_2020-1992_change_md5_0ab0cf.tif",
    'cv_habitat_value_mar': r"D:\repositories\roadmap2030\data\ndv_0.0_cv_habitat_value_marESA2020-1992_change_md5_1643a7.tif",
    'n_export_marine': r"D:\repositories\roadmap2030\data\ndv_0.0_n_export_marineESA_2020-1992_change_val_md5_18a2b3.tif",
    'realized_pollination_on_ag_mar': r"D:\repositories\roadmap2030\data\ndv_0.0_realized_pollination_on_ag_marESA_2020-1992_fullchange_md5_8e63e2.tif",
    'sed_deposition_marine': r"D:\repositories\roadmap2030\data\ndv_0.0_sed_deposition_marineESA_2020-1992_change_md5_d23c49.tif",
}


VECTOR_PATH = r"D:\repositories\roadmap2030\data\drive-download-20250117T210959Z-001\ES_combined.shp"

OUTPUT_DIR = './results'
CLIPPED_DIR = os.path.join(OUTPUT_DIR, 'clipped')
for dirpath in [OUTPUT_DIR, CLIPPED_DIR]:
    os.makedirs(dirpath, exist_ok=True)


def create_subset(gdf, name, target_vector_path):
    subset_gdf = gdf[gdf["Name"] == name]
    subset_gdf.to_file(target_vector_path, driver="GPKG")


def clip_raster(base_raster_path, summary_vector_path, temp_clip_path):
    target_pixel_size = geoprocessing.get_raster_info(base_raster_path)['pixel_size']
    vector_bb = geoprocessing.get_vector_info(summary_vector_path)['bounding_box']
    geoprocessing.warp_raster(
        base_raster_path, target_pixel_size, temp_clip_path,
        'near', target_bb=vector_bb, vector_mask_options={
            'mask_vector_path': summary_vector_path,
            'all_touched': True})


def main():
    """Entry point."""
    task_graph = taskgraph.TaskGraph(OUTPUT_DIR, os.cpu_count(), 15.0)
    wkt_projection = geoprocessing.get_raster_info(next(iter(BASE_RASTER_LOOKUP.values())))['projection_wkt']
    gdf = gpd.read_file(VECTOR_PATH)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    unique_names = gdf["Name"].unique()
    gdf = gdf.to_crs(wkt_projection)
    for place_name in unique_names:
        if not isinstance(place_name, str) or place_name.strip() == '':
            continue
        subset_vector_path = os.path.join(OUTPUT_DIR, f"{place_name}.gpkg")
        subset_task = task_graph.add_task(
            func=create_subset,
            args=(gdf, place_name, subset_vector_path),
            ignore_path_list=[subset_vector_path],
            target_path_list=[subset_vector_path],
            task_name=f'extract {place_name}')

        for raster_basename, raster_path in BASE_RASTER_LOOKUP.items():
            clipped_raster_path = os.path.join(
                CLIPPED_DIR, f'{place_name}_{raster_basename}.tif')
            clipped_task = task_graph.add_task(
                func=clip_raster,
                args=(raster_path, subset_vector_path, clipped_raster_path),
                ignore_path_list=[subset_vector_path],
                dependent_task_list=[subset_task],
                target_path_list=[clipped_raster_path],
                task_name=f'clipping {raster_path} to {subset_vector_path}')

    task_graph.join()
    print('all done!')


if __name__ == '__main__':
    main()
