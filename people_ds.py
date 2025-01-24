import datetime
import numpy
import csv
import collections
from osgeo import gdal
import logging
import sys
from ecoshard import taskgraph
import os
from ecoshard import geoprocessing
from ecoshard.geoprocessing import routing

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

POPULATION_RASTERS = [
    './data/lspop2019_compressed_md5_d0bf03bd0a2378196327bbe6e898b70c.tif',
    './data/floodplains_masked_pop_30s_md5_c027686bb9a9a36bdababbe8af35d696.tif',]

ANALYSIS_TUPLES = {
    'non-arpa': (
        './data/non-arpa-projected-in-m.gpkg',
        './data/hybas_sa_lev05_intersect_non-arpa.gpkg',
        './data/dem_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif'),
    'Arctic_si': (
        './data/Arctic.gpkg',
        './data/hybas_si_lev05_intersect_Arctic_si.gpkg',
        './data/dem_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_si_lev05_intersect_Arctic_sii.tif'),
    'Arctic_ar': (
        './data/Arctic.gpkg',
        './data/hybas_ar_lev05_intersect_Arctic_ar.gpkg',
        './data/dem_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_ar_lev05_intersect_Arctic_ar.tif'),
    'arpa': (
        './data/arpa-projected-in-m.gpkg',
        './data/hybas_sa_lev05_intersect_arpa.gpkg',
        './data/dem_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_arpa.tif'),
    'Colombia': (
        './data/Colombia.gpkg',
        './data/hybas_sa_lev05_intersect_Colombia.gpkg',
        './data/dem_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_Colombia.tif'),
    'NGP': (
        './data/NGP.gpkg',
        './data/hybas_na_lev05_intersect_NGP.gpkg',
        './data/dem_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_na_lev05_intersect_NGP.tif'),
    'Peru': (
        './data/Peru.gpkg',
        './data/hybas_sa_lev05_intersect_Peru.gpkg',
        './data/dem_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_Peru.tif'),
    'RGRB': (
        './data/RGRB.gpkg',
        './data/hybas_na_lev05_intersect_RGRB.gpkg',
        './data/dem_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_na_lev05_intersect_RGRB.tif'),
    'Tapajos': (
        './data/Tapajos.gpkg',
        './data/hybas_sa_lev05_intersect_Tapajos.gpkg',
        './data/dem_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_Tapajos.tif'),

}

OUTPUT_DIR = './results'
CLIPPED_DIR = os.path.join(OUTPUT_DIR, 'clipped')
for dirpath in [OUTPUT_DIR, CLIPPED_DIR]:
    os.makedirs(dirpath, exist_ok=True)


def calc_flow_dir(analysis_id, base_dem_raster_path, aoi_vector_path, target_flow_dir_path):
    local_workspace_dir = os.path.dirname(target_flow_dir_path)
    os.makedirs(local_workspace_dir, exist_ok=True)
    clipped_dem_path = os.path.join(local_workspace_dir, f'{analysis_id}_dem.tif')
    dem_info = geoprocessing.get_raster_info(base_dem_raster_path)
    aoi_ds = gdal.OpenEx(aoi_vector_path, gdal.OF_VECTOR | gdal.GA_ReadOnly)
    aoi_layer = aoi_ds.GetLayer()
    minx, maxx, miny, maxy = aoi_layer.GetExtent()
    aoi_bb = [minx, miny, maxx, maxy]
    bounding_box = geoprocessing.merge_bounding_box_list(
        [dem_info['bounding_box'], aoi_bb], 'intersection')
    geoprocessing.warp_raster(
        base_dem_raster_path, dem_info['pixel_size'], clipped_dem_path,
        'nearest', target_bb=bounding_box, vector_mask_options={
            'mask_vector_path': aoi_vector_path,
            'all_touched': True,
            'target_mask_value': 0
        })
    r = gdal.OpenEx(clipped_dem_path, gdal.OF_RASTER | gdal.OF_UPDATE)
    b = r.GetRasterBand(1)
    b.SetNoDataValue(0)
    b = None
    r = None

    filled_dem_path = os.path.join(local_workspace_dir, f'{analysis_id}_dem_filled.tif')
    routing.fill_pits(
        (clipped_dem_path, 1), filled_dem_path,
        working_dir=local_workspace_dir)

    routing.flow_dir_mfd(
        (filled_dem_path, 1), target_flow_dir_path,
        working_dir=local_workspace_dir)


def rasterize(aoi_vector_path, dem_raster_path, aoi_raster_mask_path):
    geoprocessing.new_raster_from_base(
        dem_raster_path, aoi_raster_mask_path, datatype=gdal.GDT_Byte, band_nodata_list=[-1])
    geoprocessing.rasterize(
        aoi_vector_path, aoi_raster_mask_path, burn_values=[1])


def mask_by_nonzero_and_sum(analysis_id, base_raster_path, mask_raster_path, target_masked_path):
    base_raster_info = geoprocessing.get_raster_info(base_raster_path)
    mask_raster_info = geoprocessing.get_raster_info(mask_raster_path)
    nodata = base_raster_info['nodata'][0]
    def _mask_by_nonzero(base_array, mask_array):
        result = base_array.copy()
        result[mask_array == 0] = nodata
        return result

    working_dir_path = os.path.dirname(target_masked_path)
    aligned_raster_path_list = [
        os.path.join(working_dir_path, f'%s_{analysis_id}_aligned%s' % os.path.splitext(os.path.basename(path)))
        for path in [base_raster_path, mask_raster_path]]
    geoprocessing.align_and_resize_raster_stack(
        [base_raster_path, mask_raster_path],
        aligned_raster_path_list, ['near']*2,
        base_raster_info['pixel_size'], mask_raster_info['bounding_box'])

    geoprocessing.raster_calculator(
        [(path, 1) for path in aligned_raster_path_list],
        _mask_by_nonzero, target_masked_path,
        base_raster_info['datatype'], nodata,
        allow_different_blocksize=True,
        skip_sparse=True)

    array = gdal.OpenEx(target_masked_path).ReadAsArray()
    array = array[array != nodata]
    return numpy.sum(array)


def main():
    """Entry point."""
    task_graph = taskgraph.TaskGraph(OUTPUT_DIR, os.cpu_count(), reporting_interval=10.0)
    task_graph.join()

    result = collections.defaultdict(lambda: collections.defaultdict())

    for analysis_id, (aoi_vector_path, subwatershed_vector_path, dem_raster_path) in ANALYSIS_TUPLES.items():
        local_workspace_dir = os.path.join(OUTPUT_DIR, analysis_id)
        os.makedirs(local_workspace_dir, exist_ok=True)
        aoi_raster_mask_path = os.path.join(
            local_workspace_dir, f'{analysis_id}_aoi_mask.tif')
        target_projection_wkt = geoprocessing.get_raster_info(dem_raster_path)['projection_wkt']
        reprojected_aoi_vector_path = '%s_projected.gpkg' % os.path.splitext(aoi_raster_mask_path)[0]
        reproject_task = task_graph.add_task(
            func=geoprocessing.reproject_vector,
            args=(
                aoi_vector_path, target_projection_wkt, reprojected_aoi_vector_path),
            target_path_list=[reprojected_aoi_vector_path],
            task_name=f'reproject {analysis_id}')
        flow_dir_path = os.path.join(local_workspace_dir, f'{analysis_id}_mfd_flow_dir.tif')
        flow_dir_task = task_graph.add_task(
            func=calc_flow_dir,
            args=(analysis_id, dem_raster_path, subwatershed_vector_path, flow_dir_path),
            target_path_list=[flow_dir_path],
            task_name=f'calculate flow dir for {analysis_id}')

        rasterize_task = task_graph.add_task(
            func=rasterize,
            args=(reprojected_aoi_vector_path, dem_raster_path, aoi_raster_mask_path),
            dependent_task_list=[reproject_task],
            target_path_list=[aoi_raster_mask_path],
            task_name=f'{analysis_id} raster mask')

        aoi_downstream_flow_mask_path = os.path.join(
            local_workspace_dir, f'{analysis_id}_aoi_ds_coverage.tif')
        flow_accum_task = task_graph.add_task(
            func=routing.flow_accumulation_mfd,
            args=(
                (flow_dir_path, 1), aoi_downstream_flow_mask_path),
            kwargs={'weight_raster_path_band': (aoi_raster_mask_path, 1)},
            dependent_task_list=[flow_dir_task, rasterize_task],
            target_path_list=[aoi_downstream_flow_mask_path],
            task_name=f'flow accum for {analysis_id}')
        for population_raster_path in POPULATION_RASTERS:
            pop_basename = os.path.splitext(os.path.basename(population_raster_path))[0]
            print(f'processing {pop_basename}')

            masked_population_raster_path = os.path.join(
                OUTPUT_DIR, f'{analysis_id}_{os.path.basename(population_raster_path)}')
            mask_by_nonzero_task = task_graph.add_task(
                func=mask_by_nonzero_and_sum,
                args=(
                    f'{analysis_id}_{pop_basename}', population_raster_path, aoi_downstream_flow_mask_path, masked_population_raster_path),
                dependent_task_list=[flow_accum_task],
                target_path_list=[masked_population_raster_path],
                store_result=True,
                task_name=f'{analysis_id}_population mask')

            result[analysis_id][pop_basename] = mask_by_nonzero_task

    for analysis_id in result:
        for pop_basename in result[analysis_id]:
            result[analysis_id][pop_basename] = result[analysis_id][pop_basename].get()

    timestamp = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    output_filename = f'pop_results_{timestamp}.csv'
    with open(output_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['analysis_id', 'pop_basename', 'value'])
        for analysis_id in result:
            for pop_basename, val in result[analysis_id].items():
                writer.writerow([analysis_id, pop_basename, val])
    task_graph.join()
    task_graph.close()
    LOGGER.info(f'all done results at {output_filename}')


if __name__ == '__main__':
    main()
