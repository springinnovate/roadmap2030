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
import numpy

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
    './data/pop_rasters/landscan-global-2023.tif',]
    #'./data/pop_rasters/lspop2019_compressed_md5_d0bf03bd0a2378196327bbe6e898b70c.tif',
    #'./data/pop_rasters/floodplains_masked_pop_30s_md5_c027686bb9a9a36bdababbe8af35d696.tif',]

POP_PIXEL_SIZE = [0.008333333333333, -0.008333333333333]
km_in_deg = POP_PIXEL_SIZE[0] * 1000 / 900  # because it's 900m so converting to 1km
BUFFER_AMOUNTS_IN_PIXELS_M = [
    (int(numpy.round(x * km_in_deg / POP_PIXEL_SIZE[0])), x*1000)
    for x in range(0, 6)]  # try buffer zones of 0-5km

ANALYSIS_TUPLES = {
    '37_GEF_Peru': (
       './data/37_GEF_Peru.gpkg',
       './data/subwatersheds/hybas_sa_lev05_Amazon.gpkg',
       './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif_merged_compressed.tif'),
    '49_GEF_Guyana_KPMA_NRW': (
       './data/49_GEF_Guyana_KPMA_NRW.gpkg',
       './data/subwatersheds/hybas_sa_lev05_Amazon.gpkg',
       './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif_merged_compressed.tif'),
    '105_Arpa_nonoverlapping_clipped': (
       './data/105_Arpa_nonoverlapping_clipped.gpkg',
       './data/subwatersheds/hybas_sa_lev05_Amazon.gpkg',
       './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif_merged_compressed.tif'),
    '106_HECO': (
       './data/106_HECO.gpkg',
       './data/subwatersheds/hybas_sa_lev05_Amazon.gpkg',
       './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif_merged_compressed.tif'),
    '110_PatriomonioPeru': (
       './data/110_PatriomonioPeru.gpkg',
       './data/subwatersheds/hybas_sa_lev05_Amazon.gpkg',
       './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif_merged_compressed.tif'),
    '118_NBSOP2': (
       './data/118_NBSOP2.gpkg',
       './data/subwatersheds/hybas_sa_lev05_Amazon.gpkg',
       './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif_merged_compressed.tif'),
    '118_NBSOP3': (
       './data/118_NBSOP3.gpkg',
       './data/subwatersheds/hybas_sa_lev05_Amazon.gpkg',
       './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif_merged_compressed.tif'),
    '118_NBSOP4': (
        './data/118_NBSOP4.gpkg',
        './data/subwatersheds/hybas_sa_lev05_Amazon.gpkg',
        './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif_merged_compressed.tif'),
    '120_Tapajos': (
       './data/120_Tapajos.gpkg',
       './data/subwatersheds/hybas_sa_lev05_Amazon.gpkg',
       './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif_merged_compressed.tif'),
    '133_Sall': (
       './data/133_Sall.gpkg',
       './data/subwatersheds/hybas_sa_lev05_Amazon.gpkg',
       './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif_merged_compressed.tif'),
     'Amazon_aoi': (
        './data/Amazon_aoi.gpkg',
        './data/subwatersheds/hybas_sa_lev05_Amazon.gpkg',
        './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif_merged_compressed.tif'),
     'WWF_Amazon_sites': (
        './data/000_All_WWF_Amazon_sites.gpkg',
        './data/subwatersheds/hybas_sa_lev05_Amazon.gpkg',
        './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif_merged_compressed.tif'),
    # 'arpa_split': (
    #    './data/aois/Arpa_nonoverlapping-in-m.gpkg',
    #    './data/subwatersheds/hybas_sa_lev05_intersect_arpa.gpkg',
    #    './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_arpa.tif_merged_compressed.tif'),
    # 'non-arpa_split': (
    #    './data/aois/non-Arpa_nonoverlapping-in-m.gpkg',
    #    './data/subwatersheds/hybas_sa_lev05_intersect_non-arpa.gpkg',
    #    './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif_merged_compressed.tif'),
    # 'union-arpa-non': (
    #    './data/aois/union_arpa_non-arpa.gpkg',
    #    './data/subwatersheds/hybas_sa_lev05_intersect_arpa_non-arpa.gpkg',
    #    './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif_merged_compressed.tif'),
    # 'non-arpa': (
    #     './data/aois/non-Arpa_nonoverlapping-in-m.gpkg',
    #     './data/subwatersheds/hybas_sa_lev05_intersect_non-arpa.gpkg',
    #     './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_non-arpa.tif_merged_compressed.tif'),
    # 'Arctic_si': (
    #     './data/aois/Arctic.gpkg',
    #     './data/subwatersheds/hybas_si_lev05_intersect_Arctic_si.gpkg',
    #     './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_si_lev05_intersect_Arctic_si.tif_merged_compressed.tif'),
    # 'Arctic_ar': (
    #     './data/aois/Arctic.gpkg',
    #     './data/subwatersheds/hybas_ar_lev05_intersect_Arctic_ar.gpkg',
    #     './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_ar_lev05_intersect_Arctic_ar.tif_merged_compressed.tif'),
    # 'arpa': (
    #     './data/aois/Arpa_nonoverlapping-in-m.gpkg',
    #     './data/subwatersheds/hybas_sa_lev05_intersect_arpa.gpkg',
    #     './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_arpa.tif_merged_compressed.tif'),
    #     ERROR 1: database disk image is malformed
    # 'Colombia': (
    #     './data/aois/Colombia.gpkg',
    #     './data/subwatersheds/hybas_sa_lev05_intersect_Colombia.gpkg',
    #     './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_Colombia.tif_merged_compressed.tif'),
    # 'NGP': (
    #     './data/aois/NGP.gpkg',
    #     './data/subwatersheds/hybas_na_lev05_intersect_NGP.gpkg',
    #     './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_na_lev05_intersect_NGP.tif_merged_compressed.tif'),
    # 'Peru': (
    #     './data/aois/Peru.gpkg',
    #     './data/subwatersheds/hybas_sa_lev05_intersect_Peru.gpkg',
    #     './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_Peru.tif_merged_compressed.tif'),
    # 'RGRB': (
    #     './data/aois/RGRB.gpkg',
    #     './data/subwatersheds/hybas_na_lev05_intersect_RGRB.gpkg',
    #     './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_na_lev05_intersect_RGRB.tif_merged_compressed.tif'),
    # 'Tapajos': (
    #     './data/aois/Tapajos.gpkg',
    #     './data/subwatersheds/hybas_sa_lev05_intersect_Tapajos.gpkg',
    #     './data/dem_rasters/merged_rasters/JAXA_ALOS_AW3D30_V3_2_hybas_sa_lev05_intersect_Tapajos.tif_merged_compressed.tif'),
}

OUTPUT_DIR = './people_ds_results_900_v2'
for dirpath in [OUTPUT_DIR,]:
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
        base_dem_raster_path, POP_PIXEL_SIZE, clipped_dem_path,
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
        working_dir=local_workspace_dir,
        max_pixel_fill_count=10000)

    routing.flow_dir_mfd(
        (filled_dem_path, 1), target_flow_dir_path,
        working_dir=local_workspace_dir)

    return clipped_dem_path


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
        result[mask_array <= 0] = nodata
        return result

    working_dir_path = os.path.dirname(target_masked_path)
    aligned_raster_path_list = [
        os.path.join(working_dir_path, analysis_id, f'%s_{analysis_id}_aligned%s' % os.path.splitext(os.path.basename(path)))
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


def create_circular_kernel(kernel_path, buffer_size_in_px):
    diameter = buffer_size_in_px * 2 + 1
    kernel_array = numpy.zeros((diameter, diameter), dtype=numpy.float32)
    cx, cy = buffer_size_in_px, buffer_size_in_px

    for i in range(diameter):
        for j in range(diameter):
            if (i - cx)**2 + (j - cy)**2 <= buffer_size_in_px**2:
                kernel_array[i, j] = 1.0

    driver = gdal.GetDriverByName("GTiff")
    out_raster = driver.Create(kernel_path, diameter, diameter, 1, gdal.GDT_Float32)
    out_raster.GetRasterBand(1).WriteArray(kernel_array)
    out_raster.FlushCache()
    out_raster = None


def main():
    """Entry point."""
    task_graph = taskgraph.TaskGraph(OUTPUT_DIR, os.cpu_count(), reporting_interval=10.0)
    kernel_task_map = {}
    result = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(dict)))
    file = open('log.txt', 'w')
    clipped_dem_work_list = []
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
            ignore_path_list=[aoi_vector_path, reprojected_aoi_vector_path],
            target_path_list=[reprojected_aoi_vector_path],
            task_name=f'reproject {analysis_id}')
        flow_dir_path = os.path.join(local_workspace_dir, f'{analysis_id}_mfd_flow_dir.tif')
        flow_dir_task = task_graph.add_task(
            func=calc_flow_dir,
            args=(analysis_id, dem_raster_path, subwatershed_vector_path, flow_dir_path),
            target_path_list=[flow_dir_path],
            store_result=True,
            task_name=f'calculate flow dir for {analysis_id}')
        clipped_dem_work_list.append((analysis_id, (local_workspace_dir, flow_dir_path, reprojected_aoi_vector_path, subwatershed_vector_path, aoi_raster_mask_path, flow_dir_task)))
    for analysis_id, (local_workspace_dir, flow_dir_path, reprojected_aoi_vector_path, subwatershed_vector_path, aoi_raster_mask_path, flow_dir_task) in clipped_dem_work_list:
        clipped_dem_path = flow_dir_task.get()
        print(analysis_id)
        rasterize_task = task_graph.add_task(
            func=rasterize,
            args=(reprojected_aoi_vector_path, clipped_dem_path, aoi_raster_mask_path),
            dependent_task_list=[reproject_task],
            ignore_path_list=[reprojected_aoi_vector_path],
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

        for buffer_size_in_px, buffer_size_in_m in BUFFER_AMOUNTS_IN_PIXELS_M:
            if buffer_size_in_px > 0:
                buffered_downstream_flow_mask_path = f'%s_{buffer_size_in_m}m%s' % os.path.splitext(aoi_downstream_flow_mask_path)
                # make a kernel raster that is a circle kernel that's all 1s within buffer_size_in_px from the center
                # dimensions should be buffer_size_in_px*2+1 X buffer_size_in_px*2+1
                # no projection is necessary
                kernel_path = os.path.join(local_workspace_dir, f'kernel_{buffer_size_in_m}.tif')
                convolve_dependent_task_list = [flow_accum_task]
                if kernel_path not in kernel_task_map:
                    kernel_task = task_graph.add_task(
                        func=create_circular_kernel,
                        args=(kernel_path, buffer_size_in_px),
                        target_path_list=[kernel_path],
                        task_name=f'kernel for {kernel_path}')
                    kernel_task_map[kernel_path] = kernel_task
                    convolve_dependent_task_list.append(kernel_task)
                else:
                    convolve_dependent_task_list.append(kernel_task_map[kernel_path])
                buffer_task = task_graph.add_task(
                    func=geoprocessing.convolve_2d,
                    args=(
                        (aoi_downstream_flow_mask_path, 1), (kernel_path, 1), buffered_downstream_flow_mask_path),
                    kwargs={'n_workers': 1},
                    target_path_list=[buffered_downstream_flow_mask_path],
                    dependent_task_list=convolve_dependent_task_list,
                    task_name=f'buffer {buffered_downstream_flow_mask_path}')
                mask_by_nonzero_and_sum_dependent_task_list = [buffer_task]
            else:
                buffered_downstream_flow_mask_path = aoi_downstream_flow_mask_path
                mask_by_nonzero_and_sum_dependent_task_list = [flow_accum_task]

            for population_raster_path in POPULATION_RASTERS:
                pop_basename = os.path.splitext(os.path.basename(population_raster_path))[0]
                print(f'processing {pop_basename} for {buffer_size_in_m}m buffer')
                masked_population_raster_path = os.path.join(
                    local_workspace_dir, f'{analysis_id}_{buffer_size_in_m}m_{os.path.basename(population_raster_path)}')
                mask_by_nonzero_task = task_graph.add_task(
                    func=mask_by_nonzero_and_sum,
                    args=(
                        f'{analysis_id}_{pop_basename}_{buffer_size_in_m}', population_raster_path, buffered_downstream_flow_mask_path, masked_population_raster_path),
                    dependent_task_list=mask_by_nonzero_and_sum_dependent_task_list,
                    target_path_list=[masked_population_raster_path],
                    store_result=True,
                    task_name=f'{analysis_id}_population mask')
                file.write(
                    f'{analysis_id}_{pop_basename}, {population_raster_path}, {buffered_downstream_flow_mask_path}, {masked_population_raster_path}\n')

                result[analysis_id][pop_basename][buffer_size_in_m] = mask_by_nonzero_task

    for analysis_id in result:
        for pop_basename in result[analysis_id]:
            for buffer_size_in_m in result[analysis_id][pop_basename]:
                result[analysis_id][pop_basename][buffer_size_in_m] = result[analysis_id][pop_basename][buffer_size_in_m].get()

    timestamp = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    output_filename = os.path.join(OUTPUT_DIR, f'pop_results_{timestamp}.csv')
    with open(output_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['analysis_id', 'pop_basename', 'buffer in m', 'value'])
        for analysis_id in result:
            for pop_basename, buffer_val_dict in result[analysis_id].items():
                for buffer, val in buffer_val_dict.items():
                    writer.writerow([analysis_id, pop_basename, buffer, val])
    task_graph.join()
    task_graph.close()
    LOGGER.info(f'all done results at {output_filename}')
    file.close()

if __name__ == '__main__':
    main()
