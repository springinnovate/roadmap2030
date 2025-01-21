import datetime
import csv
import collections
import numpy
from osgeo import osr
from osgeo import gdal
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

DEM_PATH = r"D:/repositories/downstream-beneficiaries/workspace/global_dem_3s_md5_22d0c3809af491fa09d03002bdf09748/global_dem_3s"

BASE_RASTER_LOOKUP = {
    'eii': r"D:\repositories\data_platform\Nature\eii_padj_v5140524_epsg_3395.tif",
    #'sed_export_change': r"D:\repositories\roadmap2030\data\ndv_0.0_sed_export_marineESA_2020-1992_change_md5_0ab0cf.tif",
    #'cv_habitat_change': r"D:\repositories\roadmap2030\data\ndv_0.0_cv_habitat_value_marESA2020-1992_change_md5_1643a7.tif",
    #'n_export_change': r"D:\repositories\roadmap2030\data\ndv_0.0_n_export_marineESA_2020-1992_change_val_md5_18a2b3.tif",
    #'realized_pollination_on_ag_change': r"D:\repositories\roadmap2030\data\ndv_0.0_realized_pollination_on_ag_marESA_2020-1992_fullchange_md5_8e63e2.tif",
    #'sed_deposition_change': r"D:\repositories\roadmap2030\data\ndv_0.0_sed_deposition_marineESA_2020-1992_change_md5_d23c49.tif",
    #'coastal_reference': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\coastal_risk_Sc3v1_habitat_value_md5_e889c2dbc5783fc4c782fbd3b473d7de.tif",
    #'coastal_change': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\coastal_risk_tnc_esa2020_change_esa1992_md5_ea900e.tif",
    #'coastal_2020': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\coastal_risk_tnc_esa2020_value_md5_f9f644.tif",
    #'nitrogen_reference_full': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\n_export_sc3v2pnvall_compressed_md5_09bc65fe1cd54b518cde859f57513d8c.tif",
    #'nitrogen_reference': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\n_export_sc3v1pnvnoag_compressed_md5_bd5a856e0c1f76b2e8898f533ec20659.tif",
    #'nitrogen_change': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\n_export_tnc_2020-1992_change_val_md5_18a2b3.tif",
    #'nitrogen_2020': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\n_export_tnc_esa2020_compressed_md5_1d3c17.tif",
    #'pollination_change': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\pollination_on_ag_marESA_2020-1992_fullchange_md5_8e63e2.tif",
    #'pollination_reference': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\pollination_ppl_fed_on_ag_10s_Sc3v1_PNVnoag.tif",
    #'pollination_2020v2': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\pollination_ppl_fed_on_ag_10s_tnc_esa2020ag_compressed_md5_8b5ee8.tif",
    #'pollination_2020v1': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\polllination_on_ag_ESA2020mar_md5_da610a.tif",
    #'sediment_reference_full': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\sed_export_pnv_compressed_md5_a1faed.tif"
    #'sediment_reference': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\sed_export_sc3v1pnvnoag_compressed_md5_2783ee50e908a763622d3167669b60bc.tif",
    #'sediment_change': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\sed_export_tnc_ESA_2020-1992_change_md5_0ab0cf.tif",
    #'sediment_2020': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\sed_export_tnc_ESA_2020_compressed_md5_a988c0.tif"
}


VECTOR_PATH_LOOKUP = {
    'non-arpa': r"D:\repositories\roadmap2030\data\non-arpa-projected-in-m.gpkg",
    'arpa': r"D:\repositories\roadmap2030\data\arpa-projected-in-m.gpkg",
    'colombia': r"D:\repositories\roadmap2030\data\Colombia.gpkg",
    'peru': r"D:\repositories\roadmap2030\data\Peru.gpkg",
    'tapajos': r"D:\repositories\roadmap2030\data\Tapajos.gpkg",
    'NGP': r"D:\repositories\roadmap2030\data\NGP.gpkg",
    'RGBR': r"D:\repositories\roadmap2030\data\RGRB.gpkg",
    'Arctic': r"D:\repositories\roadmap2030\data\Arctic.gpkg"
}

OUTPUT_DIR = './results'
CLIPPED_DIR = os.path.join(OUTPUT_DIR, 'clipped')
for dirpath in [OUTPUT_DIR, CLIPPED_DIR]:
    os.makedirs(dirpath, exist_ok=True)


def vector_area_in_ha(vector_path):
    dataset = gdal.OpenEx(vector_path, gdal.OF_VECTOR)
    layer = dataset.GetLayer()

    source_srs = layer.GetSpatialRef()
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(54009)

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
    LOGGER.info(f'creating subset of {name}')
    subset_gdf = gdf[gdf["Name"] == name]
    subset_gdf.to_file(target_vector_path, driver="GPKG")
    LOGGER.info(f'done with subset of {name}')


def clip_raster(base_raster_path, summary_vector_path, temp_clip_path):
    base_raster_info = geoprocessing.get_raster_info(base_raster_path)
    summary_vector_info = geoprocessing.get_vector_info(summary_vector_path)
    target_pixel_size = base_raster_info['pixel_size']
    base_vector_bb = summary_vector_info['bounding_box']

    target_bb = geoprocessing.transform_bounding_box(
        base_vector_bb, summary_vector_info['projection_wkt'],
        base_raster_info['projection_wkt'])

    geoprocessing.warp_raster(
        base_raster_path, target_pixel_size, temp_clip_path,
        'near', target_bb=target_bb, vector_mask_options={
            'mask_vector_path': summary_vector_path,
            'all_touched': True})

def get_stats(raster_path):
    r = gdal.OpenEx(raster_path)
    b = r.GetRasterBand(1)
    array = b.ReadAsArray()
    nodata = b.GetNoDataValue()
    array = array[array != nodata]

    stats = {
        'min': numpy.min(array),
        'max': numpy.max(array),
        'sum': numpy.sum(array),
        'mean': numpy.mean(array)
    }
    return stats

def dump_results_to_csv(results, vector_path_lookup, csv_path):
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["vector_id", "raster_name", "area_ha", "min", "max", "mean", "sum"])
        for vector_id, info_dict in results.items():
            area_ha = info_dict.get("area_ha", None)
            if area_ha is not None:
                writer.writerow([
                    vector_id,
                    "",         # raster_name is empty
                    area_ha,
                    "", "", "", ""  # no min/max/mean/sum for area
                ])

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

                writer.writerow([
                    vector_id,
                    raster_basename,  # raster_name
                    "",               # area_ha is empty here
                    r_min,
                    r_max,
                    r_mean,
                    r_sum
                ])


def main():
    """Entry point."""
    print(os.cpu_count())
    task_graph = taskgraph.TaskGraph(OUTPUT_DIR, os.cpu_count(), reporting_interval=10.0)
    results = collections.defaultdict(lambda: collections.defaultdict(dict))
    for vector_id, vector_path in VECTOR_PATH_LOOKUP.items():
        results[vector_id]['area_ha'] = vector_area_in_ha(vector_path)
        LOGGER.info(f'processing {vector_id}')
        for raster_basename, raster_path in BASE_RASTER_LOOKUP.items():
            if not os.path.exists(raster_path):
                raise RuntimeError(f'{raster_path} not found')
            LOGGER.info(f'clipping {raster_basename} to {vector_id}')
            clipped_raster_path = os.path.join(
                CLIPPED_DIR, f'{vector_id}_{raster_basename}.tif')
            clipped_task = task_graph.add_task(
                func=clip_raster,
                args=(raster_path, vector_path, clipped_raster_path),
                ignore_path_list=[vector_path],
                target_path_list=[clipped_raster_path],
                task_name=f'clipping {raster_path} to {vector_path}')

            stats_task = task_graph.add_task(
                func=get_stats,
                args=(clipped_raster_path,),
                dependent_task_list=[clipped_task],
                store_result=True,
                task_name=f'stats for {raster_path}')
            results[vector_id][raster_basename]['stats'] = stats_task

    for vector_id in VECTOR_PATH_LOOKUP:
        for raster_basename in BASE_RASTER_LOOKUP:
            stats_task = results[vector_id][raster_basename]['stats']
            results[vector_id][raster_basename] = {}
            for fn_str, value in stats_task.get().items():
                results[vector_id][raster_basename][fn_str] = value

    task_graph.join()
    print(results)
    timestamp = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    output_filename = f'results_{timestamp}.csv'
    dump_results_to_csv(results, VECTOR_PATH_LOOKUP, output_filename)
    print(f'all done -- results in {output_filename}!')


if __name__ == '__main__':
    main()
