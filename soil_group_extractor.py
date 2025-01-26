import glob
import time
import argparse
import logging
import os
import sys

import ee
import geopandas as gpd
import google.auth


logging.basicConfig(
    level=logging.WARNING,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)


DATASET_ID = 'projects/soilgrids-isric/soilgrids_global'
DATASET_CRS = 'EPSG:4326'
DATASET_SCALE = 1000
EXPORT_DRIVE_FOLDER = 'gee_exports'


def authenticate():
    try:
        ee.Initialize()
        return
    except Exception:
        pass

    try:
        gee_key_path = os.environ['GEE_KEY_PATH']
        credentials = ee.ServiceAccountCredentials(None, gee_key_path)
        ee.Initialize(credentials)
        return
    except Exception:
        pass

    try:
        ee.Authenticate()
        ee.Initialize()
        return
    except Exception:
        pass

    ee.Initialize()


def main():
    parser = argparse.ArgumentParser(description=(
        'Extract CHIRPS daily precipitation and optionally classify SoilGrids.'))
    parser.add_argument(
        '--aoi_vector_paths', nargs='+', help='Paths to vector/shapefiles of areas of interest', required=True)
    parser.add_argument(
        '--status', action='store_true', help='To check task status')
    parser.add_argument(
        '--dataset_scale', type=float, default=DATASET_SCALE, help=(
            f'Override the base scale of {DATASET_SCALE}m to whatever you desire.'))
    parser.add_argument(
        '--check_tasks', action='store_true', help="do this if need to guard against duplicate runs")
    args = parser.parse_args()
    authenticate()

    existing_descriptions = set()
    if args.check_tasks:
        existing_tasks = ee.batch.Task.list()
        allowed_states = {"READY", "RUNNING", "COMPLETED"}
        for t in existing_tasks:
            cfg = t.config
            if t.status()['state'] in allowed_states and cfg and 'description' in cfg:
                existing_descriptions.add(cfg['description'])

    vector_path_list = [path for path_pattern in args.aoi_vector_paths for path in glob.glob(path_pattern)]

    if args.status:
        credentials, project = google.auth.default()
        print(f"Authenticated account: {credentials.service_account_email if hasattr(credentials, 'service_account_email') else 'Unknown account'}")
        for task in ee.batch.Task.list():
            LOGGER.info(task)
        return

    task_list = []
    for aoi_vector_path in vector_path_list:
        aoi_vector = gpd.read_file(aoi_vector_path).to_crs('EPSG:4326')
        total_bounds = aoi_vector.total_bounds
        bounding_box = [
            [total_bounds[0], total_bounds[1]],
            [total_bounds[0], total_bounds[3]],
            [total_bounds[2], total_bounds[3]],
            [total_bounds[2], total_bounds[1]],
            [total_bounds[0], total_bounds[1]],
        ]
        aoi_basename = os.path.basename(os.path.splitext(aoi_vector_path)[0])

        clay_img = ee.Image('projects/soilgrids-isric/clay_mean').select('clay_0-5cm_mean')
        sand_img = ee.Image('projects/soilgrids-isric/sand_mean').select('sand_0-5cm_mean')
        hydrologic_group = sand_img.expression(
            '(sand > 70 && clay < 10) ? 1 : '
            '(sand > 50 && sand <= 70 && clay < 20) ? 2 : '
            '(sand > 20 && sand <= 50 && clay < 40) ? 3 : 4',
            {
                'sand': sand_img.divide(10),  # SoilGrids is in g/kg; convert to %
                'clay': clay_img.divide(10)
            }
        )
        soil_hsg_description = f'soil_hydrologic_groups_{aoi_basename}'
        if soil_hsg_description not in existing_descriptions:
            soil_task = ee.batch.Export.image.toCloudStorage(
                image=hydrologic_group,
                description=soil_hsg_description,
                bucket='ecoshard-root',
                fileNamePrefix=f'roadmap2030/{soil_hsg_description}.tif',
                region=bounding_box,
                scale=args.dataset_scale,
                crs=DATASET_CRS,
                maxPixels=1e13,
                fileFormat='GeoTIFF'
            )
            soil_task.start()
            task_list.append(('SOIL_HSG', soil_task))

    for id_str, task in task_list:
        while True:
            status = task.status()['state']
            if status in ['COMPLETED', 'FAILED', 'CANCELLED']:
                break
            LOGGER.info(f'{id_str}: {status}, waiting 10 more seconds')
            time.sleep(10)
        LOGGER.info(f"Task finished with status: {task.status()['state']}")


if __name__ == '__main__':
    main()
