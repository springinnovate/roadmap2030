"""Dynamic world landcover map puller.
python3 dynamic_world_extractor.py --aoi_vector_path ./NGP_intersected_hybas_na_lev05_v1c.shp --date_ranges 2024-01-01--2024-01-31
"""
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


DATASET_ID = 'JAXA/ALOS/AW3D30/V3_2'
DATASET_CRS = 'EPSG:4326'
DATASET_SCALE = 30
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
        'Extract ALOS DSM.'))
    parser.add_argument(
        '--aoi_vector_paths', nargs='+', help='Paths to vector/shapefiles of areas of interest', required=True)
    parser.add_argument(
        '--status', action='store_true', help='To check task status')
    parser.add_argument(
        '--dataset_scale', type=float, default=DATASET_SCALE, help=(
            f'Override the base scale of {DATASET_SCALE}m to '
            f'whatever you desire.'))
    parser.add_argument(
        '--check_tasks', action='store_true', help="do this if need to guard against duplicate runs")
    parser.add_argument(
        '--create_shared_folder', action='store_true', help='create shared folder')
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

    LOGGER.info(f'reading the vectors at {args.aoi_vector_paths}')
    vector_path_list = [path for path_pattern in args.aoi_vector_paths for path in glob.glob(path_pattern)]
    LOGGER.info(f'processing {len(vector_path_list)} vectors')
    if args.status:
        credentials, project = google.auth.default()
        print(f"Authenticated account: {credentials.service_account_email if hasattr(credentials, 'service_account_email') else 'Unknown account'}")
        # Loop through each task to print its status
        for task in ee.batch.Task.list():
            LOGGER.info(task)
        return

    task_list = []
    for aoi_vector_path in vector_path_list:
        LOGGER.info(f'reading file {aoi_vector_path}')
        aoi_vector = gpd.read_file(aoi_vector_path).to_crs('EPSG:4326')
        total_bounds = aoi_vector.total_bounds  # [minx, miny, maxx, maxy]
        bounding_box = [
            [total_bounds[0], total_bounds[1]],
            [total_bounds[0], total_bounds[3]],
            [total_bounds[2], total_bounds[3]],
            [total_bounds[2], total_bounds[1]],
            [total_bounds[0], total_bounds[1]],
        ]
        LOGGER.info(f'bounds are {bounding_box}')
        dataset = (
            ee.ImageCollection(DATASET_ID)
            .select('DSM')  # Select both DSM and MSK bands
            .filterBounds(ee.Geometry.Rectangle(list(total_bounds)))
        )

        # Combine images using median
        dem_image = dataset.reduce(ee.Reducer.median())

        aoi_basename = os.path.basename(os.path.splitext(aoi_vector_path)[0])
        local_description = f'{DATASET_ID}_{aoi_basename}'
        local_description = local_description.replace('/', '_')

        LOGGER.info(f'Description: "{local_description}"')
        if local_description in existing_descriptions:
            LOGGER.info(f"Task '{local_description}' already in queue. Skipping.")
            continue

        task = ee.batch.Export.image.toCloudStorage(
            image=dem_image,
            description=local_description,
            bucket='ecoshard-root',
            fileNamePrefix=f'roadmap2030/{local_description}.tif',
            region=bounding_box,
            scale=args.dataset_scale,
            crs=DATASET_CRS,
            maxPixels=1e13,
            fileFormat='GeoTIFF'
        )
        task.start()
        task_list.append((f'{aoi_basename}', task))

    for aoi_basename, task in task_list:
        while True:
            status = task.status()['state']
            if status in ['COMPLETED', 'FAILED', 'CANCELLED']:
                break
            LOGGER.info(f'{aoi_basename}: {status}, waiting 10 more seconds')
            time.sleep(10)

        LOGGER.info(f"Task finished with status: {task.status()['state']}")


if __name__ == '__main__':
    main()
