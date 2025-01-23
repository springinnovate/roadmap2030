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

logging.basicConfig(
    level=logging.WARNING,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)


DATASET_ID = 'GOOGLE/DYNAMICWORLD/V1'
DATASET_CRS = 'EPSG:4326'
DATASET_SCALE = 10
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


def parse_monthly_ranges(years):
    return [f"{year}-{str(month).zfill(2)}-01--{year}-{str(month).zfill(2)}-{28 if month == 2 else 30 if month in [4, 6, 9, 11] else 31}"
            for year in years for month in range(1, 13)]


def main():
    parser = argparse.ArgumentParser(description=(
        'Fetch CMIP6 based erosivity given a year or list of years.'))
    parser.add_argument(
        '--aoi_vector_paths', nargs='+', help='Paths to vector/shapefiles of areas of interest', required=True)
    parser.add_argument(
        '--years', nargs='+', type=int, help='List of years to fetch data for', required=True)
    parser.add_argument(
        '--status', action='store_true', help='To check task status')
    parser.add_argument(
        '--dataset_scale', type=float, default=DATASET_SCALE, help=(
            f'Override the base scale of {DATASET_SCALE}m to '
            f'whatever you desire.'))
    parser.add_argument(
        '--guard_tasks', action='store_true', help="do this if need to recreate an image")

    args = parser.parse_args()
    LOGGER.info('about to authenticate')
    authenticate()
    LOGGER.info('authenticate!')

    existing_descriptions = set()
    if args.guard_tasks:
        existing_tasks = ee.batch.Task.list()
        allowed_states = {"READY", "RUNNING", "COMPLETED"}
        for t in existing_tasks:
            cfg = t.config
            if t.status()['state'] in allowed_states and cfg and 'description' in cfg:
                existing_descriptions.add(cfg['description'])

    vector_path_list = [path for path_pattern in args.aoi_vector_paths for path in glob.glob(path_pattern)]
    LOGGER.info(f'processing {len(vector_path_list)} vectors')
    if args.status:
        # Loop through each task to print its status
        for task in ee.batch.Task.list():
            LOGGER.info(task)
        return

    # Generate date ranges for each year
    #parse_monthly_ranges(args.years)
    date_ranges = [f'{year}-01-01--{year}-12-31' for year in args.years]
    task_list = []
    for aoi_vector_path in vector_path_list:
        aoi_vector = gpd.read_file(aoi_vector_path).to_crs('EPSG:4326')
        LOGGER.debug(aoi_vector.bounds)
        return
        bounds = aoi_vector.bounds.iloc[0]
        bounding_box = [[bounds.minx, bounds.miny], [bounds.minx, bounds.maxy],
                        [bounds.maxx, bounds.maxy], [bounds.maxx, bounds.miny],
                        [bounds.minx, bounds.miny]]
        LOGGER.debug(bounding_box)

        for date_range in date_ranges:
            start_date, end_date = date_range.split('--')

            dataset = (ee.ImageCollection(DATASET_ID)
                       .select('label')
                       .filterBounds(ee.Geometry.Rectangle([bounds.minx, bounds.miny, bounds.maxx, bounds.maxy]))
                       .filterDate(start_date, end_date))

            average_landcover_image = dataset.reduce(ee.Reducer.mode())
            aoi_basename = os.path.basename(os.path.splitext(aoi_vector_path)[0])
            local_description = f'{DATASET_ID}_{aoi_basename}_{start_date}--{end_date}'
            local_description = local_description.replace('/', '_')

            LOGGER.info(f'Description: "{local_description}"')
            if not args.force and local_description in existing_descriptions:
                LOGGER.info(f"Task '{local_description}' already in queue. Skipping.")
                continue

            task = ee.batch.Export.image.toDrive(
                image=average_landcover_image,
                description=local_description,
                folder=EXPORT_DRIVE_FOLDER,
                scale=args.dataset_scale,
                crs=DATASET_CRS,
                region=bounding_box,
                maxPixels=1e10,
                fileFormat='GeoTIFF',
            )
            task.start()
            task_list.append((f'{start_date}-{end_date}', task))

    for date_str, task in task_list:
        while True:
            status = task.status()['state']
            if status in ['COMPLETED', 'FAILED', 'CANCELLED']:
                break
            LOGGER.info(f'{date_str}: {status}, waiting 10 more seconds')
            time.sleep(10)

        LOGGER.info(f"Task finished with status: {task.status()['state']}")


if __name__ == '__main__':
    main()
174332704