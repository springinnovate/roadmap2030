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

DATASET_ID = 'MODIS/061/MOD16A2'
DATASET_BAND = 'ET'  # Evapotranspiration band name in MOD16A2
DATASET_CRS = 'EPSG:4326'
DATASET_SCALE = 500  # MOD16A2 is typically 500m
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
    monthly_ranges = []
    for year in years:
        for month in range(1, 13):
            # Handle days in month (simple version)
            if month == 2:
                day_end = 28
            elif month in [4, 6, 9, 11]:
                day_end = 30
            else:
                day_end = 31
            start_str = f"{year}-{str(month).zfill(2)}-01"
            end_str = f"{year}-{str(month).zfill(2)}-{day_end}"
            monthly_ranges.append(f"{start_str}--{end_str}")
    return monthly_ranges


def main():
    parser = argparse.ArgumentParser(description=(
        'Fetch monthly MOD16A2 Evapotranspiration in mm/month for a given list of years.'))
    parser.add_argument(
        '--aoi_vector_paths', nargs='+', help='Paths to vector/shapefiles of areas of interest', required=True)
    parser.add_argument(
        '--years', nargs='+', type=int, help='List of years to fetch data for', required=True)
    parser.add_argument(
        '--status', action='store_true', help='To check task status')
    parser.add_argument(
        '--dataset_scale', type=float, default=DATASET_SCALE, help=(
            f'Override the base scale of {DATASET_SCALE}m if desired.'))
    parser.add_argument(
        '--check_tasks', action='store_true', help="Protect against duplicate tasks")

    args = parser.parse_args()
    LOGGER.info('about to authenticate')
    authenticate()
    LOGGER.info('authenticated!')

    existing_descriptions = set()
    if args.check_tasks:
        existing_tasks = ee.batch.Task.list()
        allowed_states = {"READY", "RUNNING", "COMPLETED"}
        for t in existing_tasks:
            cfg = t.config
            if t.status()['state'] in allowed_states and cfg and 'description' in cfg:
                existing_descriptions.add(cfg['description'])

    vector_path_list = [
        path for path_pattern in args.aoi_vector_paths
        for path in glob.glob(path_pattern)
    ]
    LOGGER.info(f'processing {len(vector_path_list)} vectors')

    if args.status:
        for task in ee.batch.Task.list():
            LOGGER.info(task)
        return

    # Generate monthly date ranges for each year
    monthly_date_ranges = parse_monthly_ranges(args.years)

    task_list = []
    for aoi_vector_path in vector_path_list:
        aoi_vector = gpd.read_file(aoi_vector_path).to_crs(DATASET_CRS)
        total_bounds = aoi_vector.total_bounds  # [minx, miny, maxx, maxy]
        bounding_box = [
            [total_bounds[0], total_bounds[1]],
            [total_bounds[0], total_bounds[3]],
            [total_bounds[2], total_bounds[3]],
            [total_bounds[2], total_bounds[1]],
            [total_bounds[0], total_bounds[1]],
        ]

        for date_range in monthly_date_ranges:
            start_date, end_date = date_range.split('--')

            # Pull the ET data for the monthly range
            et_collection = (ee.ImageCollection(DATASET_ID)
                             .select(DATASET_BAND)
                             .filterBounds(ee.Geometry.Rectangle(list(total_bounds)))
                             .filterDate(start_date, end_date))

            # Each ET image is in kg/m^2/8day with a 0.1 scale factor => multiply by 0.1
            # 1 kg/m^2 = 1 mm, so final ET in mm is just sum * 0.1
            monthly_et = et_collection.map(
                lambda img: img.multiply(0.1)
            ).sum().rename('ET_mm_month')

            aoi_basename = os.path.basename(os.path.splitext(aoi_vector_path)[0])
            local_description = f'MOD16A2_{aoi_basename}_{start_date}--{end_date}'
            local_description = local_description.replace('/', '_')

            LOGGER.info(f'Description: "{local_description}"')
            if local_description in existing_descriptions:
                LOGGER.info(f"Task '{local_description}' already in queue. Skipping.")
                continue

            export_task = ee.batch.Export.image.toCloudStorage(
                image=monthly_et,
                description=local_description,
                bucket='ecoshard-root',
                fileNamePrefix=f'roadmap2030/{local_description}.tif',
                region=bounding_box,
                scale=args.dataset_scale,
                crs=DATASET_CRS,
                maxPixels=1e13,
                fileFormat='GeoTIFF'
            )
            export_task.start()
            task_list.append((f'{start_date}-{end_date}', export_task))

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
