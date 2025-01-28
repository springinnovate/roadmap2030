"""Dynamic world landcover map puller.
python3 chirps_extractor.py --aoi_vector_path ./NGP_intersected_hybas_na_lev05_v1c.shp --years 20204
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
from googleapiclient.errors import HttpError
from googleapiclient.discovery import build
from google.oauth2.service_account import Credentials


logging.basicConfig(
    level=logging.WARNING,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)


DATASET_ID = 'UCSB-CHG/CHIRPS/DAILY'
DATASET_SHORT_ID = 'CHIRPS'
DATASET_CRS = 'EPSG:4326'
DATASET_SCALE = 5566
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


def delete_all_files_in_folder(folder_id):
    gee_key_path = os.environ['GEE_KEY_PATH']
    SCOPES = ['https://www.googleapis.com/auth/drive']

    # Authenticate the service account
    credentials = Credentials.from_service_account_file(gee_key_path, scopes=SCOPES)
    service = build('drive', 'v3', credentials=credentials)
    try:
        # List all files in the folder
        query = f"'{folder_id}' in parents and trashed = false"
        response = service.files().list(q=query).execute()
        files = response.get('files', [])

        for file in files:
            file_id = file['id']
            file_name = file['name']
            # Delete the file
            service.files().delete(fileId=file_id).execute()
            print(f"Deleted file: {file_name} (ID: {file_id})")

        print("All files in the folder have been deleted.")

    except HttpError as error:
        print(f"An error occurred: {error}")


def make_shared_folder():
    # Path to your service account key file
    gee_key_path = os.environ['GEE_KEY_PATH']
    SCOPES = ['https://www.googleapis.com/auth/drive']

    # Authenticate the service account
    credentials = Credentials.from_service_account_file(gee_key_path, scopes=SCOPES)
    drive_service = build('drive', 'v3', credentials=credentials)

    # Create a folder in the service account's Drive
    folder_metadata = {
        'name': EXPORT_DRIVE_FOLDER,
        'mimeType': 'application/vnd.google-apps.folder'
    }
    folder = drive_service.files().create(body=folder_metadata, fields='id').execute()
    folder_id = folder.get('id')

    print(f"Folder created with ID: {folder_id}")

    # Share the folder with your personal Google account
    user_permission = {
        'type': 'user',
        'role': 'writer',
        'emailAddress': 'richpsharp@gmail.com'
    }
    drive_service.permissions().create(
        fileId=folder_id,
        body=user_permission,
        fields='id'
    ).execute()

    print("Folder shared successfully!")

def parse_monthly_ranges(years):
    return [f"{year}-{str(month).zfill(2)}-01--{year}-{str(month).zfill(2)}-{28 if month == 2 else 30 if month in [4, 6, 9, 11] else 31}"
            for year in years for month in range(1, 13)]


def main():
    parser = argparse.ArgumentParser(description=(
        'Extract ALOS DSM.'))
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
        '--check_tasks', action='store_true', help="do this if need to guard against duplicate runs")
    parser.add_argument(
        '--create_shared_folder', action='store_true', help='create shared folder')
    args = parser.parse_args()
    authenticate()

    if args.create_shared_folder:
        make_shared_folder()
        return

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
        # Loop through each task to print its status
        for task in ee.batch.Task.list():
            LOGGER.info(task)
        return

    # Generate date ranges for each year
    date_ranges = parse_monthly_ranges(args.years)
    task_list = []
    for aoi_vector_path in vector_path_list:
        aoi_vector = gpd.read_file(aoi_vector_path).to_crs('EPSG:4326')
        aoi_vector = gpd.read_file(aoi_vector_path).to_crs('EPSG:4326')
        total_bounds = aoi_vector.total_bounds  # [minx, miny, maxx, maxy]
        bounding_box = [
            [total_bounds[0], total_bounds[1]],
            [total_bounds[0], total_bounds[3]],
            [total_bounds[2], total_bounds[3]],
            [total_bounds[2], total_bounds[1]],
            [total_bounds[0], total_bounds[1]],
        ]
        for date_range in date_ranges:
            start_date, end_date = date_range.split('--')

            dataset = (
                ee.ImageCollection(DATASET_ID)
                .select('precipitation')
                .filterDate(start_date, end_date)
                .filterBounds(ee.Geometry.Rectangle(list(total_bounds)))
            )

            # Combine images using median
            precip_image = dataset.reduce(ee.Reducer.median())

            aoi_basename = os.path.basename(os.path.splitext(aoi_vector_path)[0])
            precip_description = f'{DATASET_SHORT_ID}_precipitation_{aoi_basename}_{start_date}--{end_date}'
            precip_description = precip_description.replace('/', '_')

            LOGGER.info(f'Description: "{precip_description}"')
            if precip_description in existing_descriptions:
                LOGGER.info(f"Task '{precip_description}' already in queue. Skipping.")
                continue

            task = ee.batch.Export.image.toCloudStorage(
                image=precip_image,
                description=precip_description,
                bucket='ecoshard-root',
                fileNamePrefix=f'roadmap2030/{precip_description}.tif',
                region=bounding_box,
                scale=args.dataset_scale,
                crs=DATASET_CRS,
                maxPixels=1e13,
                fileFormat='GeoTIFF'
            )
            task.start()
            task_list.append((f'{start_date}-{end_date}', task))

            # do n_events
            # Count days > 1mm
            def mask_precip(image):
                return image.gt(1).rename('day_over_1mm')

            count_image = dataset.map(mask_precip).sum().rename('n_events').toInt16()
            n_events_description = f'{DATASET_SHORT_ID}_n_events_{aoi_basename}_{start_date}--{end_date}'
            n_events_description = n_events_description.replace('/', '_')

            LOGGER.info(f'Description: "{n_events_description}"')
            if n_events_description in existing_descriptions:
                LOGGER.info(f"Task '{n_events_description}' already in queue. Skipping.")
                continue

            task = ee.batch.Export.image.toCloudStorage(
                image=count_image,
                description=n_events_description,
                bucket='ecoshard-root',
                fileNamePrefix=f'roadmap2030/{n_events_description}.tif',
                region=bounding_box,
                scale=args.dataset_scale,
                crs=DATASET_CRS,
                maxPixels=1e13,
                fileFormat='GeoTIFF'
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
