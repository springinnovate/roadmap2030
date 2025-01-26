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


def make_shared_folder():
    from googleapiclient.discovery import build
    from google.oauth2.service_account import Credentials

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


def main():
    authenticate()

    existing_tasks = ee.batch.Task.list()
    allowed_states = {"READY", "RUNNING"}
    for t in existing_tasks:
        if t.status()['state'] in allowed_states:
            t.cancel()
            print(t.id)


if __name__ == '__main__':
    main()
