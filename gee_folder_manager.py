"""Manage google drive programatically"""
import argparse
import logging
import os
import sys

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


def delete_file_callback(request_id, response, exception):
    if exception is not None:
        print(f"Failed to delete file (request_id: {request_id}): {exception}")
    elif isinstance(response, dict) and 'id' in response:
        print(f"Deleted file: {response['id']}")
    else:
        print(f"Unexpected response for request_id {request_id}: {response}")


def delete_all_files_in_folder(folder_id):
    gee_key_path = os.environ['GEE_KEY_PATH']
    SCOPES = ['https://www.googleapis.com/auth/drive']

    # Authenticate the service account
    credentials = Credentials.from_service_account_file(gee_key_path, scopes=SCOPES)
    service = build('drive', 'v3', credentials=credentials)
    try:
        query = f"'{folder_id}' in parents and trashed = false"
        response = service.files().list(q=query).execute()
        files = response.get('files', [])

        batch = service.new_batch_http_request(callback=delete_file_callback)
        for file in files:
            batch.add(service.files().delete(fileId=file['id']))
        batch.execute()
        print("Batch deletion completed.")

        print("All files in the folder have been deleted.")

    except HttpError as error:
        print(f"An error occurred: {error}")


def make_shared_folder(folder_name):
    # Path to your service account key file
    gee_key_path = os.environ['GEE_KEY_PATH']
    SCOPES = ['https://www.googleapis.com/auth/drive']

    # Authenticate the service account
    credentials = Credentials.from_service_account_file(gee_key_path, scopes=SCOPES)
    drive_service = build('drive', 'v3', credentials=credentials)

    # Create a folder in the service account's Drive
    folder_metadata = {
        'name': folder_name,
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
    parser = argparse.ArgumentParser(description='folder manager')
    parser.add_argument('--delete', metavar='FOLDER_ID', type=str, help='Specify the folder ID to delete.')
    parser.add_argument('--create', metavar='FOLDER_NAME', type=str, help='Specify the name of the folder to create.')
    args = parser.parse_args()

    if args.delete:
        delete_all_files_in_folder(args.delete)
    if args.create:
        make_shared_folder(args.create)

    LOGGER.info('all done')


if __name__ == '__main__':
    main()
