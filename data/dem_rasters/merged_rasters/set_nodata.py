import os
from osgeo import gdal

def set_nodata_value(directory, nodata_value=0):
    """
    Set the NoData value for all .tif files in the specified directory.

    Args:
        directory (str): Path to the directory containing .tif files.
        nodata_value (int or float): The NoData value to set.
    """
    for filename in os.listdir(directory):
        if filename.endswith(".tif"):
            file_path = os.path.join(directory, filename)
            print(f"Processing {file_path}...")

            # Open the dataset
            dataset = gdal.Open(file_path, gdal.GA_Update)
            if dataset is None:
                print(f"Failed to open {file_path}")
                continue

            # Get the first band and set NoData value
            band = dataset.GetRasterBand(1)
            band.SetNoDataValue(nodata_value)
            band.FlushCache()  # Ensure changes are written to disk

            # Close the dataset
            dataset = None
            print(f"Set NoData value to {nodata_value} for {file_path}")
    print("Finished processing all .tif files.")

# Example usage
if __name__ == "__main__":
    current_directory = os.getcwd()
    set_nodata_value(current_directory, nodata_value=0)
