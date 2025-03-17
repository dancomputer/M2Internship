import geopandas as gpd
from shapely.geometry import Polygon
from PIL import Image
import numpy as np

def open_tiff_as_geodataframe(tiff_path, crs=None):
    """
    Opens a TIFF file and converts it into a GeoDataFrame.  Each pixel
    becomes a polygon in the GeoDataFrame.  Assumes the TIFF represents
    a raster where each pixel has a single value (e.g., a single band
    grayscale or categorical raster).  Handles multi-band TIFFs by
    taking only the first band.  Provides options for handling NoData values.

    Args:
        tiff_path (str): Path to the TIFF file.
        crs (str, optional):  Coordinate Reference System (CRS) string.
            If provided, sets the CRS of the resulting GeoDataFrame.
            Defaults to None (no CRS set).  A common example would be
            'EPSG:4326'.

    Returns:
        geopandas.GeoDataFrame: A GeoDataFrame where each row represents a
        pixel. The 'geometry' column contains the polygon, and 'value' column
        contains the pixel's value.  Returns None if the TIFF cannot be opened.

    Raises:
        ImportError: If `rasterio` is not installed.
        Exception: For errors during TIFF processing.
    """
    try:
        import rasterio
    except ImportError:
        raise ImportError("rasterio is required to open TIFF files.  Install it with 'pip install rasterio'")

    try:
        with rasterio.open(tiff_path) as src:
            # Read the first band.  Handles single and multi-band TIFFs.
            band1 = src.read(1)
            transform = src.transform
            height, width = band1.shape

            # Create lists to store geometries and values
            geometries = []
            values = []

            for row in range(height):
                for col in range(width):
                    # Get the pixel value
                    value = band1[row, col]

                    #nodata handling
                    if value == src.nodata:
                      continue

                    # Calculate the polygon coordinates using the affine transform
                    x_min, y_max = transform * (col, row)
                    x_max, y_min = transform * (col + 1, row + 1)

                    # Create a polygon for the pixel
                    polygon = Polygon([(x_min, y_max), (x_max, y_max), (x_max, y_min), (x_min, y_min)])
                    geometries.append(polygon)
                    values.append(value)


            # Create the GeoDataFrame
            gdf = gpd.GeoDataFrame({'value': values, 'geometry': geometries}, crs=crs)
            return gdf

    except Exception as e:
        print(f"Error opening or processing TIFF file: {e}")
        return None

