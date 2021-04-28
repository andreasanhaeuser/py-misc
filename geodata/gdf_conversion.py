"""Convert from and to GeoDataFrame."""

import numpy as np

def convert_point_gdf_to_array(gdf, value_name):
    x = np.array(gdf.geometry.x)
    y = np.array(gdf.geometry.y)
    values = np.array(gdf[value_name])

    return x, y, values
