"""shp file data reader."""

# standard modules
import os

# PyPI modules
import numpy as np
import geopandas as gp
import fiona
import pyproj
import shapely

# misc
from misc.geodata.bounds import convert_bounds_to_gdf

def read_as_gdf(filename, bounds=None, crs=None):
    """Return shp file as GeoDataFrame.

        Parameters
        ----------
        filename : str
        bounds : None or tuple
        crs : None or CRS

        Returns
        -------
        GeoDataFrame
    """
    if not os.path.isfile(filename):
        raise OSError('File does not exist: %s' % filename)

    # initialize
    shapes = []
    with fiona.open(filename) as fid:
        string = fid.crs['init']
        crs_file = pyproj.CRS.from_string(string)

        if bounds is not None:
            bounds = convert_bounds_to_gdf(bounds)
            bounds = bounds.to_crs(crs_file)
            poly_bounds = bounds.unary_union.convex_hull

        for shape in fid:
            # apply bounds filter
            if bounds is not None:
                coordinates = shape['geometry']['coordinates']
                if not polygons_intersect(poly_bounds, coordinates):
                    continue

            # the object is good to use if this line is reached
            shapes.append(shape)

    if any(shapes):
        gdf = gp.GeoDataFrame.from_features(shapes)
    else:
        gdf = gp.GeoDataFrame(geometry=[])
    gdf.crs = crs_file

    # convert crs
    if crs is not None:
        gdf = gdf.to_crs(crs)
        gdf.crs = crs

    return gdf

################################################################
# helpers                                                      #
################################################################
def polygons_intersect(polygon, coordinates):
    """Return a bool."""
    shape = np.shape(coordinates)

    # TODO: The above sometimes throws an error.
    # Use the following few lines for debugging
#    import warnings
#    with warnings.catch_warnings():
#        warnings.simplefilter('error', Warning)
#        try:
#            shape = np.shape(coordinates)
#        except Exception:
#            1/0

    if not any(shape):
        return False

    if len(shape) != 2:
        recursive = True
    elif shape[-1] != 2:
        recursive = True
    elif shape[0] < 3:
        recursive = True
    else:
        recursive = False

    if recursive:
        intersect = False
        for c in coordinates:
            intersect = polygons_intersect(polygon, c)
            if intersect:
                break
        return intersect

    poly_other = shapely.geometry.Polygon(coordinates)

    return polygon.intersects(poly_other)
