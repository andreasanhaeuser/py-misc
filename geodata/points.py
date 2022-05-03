"""Points utils."""

# PyPI modules
from shapely.geometry import Point
import fiona
import geopandas as gp

def convert_to_gdf(point, crs=None):
    """Try to convert the input to a GeoDataFrame.

        Paramters
        ---------
        point : anything that can be interpreted as a Point
            - a tuple (x, y)
            - a Point
            - a GeoSeries
            - a GeoDataFrame
        crs : coordinate reference system, optional
            regular lonlat will be assumed as default
            
        Returns
        -------
        point : GeoDataFrame
    """
    # GeoDataFrame -> nothing to do
    # ------------------------------------------------
    if isinstance(point, gp.GeoDataFrame):
        return point

    # GeoSeries -> GeoDataFrame
    # ------------------------------------------------
    if isinstance(point, gp.GeoSeries):
        if crs is None:
            crs = fiona.crs.from_epsg(4326) 

        gdf = gp.GeoDataFrame(geometry=point)
        gdf.crs = crs

        return gdf

    # Point -> GeoSeries
    # ------------------------------------------------
    if isinstance(point, Point):
        geoseries = gp.GeoSeries(point)
        return convert_to_gdf(geoseries)

    # tuple -> Polygon
    # ------------------------------------------------
    x, y = point
    point = Point(x, y)
    return convert_to_gdf(point)

def convert_to_tuple(point):
    if isinstance(point, gp.GeoDataFrame):
        return convert_to_tuple(point.geometry)

    if isinstance(point, gp.GeoSeries):
        x = point.x.iloc[0]
        y = point.y.iloc[0]
        return (x, y)

    assert len(point) == 2
    return tuple(point)
