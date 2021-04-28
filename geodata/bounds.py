"""Bounds utils."""

# PyPI modules
import geopandas as gp
from shapely.geometry import Polygon

def convert_bounds_to_gdf(bounds, crs=None):
    """Convert the input to a GeoDataFrame.

        Paramters
        ---------
        bounds : anything that can be interpreted as a bounding polygon
            - lonmin, latmin, lonmax, latmax
            - a Polygon
            - a GeoSeries
            - a GeoDataFrame
        crs : coordinate reference system, optional
            regular lonlat will be assumed as default


        Returns
        -------
        bounds : GeoDataFrame
    """
    # GeoDataFrame -> nothing to do
    # ------------------------------------------------
    if isinstance(bounds, gp.GeoDataFrame):
        return bounds

    # GeoSeries -> GeoDataFrame
    # ------------------------------------------------
    if isinstance(bounds, gp.GeoSeries):
        if crs is None:
            crs = 'epsg:4326'

        gdf = gp.GeoDataFrame(geometry=bounds)
        gdf.crs = crs

        return gdf

    # Polygon -> GeoSeries
    # ------------------------------------------------
    if isinstance(bounds, Polygon):
        geoseries = gp.GeoSeries(bounds)
        return convert_bounds_to_gdf(geoseries)

    # limits -> Polygon
    # ------------------------------------------------
    if len(bounds) != 4:
        raise IndexError('bounds must be of length for, got %i' % len(bounds))

    xmin, ymin, xmax, ymax = bounds

    if not (xmax >= xmin and ymax >= ymin):
        raise ValueError('xmin, ymin, xmax, ymax : %s' % str(bounds))

    coordinates = (
            (xmin, ymin),
            (xmax, ymin),
            (xmax, ymax),
            (xmin, ymax),
            )
    polygon = Polygon(coordinates)
    return convert_bounds_to_gdf(polygon)

def convert_bounds_to_tuple(bounds):
    """Return as (xmin, ymin, xmax, ymax)."""
    if isinstance(bounds, gp.GeoDataFrame):
        return bounds.unary_union.bounds
        # original code replaced on 2020-10-20:
        # return tuple(bounds.bounds.iloc[0])

    if isinstance(bounds, gp.GeoSeries):
        return bounds.unary_union.bounds
        # original code replaced on 2020-10-20:
        # return tuple(bounds.bounds.iloc[0])

    bounds = tuple(bounds)
    assert len(bounds) == 4
    assert bounds[2] >= bounds[0]
    assert bounds[3] >= bounds[1]

    return bounds
