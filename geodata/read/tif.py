"""Tif file data reader."""

# standard modules
import os
import affine
from copy import deepcopy as copy
import warnings
from collections import Iterable

# PyPI modules
import numpy as np
from shapely.geometry import mapping
from shapely.geometry import Point
import geopandas as gp
import rasterio as rio
from rasterio.mask import mask
from pyproj import Proj

# misc
from misc.geodata.bounds import convert_bounds_to_gdf
from misc.math import running_mean

# local
from . import tiled_tif

_data_key = 'data'

def read(
        filename, bounds=None, output_format='array', nan_thresh=None,
        scale_factor=1, layer=0, data_key=None, Nx=None, Ny=None,
        crs=None,
        ):
    """Read tif file.

        Parameters
        ----------
        filename : str
        bounds : GeoDataFrame
        output_format : {'array', 'gdf'}
        nan_thresh : float, optional
        scale_factor : float, optional
        data_key : None or str, optional


        Returns
        -------
        if output_format == 'array'
            data : array
            meta : dict
                meta information
        else
            GeoDataFrame
    """
    if output_format == 'array':
        return read_as_array(
                filename, bounds, nan_thresh, scale_factor=scale_factor,
                Nx=Nx, Ny=Ny, layer=layer, crs=crs,
                )
    elif output_format == 'gdf':
        return read_as_gdf(
                filename, bounds, nan_thresh, scale_factor=scale_factor,
                data_key=data_key, Nx=Nx, Ny=Ny, layer=layer, crs=crs,
                )
    else:
        raise ValueError('Unknown output format: %s' % output_format)

def read_as_gdf(
        filename, bounds=None, nan_thresh=None, data_key=None, scale_factor=1,
        Nx=None, Ny=None, layer=0, crs=None,
        ):
    """Return as GeoDataFrame."""
    if data_key is None:
        data_key = _data_key

    data, meta = read_as_array(
            filename, bounds, nan_thresh, scale_factor, Nx=Nx, Ny=Ny,
            layer=layer, crs=crs,
            )

    # get coordinates
    # =====================================================
    x = meta['x_native']
    y = meta['y_native']

    # build GeoDataFrame
    # =====================================================
    gdf = gp.GeoDataFrame()
    gdf.crs = meta['crs']
    gdf[data_key] = data[0].flatten()

    # points
    xm, ym = np.meshgrid(x, y)
    xf = xm.flatten()
    yf = ym.flatten()
    N = len(xf)
    gdf['geometry'] = [Point(xf[n], yf[n]) for n in range(N)]
    # =====================================================

    return gdf

def read_as_array(
        filename, bounds=None, nan_thresh=None, scale_factor=1, Nx=None,
        Ny=None, layer=0, crs=None,
        ):
    """Return tif data as array and meta as dict.

        Parameters
        ----------
        filename : str or iterable of such
        bounds : GeoDataFrame
        nan_thresh : float

        Returns
        -------
        data : array
        meta : dict
            meta information
    """
    # bounds
    # ---------------------------------------
    if bounds is None:
        bounds = get_global_bounds()
    bounds = convert_bounds_to_gdf(bounds)
    # ---------------------------------------

    # load raw
    # ---------------------------------------
    data, meta = load_file_raw(filename, bounds, layer)

    # handle invalid numbers
    # =====================================================
    if nan_thresh is not None:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            invalid = np.abs(data) > nan_thresh
        data[invalid] = np.nan

    # scale values
    # ---------------------------------------
    if scale_factor != 1:
        data = np.array(data, dtype=float)
        data *= scale_factor

    # reduce resolution
    # ---------------------------------------
    data, meta = reduce_resolution(data, meta, Nx, Ny)

    # crs
    # ---------------------------------------
    meta = change_crs(meta, crs)

    return data, meta

################################################################
# helpers - load                                               #
################################################################
def load_file_raw(filename, bounds, layer=0):
    """Return data as it is in the file."""
    # tiled data
    # =================================================
    if not isinstance(filename, str):
        if not isinstance(filename, Iterable):
            raise TypeError(
                    'filename must be str or iterable of such, got %s'
                    % str(type(filename)),
                    )

        return tiled_tif.read_tiled_tif(filename, bounds)
    # =================================================

    # filename ========================================
    if not os.path.isfile(filename):
        raise OSError('File does not exist: %s' % filename)

    # load data ===========================================
    with rio.open(filename, 'r') as fid:
        crs_file = fid.crs
        try:
            polygon = bounds.to_crs(crs_file).geometry.iloc[0]
        except Exception as exception:
            if 'init' in crs_file:
                crs_file = crs_file['init']
            polygon = bounds.to_crs(crs_file).geometry.iloc[0]
        except Exception as exception:
            raise exception

        extent = mapping(polygon)
        data, transform = mask(fid, shapes=[extent], crop=True)
        meta = fid.meta.copy()

    data = data[layer]

    meta.update({
            'driver' : 'GTiff',
            'height' : data.shape[0],
            'width' : data.shape[1],
            'transform' : transform,
            }
        )

    # raw coordinates
    # =====================================================
    x, y = compute_coordinates_from_transform(meta)
    meta['x_native'] = x
    meta['y_native'] = y

    return data, meta

################################################################
# helpers - bounds                                             #
################################################################
def get_global_bounds():
    xmin = -179.999999999999999999
    xmax = 180.0
    ymin = -90.0
    ymax = 90.0
    bounds = (xmin, ymin, xmax, ymax)
    return bounds

################################################################
# helpers - misc                                               #
################################################################
def change_crs(meta, crs=None, inplace=True):
    if not inplace:
        meta = copy(meta)

    if crs is None:
        return meta

    lon = meta['x_native']
    lat = meta['y_native']
    mlat, mlon = np.meshgrid(lat, lon)

    # convert to crs
    # =====================================================
    proj = Proj(crs)
    x, y = proj(mlon, mlat)

    # manually convert m -> km
    if 'units=km' in crs:
        x /= 1000.
        y /= 1000.
    # =====================================================

    meta['x'] = x
    meta['y'] = y
    meta['crs'] = crs

    return meta

def reduce_resolution(data, meta, Nx=None, Ny=None):
    shape = np.shape(data)
    size_y, size_x = shape

    # determine increment
    # =====================================================
    inc_x = 1
    if Nx is not None:
        inc_x = int(size_x // Nx)
        inc_x = max(inc_x, 1)

    inc_y = 1
    if Ny is not None:
        inc_y = int(size_y // Ny)
        inc_y = max(inc_y, 1)
    # =====================================================

    # reduce data
    # =====================================================
    # data = data[::inc_y, ::inc_x]
    rm = running_mean.running_mean_equidistant_no_weights

    if inc_y > 1:
        hwy = inc_y //2
        mean = rm(data, hwy, hwy, axis=0)
        data = mean[::inc_y]
    if inc_x > 1:
        hwx = inc_x //2
        mean = rm(data, hwx, hwx, axis=1)
        data = mean[:, ::inc_x]

    meta['x_native'] = meta['x_native'][::inc_x]
    meta['y_native'] = meta['y_native'][::inc_y]
    # =====================================================

    # adjust transform
    # =====================================================
    tr_values = np.array(meta['transform'].column_vectors)
    tr_values[0, 0] *= inc_x
    tr_values[1, 1] *= inc_y
    meta['transform'] = affine.Affine(*tr_values.T.flatten())
    # =====================================================

    # adjust other meta
    # =====================================================
    meta['height'] = np.shape(data)[0]
    meta['width'] = np.shape(data)[1]
    # =====================================================

    return data, meta

def compute_coordinates_from_transform(meta):
    Nx = meta['width']
    Ny = meta['height']
    transform = meta['transform']

    # i : indices in x-direction
    # j : inideces in y-direction
    i = np.arange(Nx)
    j = np.arange(Ny)

    # apply affine transformation from index to coordinate
    # `x0` and `y0` remain unused
    # (the zeros can be replaced by any value without altering the output)
    x, y0 = transform * (i, 0)
    x0, y = transform * (0, j)
    return x, y
