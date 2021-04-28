"""Helpers to coordinate reference systems (CRS)."""

# standard
import warnings

# PyPI
import numpy as np
from pyproj import Proj
from pyproj import CRS

# misc
from misc.physics import units as unit_utils
from misc.text.string_utils import lcrop

def lonlat_to_utm(lon, lat, crs=None, units='m', false_northing=True):
    """Return x, y."""
    if crs is None:
        proj4 = get_utm_proj4string(lon, lat, units)
        crs = CRS.from_proj4(proj4)

    crs = cast_to_crs(crs)
    proj = cast_to_proj(crs)

    # apply projection
    x, y = proj(lon, lat)

    # remove false northing
    units = retrieve_units_from_crs(crs)
    if false_northing:
        pass
    elif not hasattr(crs, 'utm_zone'):
        pass
    elif not crs.utm_zone.endwith('S'):
        pass
    elif units == 'm':
        y -= 10**7
    elif units == 'km':
        y -= 10**4
    else:
        raise ValueError('Unknown units: %s' % get_utm_units(crs))

    return x, y

def utm_to_lonlat(x, y, crs):
    """Return lon, lat."""
    proj = cast_to_proj(crs)
    return proj(x, y, inverse=True)

def get_utm_projection(lon, lat, units='m', fractional=False):
    proj4string = get_utm_proj4string(lon, lat, units, fractional)
    return Proj(proj4string)

def get_utm_proj4string(lon, lat, units='m', fractional=False):
    """Return a proj4sring."""
    zone, hem = get_utm_zone(lon, lat, fractional)
    return utm_proj4string_from_zone(zone, hem, units)

def utm_proj4string_from_zone(zone, hem, units='m'):
    """Return a proj4sring.

        Parameters
        ----------
        zone : int
        hem : 'N' or 'S'

        Returns
        -------
        proj4string : str
    """
    if hem.lower() in ('s', 'south'):
        hem_long = 'south'
    elif hem.lower() in ('n', 'north'):
        hem_long = 'north'
    else:
        raise ValueError('Unkown UTM hemisphere: %s' % hem)

    fmt = '+proj=utm +datum=WGS84 +no_defs +zone=%s +%s +units=%s'
    proj4string = fmt % (str(zone), hem_long, units)
    return proj4string


################################################################
# helpers                                                      #
################################################################
def cast_to_proj(x):
    if isinstance(x, Proj):
        return x
    return Proj(x)

def cast_to_crs(x):
    if isinstance(x, CRS):
        return x
    elif isinstance(x, str):
        return CRS.from_proj4(x)
    raise TypeError('Cannot handle type %s' % type(x))

def get_utm_zone(lon, lat, fractional=False):
    """Return UTM zone as int an and str.

        Does not take into account any exceptions in the UTM zone definitions
        (neither those around Norway, nor those in the polar regions).

        Parameters
        ----------
        lon : float
            (deg) geographical longitude
        lat : float
            (deg) geographical latitute

        Returns
        -------
        zone : int
            a value from 1..60
        hem : 'N' or 'S'
    """
    assert -180 <= lat <= 360
    assert -90 <= lat <= 90

    # compute UTM zone
    # ------------------------------------------------
    # fraction of plant's circumference into eastward direction,
    # starting at -180W:
    fraction = ((lon + 180) / 360) % 1          # 0..1
    zone_fractional = fraction * 60 + 0.5       # 0.5..60.5
    if fractional:
        zone = zone_fractional
    else:
        zone = int(np.round(zone_fractional))

        # Check if in bounds
        warn_fmt = 'Computed zone %i for lon %f. Setting to zone %i'
        if zone == 0:
            newzone = 60
            warnings.warn(warn_fmt % (zone, lon, newzone))
            zone = newzone
        if zone == 61:
            warnings.warn(warn_fmt % (zone, lon, newzone))
            zone = 1
    # ------------------------------------------------

    # compute UTM hemisphere
    # ------------------------------------------------
    if lat > 0:
        hem = 'N'
    else:
        hem = 'S'
    # ------------------------------------------------

    return zone, hem

def retrieve_units_from_proj4(s):
    """Return a str."""
    pattern = 'units='
    L = len(pattern)
    if pattern not in s:
        return 'm'

    ibeg = s.index(pattern) + L
    words = s[ibeg:].split()
    units = words[0]
    return units

def retrieve_units_from_crs(crs):
    pattern = '+units='
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        proj4 = crs.to_proj4()
    return retrieve_units_from_proj4(proj4)

################################################################
# deprecated                                                   #
################################################################
def lonlat_to_utm_old(lon, lat, crs=None, units='m'):
    """Return x, y."""
    raise Exception('Deprecated.')
    if crs is None:
        crs = get_utm_proj4string(lon, lat, units)

    # create projection
    if isinstance(crs, Proj):
        proj = crs
    else:
        proj = Proj(crs)

    # apply projection
    x, y = proj(lon, lat)

    # scale
    units = retrieve_units_from_crs(crs)
    scale_factor = unit_utils.get_scale_factor(units)
    x_out = x / scale_factor
    y_out = y / scale_factor

    return x_out, y_out

def utm_to_lonlat_old(x, y, crs):
    raise Exception('Deprecated.')
    """Return lon, lat."""
    crs = cast_to_crs(crs)
    proj = cast_to_proj(crs)

    # scale to SI base units
    units = retrieve_units_from_crs(crs)
    scale_factor = unit_utils.get_scale_factor(units)
    x_si = x * scale_factor
    y_si = y * scale_factor

    return proj(x_si, y_si, inverse=True)

