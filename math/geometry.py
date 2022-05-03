#!/usr/bin/python3
"""Geometry functions.

    Author
    ------
    Andreas Anhaeuser (AA)
    Insitute for Geophysics and Meteorology
    University of Cologne, Germany
    <andreas.anhaeuser@posteo.net>
"""

# standard modules
import collections
import warnings 

# PyPI modules
import numpy as np
from shapely.geometry import Polygon

# The spheroid values are from the World Geodetic System:
_earth_radius       = 6.3710e6    # (m) Earth's mean radius
_earth_polar_radius = 6.356752e6  # (m) Earth's polar radius
_earth_eq_radius    = 6.378137e6  # (m) Earth's equatorial radius

def distance_on_sphere(
        longitude_1,
        latitude_1,
        longitude_2,
        latitude_2,
        units='deg',
        radius=_earth_radius,
        ):
    """Return distance between points on a speherical surface.
    
        Parameters
        ----------
        longitude_1 : float or array
            arbitrary number of axes
        latitude_1 : float or array
            same shape as longtitude_1
        longitude_2 : float or array
            arbitrary number of axes
        latitude_2 : float or array
            same shape as longtitude_2
        units : str, optional
            'deg' or 'rad', refers to the first four parameters, default is
            'deg'.
        radius : float or array, optional
            radius of the sphere. Default: Earth radius in metres. 
            If `radius` is not scalar its dimensions must be compatible for
            multiplication with an array of shape (shape1 + shape2), where
            shape1 is the shape of `longitude_1` and shape2 is the shape of
            `longitude_2`.

        Returns
        -------
        distance : float or array
            the shape is (shape1 + shape2), where shape1 is the shape of
            `longitude_1` and shape2 is the shape of `longitude_2`.
        
        History
        -------
        2021-11-17: (AA) str->float conversion
        2018-12-01: (AA) Extention to arrays of arbitrary shape
        2014-2016 : (AA)
    """
    ############################################################
    # cast str -> float                                        #
    ############################################################
    if isinstance(longitude_1, str):
        longitude_1= float(longitude_1)
    if isinstance(latitude_1, str):
        latitude_1= float(latitude_1)
    if isinstance(longitude_2, str):
        longitude_2= float(longitude_2)
    if isinstance(latitude_2, str):
        latitude_2= float(latitude_2)

    ###################################
    # INPUT CHECK                     #
    ###################################
    # ============ units  ================================ #
    known_units = ('deg', 'degrees', 'rad', 'radians')
    message = "units must be 'deg' oder 'rad'."

    if not isinstance(units, str):
        raise TypeError(message)

    if units.lower() not in known_units:
        raise ValueError(message)

    # ========== shapes  ================================= #
    shape_1 = np.shape(longitude_1)
    shape_2 = np.shape(longitude_2)

    if np.shape(latitude_1) != shape_1:
        raise IndexError('Shapes of longitude_1 and latitude_1 must agree.')

    if np.shape(latitude_2) != shape_2:
        raise IndexError('Shapes of longitude_2 and latitude_2 must agree.')

    ###################################
    # CONVERSIONS                     #
    ###################################
    if units.lower()[:1] == 'd':
        # deg --> rad
        lon1 = np.radians(longitude_1)
        lat1 = np.radians(latitude_1)
        lon2 = np.radians(longitude_2)
        lat2 = np.radians(latitude_2)
    else:
        lon1 = longitude_1
        lat1 = latitude_1
        lon2 = longitude_2
        lat2 = latitude_2

    ###################################################
    # ALGEBRAICALLY EXACT FUNCTION                    #
    ###################################################
    dlon = np.add.outer(-lon1, lon2)
    summand1 = np.multiply.outer(np.sin(lat1), np.sin(lat2))
    summand2 = np.multiply.outer(np.cos(lat1), np.cos(lat2)) * np.cos(dlon)
    S = summand1 + summand2

    # ========== correct for numerical errors  =========== #
    # Numerical errors can lead to S > 1.
    # This is mathematically not possible.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        too_large = S > 1.

    if np.sum(too_large) > 0:
        if not isinstance(too_large, np.ndarray):
             S = 1.
        else:
             S[too_large] = 1.
    # ==================================================== #

    # central angle (separation angle between the points)
    ds = np.arccos(S)

    ###################################################
    # Handle small angular separation                 #
    ###################################################
    # For small central angles, the above algeraically exact function causes
    # numerical errors.
    # --> Use better suited haversine function for small angles.
        
    # check for small angles
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        ds_is_small = ds < 1e-3

    # ========== correct small angles  =================== #
    if np.sum(ds_is_small) > 0:
        # compute the haversine small angle
        dlat = np.add.outer(-lat1, lat2)
        summand1 = (np.sin(dlat / 2.))**2
        prod = np.multiply.outer(np.cos(lat1), np.cos(lat2))
        summand2 = prod * (np.sin(dlon / 2.))**2
        ds_hav = 2 * np.arcsin(np.sqrt(summand1 + summand2))

        if not isinstance(too_large, np.ndarray):
            ds = ds_hav
        else:
            ds[ds_is_small] = ds_hav[ds_is_small]
    # ==================================================== #
                
    distance = radius * ds
    return distance

def nearest_neighbour_on_sphere(
        lon_point,
        lat_point,
        lon_array,
        lat_array,
        units='deg',
        lonlatmode='mesh',
        values=None,
        radius=_earth_radius,
        skip_nans=False,
        ):
    """Return the nearest neighbour on a sphere.
    
        lon_array and lat_array can be given in different formats, which is
        specified by lonlatmode:
        1) If you give them as a field (meshgrid), where lon[i, j] and 
           lat[i, j] are the coordinates at point [i, j], set lonlatmode to
           'mesh'. In this case, index is the pair (i, j). This works also with
           arrays that have more than 2 axes for some reason. Then index is
           returned as (i, j, k, ...)
        2) If they are a 1d-list of equal length, where lon[i] and lat[i] give
           the coordinates of point i, set lonlatmode to 'list'. In this case,
           index is the tuple (i,)
        3) If they are axes, so that the coordinates of point [i, j] are given
           by lon[i] and lat[j], then set lonlatmode to 'axes'. In this case,
           index is the pair (i, j)

        Parameters
        ----------
        lon_point : float
            longitude of the central point
        lat_point : float
            latitude of the central point
        lon_array : array of floats
            array with arbitrary number of axes
        lat_array : array of floats
            array with arbitrary number of axes
            if lonlatmode is not 'axes', then dimensions must be same as of
            lon_array
        units : 'deg' or 'rad'
            units of longitudes and latitudes
        lonlatmode : 'mesh', 'list' or 'axes'
            determines how lon_array and lat_array are interpreted.
        values : None or array or list, optional
            values at the grid points. Same length as lon_array and lat_array.
            array with arbitrary number of axes. Dimensions must be same as of
            lon_array. Default: None
        skip_nans : bool, optional
            if True, then positions in the grid are ignored if the
            corresponding element in values is nan. Default: False.
      
        Returns
        -------
        nlon :float or None
            longitude of the nearest point. None if no such point exists
        nlat : float or None
            latitude of the nearest point. None if no such point exists
        idx : tuple of ints or None
            index of nlon and nlat in lon_array and lat_array
        value : object or None
            the values[idx] (if applicable)
        dist : float
            distance of (nlon, nlat) to (lon_point, lat_point)

        History
        -------
        2014-2016 : (AA)
    """ 
    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    # units:
    if not isinstance(units, str):
        raise TypeError("units must be 'deg' oder 'rad'.")
    if units.lower() not in ['deg', 'degrees', 'rad', 'radians']:
        raise ValueError("units must be 'deg' oder 'rad'.")

    # lonlatmode:
    if lonlatmode.lower() in ['m', 'mesh', 'field', '2d', 'list', 'l']:
        if not np.shape(lon_array) == np.shape(lat_array):
            raise TypeError('input arrays must be of same shape.')
        lon, lat = lon_array, lat_array
    elif lonlatmode.lower() in ['axes', 'a']:
        lon, lat = np.meshgrid(lon_array, lat_array)
    else:
        raise ValueError('unsupported value for lonlatmode.')

    ###################################################
    # TRIVIAL CASE                                    #
    ###################################################
    if not np.any(lon):
        return (None, None, None, None, np.nan)

    ###################################################
    # FIND NEAREST NEIGHBOUR                          #
    ###################################################
    D = distance_on_sphere(
            lon_point, lat_point, lon, lat, units=units,
            radius=radius)

    # exclude nan's:
    D[np.isnan(D)] = np.inf

    # exclude positions where values is nan:
    if skip_nans and np.any(values):
        D[np.isnan(values)] = np.inf

    # find nearest:
    idx = np.unravel_index(D.argmin(), D.shape)

    # check whether it is not infinity (this is the case if there is no
    # suitable point):
    if not D[idx] < np.inf:
        idx = None

    ###################################################
    # RETURN                                          #
    ###################################################
    # special case: no nearest neighbour found
    if idx is None:
        return (np.nan, np.nan, None, None, np.nan)

    # regular case
    if values is None:
        return (lon[idx], lat[idx], idx, None, D[idx])
    else:
        return (lon[idx], lat[idx], idx, values[idx], D[idx])        

def radius_on_spheroid(
        latitude, radius_e=_earth_eq_radius,
        radius_p=_earth_polar_radius):
    """Return the local radius on an oblate ellipsoid of revolution.

        Parameters
        ----------
        latitude : number (deg)
                   geographical latitude of the point.
        radius_e : number > 0
                   equatorial radius of the spheroid.
                   Default: Value for the Earth in meters.
        radius_p : number > 0
                   polar radius of the spheroid.
                   Default: Value for the Earth in meters.

        Returns
        -------
        float
            Distance of the point to the center of the spheroid.
            
        Note
        ----
        This is not equivalent to the radius of curvature at this point.

        History
        -------
        2016 : (AA)
    """

    # convert from deg to rad:
    phi = latitude * np.pi/180.

    # abbreviations:
    a = radius_e
    b = radius_p

    numerator   = (a**2 * np.cos(phi))**2 + (b**2 * np.sin(phi))**2
    denominator = (a    * np.cos(phi))**2 + (b    * np.sin(phi))**2
    return np.sqrt(numerator / denominator)

def circular_polygon(x, y, radius, N=1024):
    """Return a Polygon whose corners lie on a circle.

        Parameters
        ----------
        x : float
            x-coordinate of the center
        y : float
            x-coordinate of the center
        radius : float
        N : int (positive)
            number of corners

        Returns
        -------
        shapely.Polygon
    """
    N = int(N)
    if N <= 0:
        raise ValueError('N must be positve, got %s' % N)

    phi_inc = 2 * np.pi / N
    phi_min = 0.
    phi_max = 2 * np.pi
    phis = np.arange(phi_min, phi_max, phi_inc)

    x_edge = x + radius * np.cos(phis)
    y_edge = y + radius * np.sin(phis)

    points = zip(x_edge, y_edge)
    polygon = Polygon(points)
    return polygon


###################################################
# TESTING                                         #
###################################################
if __name__ == '__main__':
    N1 = 10**3
    N2 = 10**4
    lons1 = np.arange(1., 100., 0.0001)[:N1]
    lats1 = np.arange(50., 60., 0.0001)[:N1]
    lons2 = np.arange(1., 100., 0.001)[:N2]
    lats2 = np.arange(50., 60., 0.001)[:N2]
    lon1 = lons1[0]
    lat1 = lats1[0]
    lon2 = lons2[-1]
    lat2 = lats2[-1]

    lons1_2d = np.random.random((2, 3))
    lats1_2d = np.random.random((2, 3))
    lons2_2d = np.random.random((4, 5))
    lats2_2d = np.random.random((4, 5))

    print('\n0d vs. 0d')
    print(distance_on_sphere_old(lon1, lat1, lon2, lat2))
    print(distance_on_sphere(lon1, lat1, lon2, lat2))

    print('\n0d vs. 1d')
    print(distance_on_sphere_old(lon1, lat1, lons2, lats2))
    print(distance_on_sphere(lon1, lat1, lons2, lats2))

    fold = lambda : distance_on_sphere_old(lon1, lat1, lons2, lats2)
    fnew = lambda : distance_on_sphere(lon1, lat1, lons2, lats2)
    test = np.prod(fold() == fnew())

    print('\n1d vs. 1d')
    print(distance_on_sphere(lons1, lats1, lons2, lats2))

    print('\n2d vs. 2d')
    f = lambda : distance_on_sphere(lons1_2d, lats1_2d, lons2_2d, lats2_2d)
    D = f()
    print(np.shape(D))
    print(D)
