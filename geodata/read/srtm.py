#!/usr/bin/python3

# standard modules
import os
import sys
from copy import deepcopy as copy
import itertools
import warnings

# PyPI modules
import numpy as np
import elevation

# IGMK modules
from misc.chronometer import Chronometer
from misc.geodata.bounds import convert_bounds_to_tuple

if __name__ != '__main__':
    from . import tif


# constants
# ========================================================
if 'DIRNAME_GAPU_DATA' in os.environ:
    _dirname_cache = os.environ['DIRNAME_GAPU_DATA'] + '/geodata/srtm'
else:
    _dirname_cache = os.path.expanduser('~/.cache/srtm')
_filename_tmp_fmt = _dirname_cache + '/cache.%i.tif'

_product = 'SRTM1'
_RE = 6371e3
_height_factor = 2
_sun_zenith_angle = 30.
_sun_azimuth_angle = 240.
# ========================================================

def get_elevation(
        bounds, filename_tmp=None, remove_tmp_file=True, reuse_tmp_file=False,
        dirname_cache=None, product=None, margin='1%', Nx=4096, Ny=4096,
        crs=None, max_download_tiles=50, 
        ):
    """

        Parameters
        ----------
        bounds : iterable of float
            (deg) lonmin, latmin, lonmax, latmax
        filename_tmp : str, optional
            temporary file
        dirname_cache : str, optional
            directory where SRTM tiles are to be stored
        remove_tmp_file : bool, optional
            (default: True)
        margin : str, optional

        Returns
        -------
        data : 2d-array
            elevation
        meta : dict
            coordinates and such
    """
    # set defaults
    # =====================================================
    if filename_tmp is None:
        filename_tmp = get_filename_tmp()
    if dirname_cache is None:
        dirname_cache = _dirname_cache

    filename_tmp = os.path.abspath(filename_tmp)
    # =====================================================

    # expand user
    # =====================================================
    filename_tmp = os.path.expanduser(filename_tmp)
    dirname_cache = os.path.expanduser(dirname_cache)
    # =====================================================

    # create directories
    # =====================================================
    dirname_tmp = os.path.dirname(filename_tmp)
    for dirname in (dirname_cache, dirname_tmp):
        if dirname == '':
            continue
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
    # =====================================================

    # load
    # =====================================================
    bounds = convert_bounds_to_tuple(bounds)

    # clip
    if not reuse_tmp_file or not os.path.isfile(filename_tmp):
        # set cache dir
        elevation.CACHE_DIR = dirname_cache
        elevation.DEFAULT_PRODUCT = product

        # download if necessary and crop section to file
        kwargs = {
                'output' : filename_tmp,
                'margin' : margin,
                'max_download_tiles' : max_download_tiles,
                }
        elevation.clip(bounds, **kwargs)

    # load file
    data, meta = tif.read(filename_tmp, bounds=bounds, Nx=Nx, Ny=Ny, crs=crs)

    # remove temporary file
    # =====================================================
    if remove_tmp_file and os.path.isfile(filename_tmp):
        os.remove(filename_tmp)
    # =====================================================

    # crop singleton first dimension
    # =====================================================
    while np.shape(data)[0] == 1 and len(np.shape(data)) > 2:
        data = data[0]
    # =====================================================

    # convert to float and replace nans
    # =====================================================
    # int -> float
    data = 1. * data

    # NaN
    idx_nan = data == meta['nodata']
    data[idx_nan] = np.nan
    # =====================================================

    return data, meta

def get_brightness(bounds, kwargs_elevation=None, kwargs_brightness=None):
    if kwargs_elevation is None:
        kwargs_elevation = {}

    if kwargs_brightness is None:
        kwargs_brightness = {}

    data, meta = get_elevation(bounds, **kwargs_elevation)
    brightness = compute_brightness(data, meta, **kwargs_brightness) 
    return brightness, meta

def compute_brightness(
        data, meta, height_factor=None, compute_shadows=False,
        sun_zenith_angle=None, sun_azimuth_angle=None,
        ):
    """Cast shadows from higher tips to lower areas.

        Paramters
        ---------
        sun_zenith_angle
            (deg)
        sun_azimuth_angle
            (deg)
    """
    # defaults ============================================
    if height_factor is None:
        height_factor = _height_factor
    if sun_zenith_angle is None:
        sun_zenith_angle = _sun_zenith_angle
    if sun_azimuth_angle is None:
        sun_azimuth_angle = _sun_azimuth_angle

    # recursive function call
    # =====================================================
    if compute_shadows:
        return compute_brightness_with_shadows(
                data, meta, height_factor=height_factor,
                sun_zenith_angle=sun_zenith_angle,
                sun_azimuth_angle=sun_azimuth_angle,
                )

    # retrieve variables
    # =====================================================
    x, y, z = get_xyz(data, meta, height_factor)
    sun_x, sun_y, sun_z = get_sun_vector(sun_zenith_angle, sun_azimuth_angle)

    # finite differences
    # =====================================================
    diff_yj = - np.gradient(y)[0]
    diff_xi = np.gradient(x)[1]

    diff_z = np.gradient(z)
    diff_zj = - diff_z[0]
    diff_zi = diff_z[1]

    # vector perpendicular to the surface (points upwards)
    # =======================================================
    perp_x = - diff_yj * diff_zi
    perp_y = - diff_xi * diff_zj
    perp_z = + diff_xi * diff_yj

    # make sure the perpendicular points upwards
    if np.sum(perp_z < 0) > np.sum(perp_z >= 0):
        perp_x *= -1
        perp_y *= -1
        perp_z *= -1

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        perp_len = np.sqrt(perp_x**2 + perp_y**2 + perp_z**2)
        cos_perp = (perp_z) / perp_len

    # illumination (without shadows)
    # =======================================================
    # angle between sun and the perpendicular
    cos = (perp_x * sun_x + perp_y * sun_y + perp_z * sun_z) / perp_len

    illumination = cos

    # fix negative values
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        illumination[illumination<0] = 0

    # brightness
    brightness = illumination * perp_z / perp_len
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        brightness[brightness<0] = 0
    # =======================================================

    return brightness

def compute_brightness_with_shadows(
        data, meta, height_factor, sun_zenith_angle, sun_azimuth_angle,
        ):
    """Cast shadows from higher tips to lower areas.

        Paramters
        ---------
        sun_zenith_angle
            (deg)
        sun_azimuth_angle
            (deg)
    """
    brightness = compute_brightness(
            data, meta, height_factor=height_factor,
            sun_zenith_angle=sun_zenith_angle,
            sun_azimuth_angle=sun_azimuth_angle,
            )
    if sun_zenith_angle == 0:
        return brightness

    # retrieve variables
    # =====================================================
    x, y, z = get_xyz(data, meta, height_factor)
    sun_x, sun_y, sun_z = get_sun_vector(sun_zenith_angle, sun_azimuth_angle)
    # =====================================================
    
    shadow_max = copy(z)
    z_min = np.nanmin(z)

    J, I = np.shape(data)
    header = 'Compute shadows'
    count = 0
    with Chronometer(J*I, header=header) as chrono:
        for j_obstacle, i_obstacle in itertools.product(range(J), range(I)):

            # chrono
            # ----------------------
            if count % 100 == 0:
                chrono.show()
            count += 1
            chrono.loop()
            # ----------------------

            # coordinates of the shadow-casting point
            # ---------------------------------------
            z_obstacle = z[j_obstacle, i_obstacle]
            y_obstacle = y[j_obstacle, i_obstacle]
            x_obstacle = x[j_obstacle, i_obstacle]
            # ---------------------------------------

            if np.abs(sun_x) >= np.abs(sun_y):
                # traverse in x-direction

                # ray of light
                # -------------------
                dx = x[j_obstacle, :] - x_obstacle
                dy = dx / sun_x * sun_y
                dz = dx / sun_x * sun_z

                y_ray = y_obstacle + dy
                z_ray = z_obstacle + dz
                # -------------------

                # index of i-positions that need to be checked
                # --------------------------------------------
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)

                    # ray of light must be above ground
                    idx = z_ray > z_min

                    # ray of light must have hit the object (be below it)
                    idx &= dz < 0

                iterator = np.arange(I)[idx]
                # --------------------------------------------

                for i in iterator:
                    z_point = z_ray[i]

                    y_point = y_ray[i]
                    y_choice = y[:, i]

                    jlos = np.where(y_choice >= y_point)[0]
                    jhis = np.where(y_choice <= y_point)[0]

                    # check if ray of light/shadow is above the point
                    j_neighbours = []
                    if np.any(jlos):
                        j_neighbours.append(jlos[-1])
                    if np.any(jhis):
                        j_neighbours.append(jhis[0])

                    for j in j_neighbours:
                        if z_point > shadow_max[j, i]:
                            shadow_max[j, i] = z_point

            else:
                # traverse in y-direction

                # ray of light
                # -----------------------
                dy = y[:, i_obstacle] - y_obstacle
                dx = dy / sun_y * sun_x
                dz = dy / sun_y * sun_z

                x_ray = x_obstacle + dx
                z_ray = z_obstacle + dz
                # -----------------------

                # index of j-positions that need to be checked
                # --------------------------------------------
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)

                    # ray of light must be above ground
                    jdx = z_ray > z_min

                    # ray of light must have hit the object (be below it)
                    jdx &= dz < 0

                iterator = np.arange(J)[jdx]
                # --------------------------------------------

                for j in iterator:
                    z_point = z_ray[j]

                    x_point = x_ray[j]
                    x_choice = x[j, :]

                    ilos = np.where(x_choice >= x_point)[0]
                    ihis = np.where(x_choice <= x_point)[0]

                    # check if ray of light/shadow is above the point
                    i_neighbours = []
                    if np.any(ilos):
                        i_neighbours.append(ilos[0])
                    if np.any(ihis):
                        i_neighbours.append(ihis[-1])

                    for i in i_neighbours:
                        if z_point > shadow_max[j, i]:
                            shadow_max[j, i] = z_point

        idx_dark = shadow_max > z
        brightness[idx_dark] = 0.
        return brightness


################################################################
# helpers                                                      #
################################################################
def get_filename_tmp():
    fmt = _filename_tmp_fmt
    rand = np.random.random() * 10**16
    filename = fmt % rand
    return filename

def get_sun_vector(theta_deg, phi_deg):
    """Return (ex, ey, ez)."""
    theta = np.radians(theta_deg)
    phi = np.radians(phi_deg)

    ex = np.sin(theta) * np.sin(phi)
    ey = np.sin(theta) * np.cos(phi)
    ez = np.cos(theta)
    
    return ex, ey, ez

def get_xyz(data, meta, height_factor):
    lons = np.radians(meta['x_native'])
    lats = np.radians(meta['y_native'])
    mlons, mlats = np.meshgrid(lons, lats)

    x = _RE * mlons * np.cos(mlats)
    y = _RE * mlats
    z = data * height_factor

    return x, y, z


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import tif
    def plot(data, meta, nsub, cmap, sun, vmin=None, vmax=None):
        ax = plt.subplot(1, 1, nsub)
        y = meta['y_native']
        x = meta['x_native'] * np.cos(np.radians(np.nanmean(y)))
        ax.axis('equal')
        plt.pcolormesh(x, y, data, cmap=cmap, vmin=vmin, vmax=vmax)

        xcenter = np.nanmean(x)
        ycenter = np.nanmean(y)
        xspan = np.ptp(x)
        yspan = np.ptp(y)
        sun_x, sun_y, sun_z = sun

        f_sun = 0.05 * xspan
        f_shadow = 2 * f_sun

        xsun = xcenter + sun_x * f_sun
        ysun = ycenter + sun_y * f_sun
        xtip = xcenter - sun_x * f_shadow
        ytip = ycenter - sun_y * f_shadow
        # plt.plot([xsun], [ysun], 'yo', markersize=10)
        # plt.plot([xcenter], [ycenter], 'y.')
        # plt.plot([xcenter, xtip], [ycenter, ytip], 'y-')

        plt.colorbar()

    # Kobe
    lon = 135.190
    lat = 34.725
    dlon = 0.5
    dlat = 0.5
    product = 'SRTM1'

    thetas_sun_deg = [30., 30., 45, 60]
    phis_sun_deg = [30., 240., 240, 240]

    # Bengaluru
    lon = 77.59369
    lat = 12.97194
    dlon = 1.0
    dlat = 1.0
    product = 'SRTM3'
    thetas_sun_deg = [60]
    phis_sun_deg = [240]
    height_factor = 10

    # Reunion
    lon = 55.53
    lat = -21.13
    dlon = 0.8
    dlat = 0.6
    product = 'SRTM1'
    theta_sun = 45
    phi_sun = 240
    height_factor = 2

    # Korsika
    lon = 9.09
    lat = 42.16
    dlon = 1.2
    dlat = 1.8
    product = 'SRTM1'
    theta_sun = 45
    phi_sun = 240
    height_factor = 2

    # Kukenan
    lon = -60.83
    lat = 5.2
    dlon = 0.5
    dlat = 0.5
    product = 'SRTM1'
    theta_sun = 45
    phi_sun = 240
    height_factor = 1

    # Wetterstein
    lon = 11.11
    lat = 47.41
    dlon = 0.5
    dlat = 0.5
    product = 'SRTM1'
    theta_sun = 45
    phi_sun = 240
    height_factor = 1

    # Bike trip
    lon = 14.06
    lat = 46.23
    dlon = 0.8
    dlat = 0.8
    product = 'SRTM1'
    theta_sun = 45
    phi_sun = 240
    height_factor = 1

    # Large section
    lon = 14.06
    lat = 46.23
    dlon = 4.0
    dlat = 3.0
    product = 'SRTM1'
    theta_sun = 45
    phi_sun = 240
    height_factor = 1

    # Triglav
    lon = 13.82
    lat = 46.39
    dlon = 0.3
    dlat = 0.3
    product = 'SRTM1'
    theta_sun = 45
    phi_sun = 240
    height_factor = 1

    lonmin = lon - dlon/2
    lonmax = lon + dlon/2
    latmin = lat - dlat/2
    latmax = lat + dlat/2

    kwargs_ele = {'product' : product}
    kwargs_ill = {
            'compute_shadows' : False,
            'sun_zenith_angle' : None,
            'sun_azimuth_angle' : None,
            'height_factor' : height_factor,
            }
            
    bounds = (lonmin, latmin, lonmax, latmax)

    alt, alt_meta = get_elevation(bounds, **kwargs_ele)
    args = alt, alt_meta, height_factor, False, theta_sun, phi_sun,
    ill = compute_brightness(*args)
    ill_meta = alt_meta

    plt.close()
    figsize = (16, 9)
    plt.figure(figsize=figsize)
    N = len(thetas_sun_deg)
    cmap = plt.get_cmap('binary_r')
    for n in range(N):
        sun_azimuth_angle = phis_sun_deg[n]
        sun_zenith_angle = thetas_sun_deg[n]
        kwargs_ill['sun_zenith_angle'] = thetas_sun_deg[n]
        kwargs_ill['sun_azimuth_angle'] = phis_sun_deg[n]
        sun = get_sun_vector(sun_zenith_angle, sun_azimuth_angle)
        # ill, ill_meta = get_brightness(bounds, kwargs_ele, kwargs_ill)

        nplot = n + 1
        plot(ill, ill_meta, nplot, cmap, sun)
    plt.tight_layout()
    plt.show()
