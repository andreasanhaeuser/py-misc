#!/usr/bin/python3
"""Compute object height using its shadow.

    The height of an object is computed using the coordinates of its base and
    the tip of its shadow together with the calendar date. This best works with
    thin objects (towers, masts, chimneys).

    The date matters a lot and must thus be known. "Guessing" the date won't
    work. Make sure you use the correct date format (use -f option if
    necessary).


    Output
    =====
    Height of the building in metres.


    Dependencies
    ============
    * python3
    * python3 module `ephem`


    Assumptions
    ===========
    * The building is vertical.
    * The the shadow tip is at the same elevation as the stack base (flat
      terrain). If, for some reason, you know the elevation difference between
      stack base and shadow tip, you can correct for it by simple
      subtraction/addition to the result.


    Precision
    =========
    The only approximation is the assumtion of a locally flat Earth surface,
    which is very close to exact for all practical purposes
    (building height << Earth's radius).

    Otherwise, the computation is astronomically and mathematically exact.

    Imprecisions thus stem solely from measurement and deviations from the
    above assumptions:
    - elevation difference between stack base and shadow tip
    - building not vertical

    In my tests, I achieved <10% error in flat terrain with little effort.


    Outline of the algorithm
    ========================
    [Caution: `altitude` is an angle here (-> celestial coordinates.)]

    1. From the shadow coordinates, its length and orientation (angle from
       North) are computed.
       (skipped if shadow length and orientation are given directly).

    2. The different positions of the sun in the sky throughout the day in
       question are computed.
       This is done by the `ephem` module. It gives the position as azimuth
       (angle from North) and altitude (angle to the horizon).

    3. The program searches for the time of the day when the azimuth of the sun
       is opposite to the shadow orientation.

    4. The altitude of the sun at this instant is then used to retrieve the
       building height by basic trigonometry.


    Testing
    =======
    Script has been tested against buildings of known height.


    Licence
    =======
    * Feel free the use, modify, develop and distribute the script without
      notice to the author.
    * The author doesn't mind being mentioned if you use the script.
    * The author denies all liability for potential software or hardware damage
      or incorrect data. You use this software on your own risk and
      responsibility.


    Author
    ======
    AA : Andreas Anhaeuser
        <andreas.anhaeuser@greenpeace.org>       (currently regularily checked)
        <andreas.anhaeuser.data_analyst@posteo.net> (stable,for the far future)


    History
    =======
    2021-03-12 (AA): Implemented -p option for direct polar coordinate input.
    2020-12-22 (AA): Created.
"""
# standard modules
import sys
import argparse
import datetime as dt
import math
from math import pi

# PyPI modules
import ephem

# defaults
# =============================
_planet_radius = 6371e3
_date_format = '%Y%m%d'
_time_resolution = 60
# =============================

################################################################
# parser for command-line arguments                            #
################################################################
description = __doc__.split('\n')[0]

FmtClass = argparse.RawDescriptionHelpFormatter
parser = argparse.ArgumentParser(
        description=description, formatter_class=FmtClass,
        )

# Show verbose help
# -----------------
# this is a bit hack-ish, I don't know how to do that in a clean way:
parser.add_argument(
        '-d', '--doc', action='version', version=__doc__,
        help='show full documentation and exit',
        )


# base coordinates
# ----------------
# lon
parser.add_argument(
        dest='lon_base_deg', type=float,
        help='(deg) longitude of the object base',
        )
# lat
parser.add_argument(
        dest='lat_base_deg', type=float,
        help='(deg) latitude of the object base',
        )


# tip coordinates
# ---------------
# either (lon, lat) or (length, azimuth)

# tip_coordinate_a
help_tca = (
        'Without -p option: longitude of the shadow tip (deg).'
        + '\nWith -p option: length of the shadow (m).'
        )
parser.add_argument(dest='tip_coordinate_a', type=float, help=help_tca)

# tip_coordinate_b
help_tcb = (
        'Without -p option: latitude of the shadow tip (deg).'
        + '\nWith -p option: azimuth of the shadow (deg)'
        + ' (angle clockwise from North)'
        )
parser.add_argument(dest='tip_coordinate_b', type=float, help=help_tcb)


# date
# ----
parser.add_argument(
        dest='date_stamp', type=str,
        help='(yyyymmdd) date of the scene',
        )


# optional arguments
# ------------------
# polar coordinates (shadow length + orientation)
polar_help = (''
        + 'Interpret tip coordinates as (length, azimuth_deg).'
        + '\nlength: Shadow length in metres.'
        + '\nazimuth_deg: angle from North (clockwise) in degrees.'
        )
parser.add_argument(
        '-p', '--polar', dest='polar_coordinates', action='store_true',
        help=polar_help,
        )

# date format
parser.add_argument(
        '-f', dest='date_format', type=str, default=_date_format,
        help='custom format for date stamp (eg. %%Y-%%m-%%d for yyyy-mm-dd)',
        )

# planet radius 
parser.add_argument(
        '-r', dest='planet_radius', type=float,
        default=_planet_radius,
        help=('(m) radius of the planet. Default: %i' % _planet_radius),
        )

# time resolution
parser.add_argument(
        '-t', dest='time_resolution', type=int, default=_time_resolution,
        help='(s) temporal resolution. Default: %i' % _time_resolution,
        )

# verbose output
parser.add_argument(
        '-v', '--verbose', dest='verbose', action='store_true',
        help='show verbose output. Otherwise only the height is shown',
        )


args = parser.parse_args()



################################################################
# helper functions                                             #
################################################################
def standard_angle(alpha):
    """Wrap an angle to [-pi..pi)."""
    alpha = ((alpha + pi) % (2*pi)) - pi
    assert -pi <= alpha < pi
    return alpha

def absolute_separation(alpha, beta):
    """Return angle separation [0..pi]."""
    sep = (alpha - beta) % (2*pi)
    if sep > pi:
        sep = - (sep - 2*pi)
    assert 0 <= sep <= pi
    return sep


################################################################
# main                                                         #
################################################################
# extract input variables
# =================================================
time = dt.datetime.strptime(args.date_stamp, args.date_format)
lat_base_deg =  args.lat_base_deg
lon_base_deg = args.lon_base_deg

# Check whether base coordinates are in reasonable range
if not (-90 <= lat_base_deg <= 90):
    message = 'Base latitude must be in [-90, 90], got %s' % str(lat_base_deg)
    raise ValueError(message)
if not (-180 <= lon_base_deg <= 300):
    print(
            'WARNING: Base longitude %s is outside [-180, 360],'
            % str(lon_base_deg)
            + " but I'm assuming you know what you are doing..."
            )


############################################################
# compute/retrieve shadow length and orientation           #
############################################################
if args.polar_coordinates:
    # Tip coordinates are interpreted as (length, orienation):
    # =================================================
    len_shadow = args.tip_coordinate_a      # (m)
    azi_shadow_deg = args.tip_coordinate_b  # (degrees)
    # =================================================

    # deg -> rad
    # =================================================
    azi_shadow = math.radians(azi_shadow_deg)   # (rad)
    # =================================================
else:
    # Tip coordinates are interpreted as (lon, lat):
    # =================================================
    lon_tip_deg = args.tip_coordinate_a     # (degrees)
    lat_tip_deg = args.tip_coordinate_b     # (degrees)
    # =================================================

    # deg -> rad
    # =================================================
    lon_base = math.radians(lon_base_deg)
    lat_base = math.radians(lat_base_deg)
    lon_tip = math.radians(lon_tip_deg)
    lat_tip = math.radians(lat_tip_deg)
    # =================================================

    # compute shadow length and orientaion
    # =================================================
    R = args.planet_radius
    lat_mean = (lat_base + lat_tip) / 2

    # local cartesian coordinates (metres)
    x = R * (lon_tip - lon_base) * math.cos(lat_mean)
    y = R * (lat_tip - lat_base)

    # length & orientation
    len_shadow = math.sqrt(x**2 + y**2)
    azi_shadow = math.atan2(x, y)
    # =================================================


# prepare search loop
# =================================================
# separation between sun rays and shadow orientations
# (the loop below minimizes this value)
closest_separation = math.inf

time_inc = dt.timedelta(seconds=args.time_resolution)

# end of loop
time_max = time + dt.timedelta(1)   # one day later

found_any = False
# =================================================

while time < time_max:
    time += time_inc

    # create the observer (i. e. the shadow casting object)
    site = ephem.Observer()
    site.lon = str(lon_base_deg)
    site.lat = str(lat_base_deg)
    site.date = time

    # create the sun
    sun = ephem.Sun(time)
    sun.compute(site)

    # retrieve the sun's position in the sky
    # altitude: angle above horizon
    # azimuth: angle clockwise from North
    alt_sun = standard_angle(sun.alt)
    azi_sun = standard_angle(sun.az)

    # discard night times
    if alt_sun < 0:
        continue

    # ray of light has azimuth opposite to sun
    azi_ray = azi_sun + pi

    # separation between ray of light and shadow
    separation = absolute_separation(azi_ray, azi_shadow)

    # check whether ray of light and shadow are aligned
    if separation < closest_separation:
        closest_separation = separation
        closest_altitude = alt_sun
        closest_time = time
        found_any = True

if not found_any:
    date_str = args.date_stamp
    coord_str = '%1.5f, %1.5f (lon, lat)' % (lon_base_deg, lat_base_deg)
    print('ERROR: Sun is never up on %s at %s' % (date_str, coord_str))
    sys.exit(1)


# compute object height
height = len_shadow * math.tan(closest_altitude)

if args.verbose:
    if args.time_resolution < 60:
        time_fmt = '%Y-%m-%d %H:%M:%S'
    else:
        time_fmt = '%Y-%m-%d %H:%M'
    time_utc = closest_time

    # compute mean solar time
    offset_hours = ((lon_base_deg/360*24) + 12) % 24 - 12
    offset = dt.timedelta(hours=offset_hours)
    time_mst = time_utc + offset

    time_utc_str = time_utc.strftime(time_fmt)
    time_mst_str = time_mst.strftime(time_fmt)

    print('-------------------------')
    print('Shadow length:  %1.1f m' % len_shadow)
    print('Shadow azimuth: %1.1f deg' % (math.degrees(azi_shadow) % 360))
    print('-------------------------')
    print('Time:           %s (UTC)' % time_utc_str) 
    print('                %s (local mean solar time)' % time_mst_str)
    print('Height:         %1.1f m' % height)
else:
    print('%1.1f' % height)
