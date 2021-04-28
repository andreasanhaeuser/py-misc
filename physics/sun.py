#!/usr/bin/python
"""Functions related to Sun's position on the sky.

    Author
    ------
    2014-2018
    Andreas Anhaeuser (AA) <andreas.anhaeuser@posteo.net>
    Institute for Geophysics and Meteorology
    University of Cologne, Germany
"""
# standard modules
import datetime as dt

# PyPI modules
import numpy as np
import ephem

def utc_sunrise_sunset(d, lon, lat, alt=0, pres=None, temp=None):
    """Return SR and SS for this day as a pair of dt.datetime.

        The time of day of the input is ignored, only the date is taken into
        account.

        The horizon is assumed to be at 0 degrees, regardless of the altitude.

        Parameters
        ----------
        d : datetime.date or datetime.datetime
        lon: float
            the observer's longitude in deg North
        lat: float
            the observer's latitude in deg East
        alt: float
            the observer's altitude above the surface in meters
        pres : float, optional
            (Pa) pressure at ground. Used to compute atmospheric light
            refraction. Note: overrides alt !
        temp : float, optional
            (K) temperature at ground. Used to compute atmospheric light
            refraction. Default: 288.15.

        Returns
        -------
        SR, SS : datetime.datetime or None
            SR and SS denote sun rise and sun set of the day specified by d.
            If (in polar regions,) there is no sun set or sun rise at this day,
            SR and/or SS are None.

        Notes
        -----
        All times and dates are treated as UTC (not local time). This means
        that SS can be before SR (this is roughly the case for places that lie
        on the hemisphere where abs(lon) > 90).

        SR and SS as returned by this function are always on the same UTC date
        as d. If such a SR and/or SS does not exist, then the respective value
        is set to None. This may be counter-intuitive in some cases, e. g. on a
        day where the previous sun rise is just before 0:00 and the following
        is just after 0:00 (i. e. 24 h and some minutes later), then on the
        actual day, there is no sun rise and None is returned:

              | day (n-1) | day n   | day (n+1)
        ------+-----------+---------+-----------
        SR    | 23:59     | ----    | 00:01
        SS    | 11:20     | 11:18   | 11:16

        Dependencies
        ------------
        Makes use of the package ephem.

        Raises
        ------
        TypeError, ValueError

        Tested
        ------
        Moderately tested. Seems bug-free, but not 100% sure.
    """
    #############################################
    # INPUT CHECK                               #
    #############################################
    if not isinstance(d, dt.date):
        raise TypeError('d must be an instance of datetime.date.')
    if not -180 <= lon <= 180:
        raise ValueError('lon must be between -180 and 180.')
    if not -90 <= lat <= 90:
        raise ValueError('lat must be between -90 and 90.')

    #############################################
    # SET TIME TO MIDNIGHT                      #
    #############################################
    if isinstance(d, dt.datetime):
        d = d.date()

    #############################################
    # CREATE EPHEM OBJECTS                      #
    #############################################
    # create Sun object:
    sun = ephem.Sun(d)

    # create observer object:
    site = ephem.Observer()
    site.lon = str(lon)
    site.lat = str(lat)
    site.elevation = alt
    if pres is None:
        site.compute_pressure()
    else:
        site.pressure = pres * 1e-2  # (convert from Pa to hPa)
    if temp is not None:
        site.temp = temp - 273.15    # (convert from deg C to K)
    site.date = d

    #############################################
    # SR AND SS                                 #
    #############################################
    try:
        SR = site.next_rising(sun).datetime()
        # make sure SR is on the same day:
        if SR.date() != d:
            SR = None
    except ephem.NeverUpError:
        SR = None
    except ephem.AlwaysUpError:
        SR = None

    try:
        SS = site.next_setting(sun).datetime()
        # make sure SS is on the same day:
        if SS.date() != d:
            SS = None
    except ephem.NeverUpError:
        SS = None
    except ephem.AlwaysUpError:
        SS = None

    return (SR, SS)

def is_day(d, lon, lat, alt=0., pres=None, temp=None):
    """Return a bool.

        Consistent with utc_sunrise_sunset within tens of microseconds.

        Parameters
        ----------
        d : datetime.date or datetime.datetime
            UTC
        lon: float
            the observer's longitude in deg North
        lat: float
            the observer's latitude in deg East
        alt: float
            the observer's altitude above the surface in meters
        pres : float, optional
            (Pa) pressure at ground. Used to compute atmospheric light
            refraction. Note: overrides alt !
        temp : float, optional
            (K) temperature at ground. Used to compute atmospheric light
            refraction. Default: 288.15.

        Returns
        -------
        bool

        Dependencies
        ------------
        Makes use of the package ephem.

        Raises
        ------
        TypeError, ValueError

        Tested
        ------
        Moderately tested. Seems bug-free, but not 100% sure.
    """
    #############################################
    # INPUT CHECK                               #
    #############################################
    if not isinstance(d, dt.datetime):
        raise TypeError('d must be datetime.datetime.')
    if not -180 <= lon <= 180:
        raise ValueError('lon must be between -180 and 180.')
    if not -90 <= lat <= 90:
        raise ValueError('lat must be between -90 and 90.')

    #############################################
    # CREATE EPHEM OBJECTS                      #
    #############################################
    # create Sun object:
    sun = ephem.Sun(d)

    # create observer object:
    site = ephem.Observer()
    site.lon = str(lon)
    site.lat = str(lat)
    site.elevation = alt
    site.date = d
    if pres is None:
        site.compute_pressure()
    else:
        site.pressure = pres * 1e-2  # (convert from Pa to hPa)
    if temp is not None:
        site.temp = temp - 273.15    # (convert from deg C to K)

    # compute sun elevation
    sun.compute(site)
    elevation = sun.alt

    # take into account extent of the sun:
    size_arcsec = sun.size
    size = size_arcsec *np.pi / (3600 *180)
    elev_top = elevation + size / 2
    return elev_top >= 0

def is_night(d, lon, lat):
    """Return a bool.

    Description, see `is_day`."""
    return not is_day(d, lon, lat)

def last_sunrise(*args, **kwargs):
    """Return a dt.datetime.

        Regardless of the altitude, the horizon is assumed to be at 0 deg (for
        an observer at zero altitude).

        Parameters
        ----------
        time : datetime.date or datetime.datetime
        lon: float
            the observer's longitude in deg North
        lat: float
            the observer's latitude in deg East
        alt: float
            the observer's altitude above the surface in meters
        pres : float, optional
            (Pa) pressure at ground. Used to compute atmospheric light
            refraction. Note: overrides alt !
        temp : float, optional
            (K) temperature at ground. Used to compute atmospheric light
            refraction. Default: 288.15.

        Returns
        -------
        datetime.datetime

        Notes
        -----
        Unlike utc_sunrise_sunset, this function always returns a
        datetime.datetime object. If necessary, it goes back several days until
        it finds a sunrise.
    """
    event_name = 'last_sr'
    return find_neighbouring_event(event_name, *args, **kwargs)

def last_sunset(*args, **kwargs):
    """Return a dt.datetime.

        For description, see last_sunrise
    """
    event_name = 'last_ss'
    return find_neighbouring_event(event_name, *args, **kwargs)

def next_sunrise(*args, **kwargs):
    """Return a dt.datetime.

        For description, see last_sunrise
    """
    event_name = 'next_sr'
    return find_neighbouring_event(event_name, *args, **kwargs)

def next_sunset(*args, **kwargs):
    """Return a dt.datetime.

        For description, see last_sunrise
    """
    event_name = 'next_ss'
    return find_neighbouring_event(event_name, *args, **kwargs)

def find_neighbouring_event(event_name, time, lon, lat, alt=0, pres=None, temp=None):
    """Return a datetime.datetime.

        Regardless of the altitude, the horizon is assumed to be at 0 deg.

        Parameters
        ----------
        event_name : str
            {'next_sr', 'next_ss', 'last_sr', 'last_ss'}
        time : datetime.date or datetime.datetime
        lon: float
            the observer's longitude in deg North
        lat: float
            the observer's latitude in deg East
        alt: float
            the observer's altitude above the surface in meters
        pres : float, optional
            (Pa) pressure at ground. Used to compute atmospheric light
            refraction. Note: overrides alt !
        temp : float, optional
            (K) temperature at ground. Used to compute atmospheric light
            refraction. Default: 288.15.

        Returns
        -------
        datetime.datetime

        Notes
        -----
        Unlike utc_sunrise_sunset, this function always returns a
        datetime.datetime object. If necessary, it goes back or forth several days until
        it finds the event.

        Tested
        ------
        Moderately tested. Seems bug-free, but not 100% sure.

        History
        -------
        2018-11-17 (AA): Bug-fix: Close to poles, last event can be up to a
                         whole year away (was: half a year).
        2018-11-17 (AA): Created this generic sub-function to last_sunrise,
                         next_sunset, etc
    """
    # Idea
    # ====
    # 1. If sunrise on this day is earlier than d, return it.
    # 2. If sunrise is later than d, go back one day.
    # 3. In polar regions, it may be necessary to go back several days to find
    #    a sun rise, hence the while-loop.

    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    event_name = event_name.lower()
    allowed_event_names = ('next_sr', 'next_ss', 'last_sr', 'last_ss')
    assert event_name in allowed_event_names

    ###################################################
    # CREATE OBSERVER OBJECT                          #
    ###################################################
    site = ephem.Observer()
    site.lon = str(lon)
    site.lat = str(lat)
    site.elevation = alt
    if pres is None:
        site.compute_pressure()
    else:
        site.pressure = pres * 1e-2
    if temp is not None:
        site.temp = temp - 273.15
    site.date = time

    ###################################################
    # EVENT FUNCTION & TIME INCREMENT                 #
    ###################################################
    # create a generic function `event_function` that points to the appropriate
    # method of `site` for the specific event in question
    if event_name == 'next_sr':
        event_function = site.next_rising
        increment = dt.timedelta(days=1)

    elif event_name == 'next_ss':
        event_function = site.next_setting
        increment = dt.timedelta(days=1)

    elif event_name == 'last_sr':
        event_function = site.previous_rising
        increment = dt.timedelta(days=-1)

    elif event_name == 'last_ss':
        event_function = site.previous_setting
        increment = dt.timedelta(days=-1)

    ###################################################
    # FIND SUN RISE                                   #
    ###################################################
    # in extreme cases (close to the pole), the event in question may be up to a
    # year away; thus the loop.

    pivot = time
    found = False
    while not found:
        # make sure it is not more than a year away
        assert abs(pivot-time) <= dt.timedelta(days=366)

        sun = ephem.Sun(pivot)
        site.date = pivot

        try:
            event = event_function(sun).datetime()
            found = True
        except ephem.NeverUpError:
            pass
        except ephem.AlwaysUpError:
            pass

        pivot = pivot + increment
    return event

def fraction_of_day(time, sunrise, sunset):
    """Return a float between 0 and 1.

        Divides the day into daytime and nighttime. Daytime ranges from 0 to
        0.5, nighttime from 0.5 to 1.0. If day and night do not have the same
        (physical) length, then the mapping from time to day_fraction is NOT
        linear, i. e. the interval from day_fraction 0.4 to 0.5 is longer or
        shorter (in seconds) than the interval from 0.5 to 0.6.

        Parameters
        ----------
        time : datetime.datetime
        sunrise : datetime.datetime
        sunset : datetime.datetime

        Returns
        -------
        float
            0. : sun rise
            0.5 : sun set
            1. sun rise
    """
    for x in [time, sunrise, sunset]:
        if not isinstance(x, dt.datetime):
            raise TypeError('time, sunrise and sunset must be ' +
                    'datetime.datetime objects.')
    if not time.date() == sunrise.date() == sunset.date():
        raise ValueError('time, sunrise and sunset must be on the same day.')

    # abbreviations:
    t = time
    SR = sunrise
    SS = sunset

    #####################################
    # CASE 1: sunrise before sunset     #
    #####################################
    if SR <= SS:
        ld = (SS - SR).total_seconds()
        ln = 86400. - ld
        if t < SR:
            diff = ln - (SR - t).total_seconds()
            df = 0.5 + 0.5 * diff / ln
        elif SR <= t < SS:
            diff = (t - SR).total_seconds()
            df = 0.5 * diff / ld
        else:
            diff = (t - SS).total_seconds()
            df = 0.5 + 0.5 * diff / ln

    #####################################
    # CASE 2: sunrise after sunset      #
    #####################################
    else:
        ln = (SR - SS).total_seconds()
        ld = 86400. - ln
        if t < SS:
            diff = ld - (SS - t).total_seconds()
            df = 0.5 * diff / ln
        elif SS <= t < SR:
            diff = (t - SS).total_seconds()
            df = 0.5 + 0.5 * diff / ld
        else:
            diff = (t - SR).total_seconds()
            df = 0.5 * diff / ln
    return df

def sun_position(d, lon, lat):
    """Calculate elevation (= altitude) and azimuth.
    
        For those unfamiliar with azimuth and altitude: They describe position
        in the sky by measuring angle around the horizon, then angle above the
        horizon.
        
        Parameters
        ----------
        d : datetime.date or datetime.datetime
        lon : float
            (deg) the observer's longitude
        lat : float
            (deg) the observer's latitude
        
        Returns
        -------
        altitude : float
            (deg)
        azimuth : float
            (deg)

        Authors
        -------
        CF : Christopher Frank <cfrank@meteo.uni-koeln.de>
        AA : Andreas Anhaeuser <andreas.anhaeuser@posteo.net>
        Institute for Geophysics and Meteorology
        University of Cologne, Germany
        http://www.geomet.uni-koeln.de/en/home/


        History
        -------
        2019-05-19 (AA): Debug (wrong variable names, missing numpy-import)
        2016-09-09 (AA): Slightly modified
        2016       (CF): Written
    """
    ###################################################
    # CREATE SUN OBJECT                               #
    ###################################################
    sun = ephem.Sun(d)

    ###################################################
    # CREATE OBSERVER OBJECT                          #
    ###################################################
    site = ephem.Observer()
    site.lon = str(lon)
    site.lat = str(lat)
    site.date = d

    sun.compute(site)
    altitude_rad = sun.alt  # radians
    azimuth_rad = sun.az    # radians
    altitude_deg = np.degrees(altitude_rad)
    azimuth_deg = np.degrees(azimuth_rad)
    return altitude_deg, azimuth_deg


###################################################
# TESTING                                         #
###################################################
if __name__ == '__main__':
    lon = 0
    lat = 85.
    import datetime_utils as dt_utils
    import matplotlib.pyplot as plt
    import sys

    time_beg = dt.datetime(2018, 1, 1)
    time_end = dt.datetime(2019, 1, 1)
    time_inc = dt.timedelta(1)
    times = dt_utils.datetime_range(time_beg, time_end, time_inc)
    N = len(times)
    sunrises = [None] * N
    for n, time in enumerate(times):
        sys.stderr.write("\r%s" % str(time))
        sunrise = last_sunrise(time, lon, lat)
        sunrises[n] = sunrise

    diffs = [(sunrises[n] - times[n]).days for n in range(N)]
    plt.plot(times, sunrises)
    plt.show()
