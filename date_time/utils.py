"""A collection of functions that act on classes of the datetime module.

    Author
    ------
    Written in 2014-2019
    Andreas Anhaeuser (AA)
    <andreas.anhaeuser@posteo.net>
    Institute for Geophysics and Meteorology
    University of Cologne, Germany
"""

# standard modules
import calendar
import datetime as dt
from collections import Iterable

# PyPI modules
import numpy as np

# internal modules
from . import timedelta
# from .intervals import Interval, DaytimePeriod, Season
from .timedelta import MonthTimeDelta, YearTimeDelta
from .intervals.interval import Interval
from .intervals.daytime_period import DaytimePeriod
from .intervals.season import Season

###################################################
# range                                           #
###################################################
def date_range(beg, end, inc=None, season_of_year=None):
    """Return a list of dt.date. Behaves in analogy to datetime_range().

        Works as datetime_range(). Refer there for documentation.

        Parameters
        ----------
        beg : datetime.date
            inclusive
        end : datetime.date
            exclusive
        inc : datetime.timedelta, optional
            increment. Must be a multiple of 1 day. Default: 1 day
        season_of_year: Season, optional
            Default: whole year

        Returns
        -------
        list of datetime.date

        History
        -------
        2018-12-20 (AA): Force `inc` to be a multiple of 1 day
        2018-01-04 (AA): Made `inc` an optional argument (default: 1 day).
        2017-12-30 (AA): Created
    """
    # ========== default ================================= #
    if inc is None:
        inc = dt.timedelta(days=1)

    # ========== input check ============================= #
    # make sure inc is multiple of a day
    if inc % dt.timedelta(days=1) != dt.timedelta():
        raise ValueError('`inc` must be a multiple of 1 day.')

    # ========== normalize bounds ======================== #
    # make them datetime.datetime objects

    # the start is always midnight
    time_beg = dt.datetime.combine(beg, dt.time())

    # the end is set to midnight only if it is not already a dt.datetime
    if isinstance(end, dt.datetime):
        time_end = end
    else:
        time_end = dt.datetime.combine(end, dt.time())

    # ========== construct the list ====================== #
    dtrange = datetime_range(
            beg=time_beg, end=time_end, inc=inc,
            season_of_year=season_of_year)

    # ========== re-cast ================================= #
    # convert to datetime.date objects
    return [t.date() for t in dtrange]

def datetime_range(
        beg, end, inc=None, season_of_year=None, daytime_period=None,
        ):
    """Return a list of dt.datetime. Works in analogy to range().

        The function works similar to the standard range() function and is
        intended to be an extension of it to datetime objects.

        The returned list contains only elements that lie within
        season_of_year and daytime_period.

        Parameters
        ----------
        beg : dt.datetime
            inclusive
        end : dt.datetime
            exclusive
        inc : dt.timedelta, optional
            Default: one day
        season_of_year : Season, optional
            Default: whole year
        daytime_period : DaytimePeriod, optional
            Default : whole day

        Tested
        ------
        Extensively tested and heaviliy used. Should be bug-free.

        History
        -------
        2020-02-28 (AA): Set default for `inc` (1 day)
        2017-12-30 (AA): Extention for inc < 0. Raise error on inc == 0.
        2014       (AA): Created
    """
    ###################################################
    # DEFAULT                                         #
    ###################################################
    if inc is None:
        inc = dt.timedelta(1)

    if season_of_year is None:
        season_of_year = Season(dt.datetime(1, 1, 1), dt.datetime(1, 1, 1))

    if daytime_period is None:
        daytime_period = DaytimePeriod(dt.time(), dt.time())

    ############################################################
    # MonthTimeDelta                                           #
    ############################################################
    if isinstance(inc, MonthTimeDelta):
        times_raw = month_range(beg, end, inc.number)
        times = []
        for time in times_raw:
            if time not in season_of_year:
                continue
            if time not in daytime_period:
                continue
            times.append(time)
        return times

    ############################################################
    # YearTimeDelta                                            #
    ############################################################
    if isinstance(inc, YearTimeDelta):
        times_raw = year_range(beg, end, inc.number)
        times = []
        for time in times_raw:
            if time not in season_of_year:
                continue
            if time not in daytime_period:
                continue
            times.append(time)
        return times

    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    if not isinstance(beg, dt.datetime):
        raise TypeError('start must be datetime.datetime.')

    if not isinstance(end, dt.datetime):
        raise TypeError('end must be datetime.datetime.')

    if not isinstance(inc, dt.timedelta):
        raise TypeError('inc must be datetime.timedelta.')

    if inc == dt.timedelta():
        raise ValueError('inc must not be 0.')

    if not isinstance(season_of_year, Season):
        raise TypeError('season_of_year must be Season.')

    if not isinstance(daytime_period, DaytimePeriod):
        raise TypeError('daytime_period must be DaytimePeriod.')

    ###################################################
    # ABBREVIATIONS                                   #
    ###################################################
    season = season_of_year
    timeperiod = daytime_period

    ###################################################
    # BUILD LIST                                      #
    ###################################################
    # case: positive increment
    if inc > dt.timedelta():
        proceed = lambda d : d < end
    # case: negative increment
    else:
        proceed = lambda d : d > end

    out = []
    d = beg
    while proceed(d):
        if (d in season) and (d in timeperiod):
            out.append(d)
        d += inc
    return out

def month_range(beg, end, inc=1):
    """Return a list of dt.datetime. Works in analogy to range().

        The function works similar to the standard range() function and is
        intended to be an extension of it to datetime objects.

        Parameters
        ----------
        beg : dt.datetime or int
            inclusive. int is interpreted as year
        end : dt.datetime
            exclusive. int is interpreted as year
        inc : int > 0
            (months) increment

        History
        -------
        2021-02-05 (AA): Implemented beg, end as int
        2019-02-07 (AA): Created
    """
    ############################################################
    # cast                                                     #
    ############################################################
    if isinstance(beg, int):
        beg = dt.date(beg, 1, 1)
    if isinstance(end, int):
        end = dt.date(end, 1, 1)

    ############################################################
    # input check                                              #
    ############################################################
    if not isinstance(beg, dt.date):
        raise TypeError('`beg` must be int or dt.date, got %s' % type(beg))
    if not isinstance(end, dt.date):
        raise TypeError('`end` must be int or dt.date, got %s' % type(end))
    if not isinstance(inc, int):
        raise TypeError('`inc` must be int, got %s' % type(end))

    if not inc > 0:
        raise ValueError('`inc` must be positive, got %i' % inc)

    ############################################################
    # build list                                               #
    ############################################################
    times = []
    time = beg
    while time < end:
        times.append(time)
        time = timedelta.add_months(time, inc)

    return times

def year_range(beg, end, inc=1):
    """Return a list of dt.datetime. Works in analogy to range().

        The function works similar to the standard range() function and is
        intended to be an extension of it to datetime objects.

        Parameters
        ----------
        beg : dt.datetime
            inclusive
        end : dt.datetime
            exclusive
        inc : int
            (years) increment

        History
        -------
        2019-02-07 (AA): Created
    """
    assert isinstance(beg, dt.date)
    assert isinstance(end, dt.date)
    assert isinstance(inc, int)
    inc = int(inc)
    assert inc > 0

    times = []
    time = beg
    while time < end:
        times.append(time)
        time = timedelta.add_years(time, inc)

    return times

def months_in_interval(interval):
    """Return all overlapping months as list of Interval.

        Parameters
        ----------
        interval : Interval

        Returns
        -------
        months : list of Interval
    """
    time_beg = interval.start.replace(
            day=1, hour=0, minute=0, second=0, microsecond = 0,
            )
    time_end = interval.end

    months = []
    for time_lo in month_range(time_beg, time_end):
        time_hi = timedelta.add_months(time_lo)
        month = Interval(time_lo, time_hi)
        months.append(month)

    # special case: time_end is exactly at the beginning of the month and
    # inclusive
    if interval.end_inclusive and (time_end not in months[-1]):
        time_lo = time_end
        time_hi = timedelta.add_months(time_lo)
        month = Interval(time_lo, time_hi)
        months.append(month)

    return months

def dates_in_year(year, inc=None, season=None):
    """Return a list of dt.date.

        Parameters
        ----------
        year : int or dt.date
        inc : int or None, optional
            default: None. Day increment
        season : Season or None, optional

        Returns
        -------
        dates: list of dt.date
            all dates in `year` (and `season`, if given)
    """
    ############################################################
    # input check                                              #
    ############################################################
    if isinstance(year, dt.date):
        year = year.year

    if not isinstance(year, int):
        raise TypeError('year must be int or dt.date, got %s' % type(year))

    if season is not None and not isinstance(season, Season):
        raise TypeError('season must be Season, got %s' % type(season))

    ############################################################
    # build list                                               #
    ############################################################
    beg = dt.date(year, 1, 1)
    end = dt.date(year+1, 1, 1)
    return date_range(beg, end, inc=inc, season_of_year=season)


###################################################
# unixtime (seconds since ...)                    #
###################################################
def datetime_to_seconds(d, reference_date=dt.datetime(1970, 1, 1)):
    """Convert datetime.datetime to seconds since reference date.

        None values are converted to nan's.

        The function also works for d given as nested lists, but only if the
        length of the sub-lists are equal (i. e. if it is matrix-like).

        Parameters
        ----------
        d : datetime.datetime or datetime.date or None, or list of such
        reference_date : datetime.datetime

        Returns
        -------
        seconds : float or array of such
            (s) time since refernce date

        Development status
        ------------------
        Extensively tested and heavily used. Should be bug-free.

        Note
        ----
        This function is intended to be 100% consistent with its reverse
        function `seconds_to_datetime`. If you find exceptions to this, the
        author would be very thankful for a note!

        See also
        --------
        datetime_to_seconds : reverse function
        """
    ref = reference_date

    # convert date to datetime:
    if ref.__class__ == dt.date:
        ref = dt.datetime(ref.year, ref.month, ref.day)

    if d.__class__ == dt.date:
        d = dt.datetime(d.year, d.month, d.day)

    # recursively call function:
    if isinstance(d, Iterable):
        if len(d) == 0:
            return []
        else:
            s = [datetime_to_seconds(dd, reference_date=ref) for dd in d]
            return np.array(s)

    if d is None:
        return np.nan
    else:
        distance = d - ref                 # (this is an instance of timedelta)
        return distance.total_seconds()

def seconds_to_datetime(seconds, reference_date=dt.datetime(1970, 1, 1)):
    """Convert seconds since reference date to dt.datetime.

        nan's are converted to None-values.

        Parameters
        ----------
        seconds : float or array of such
        reference_date : datetime.datetime

        Returns
        -------
        time : datetime.datetime or None, or list of such
            If the input is given in form of array, the output is a list of
            corresponding shape

        Development status
        ------------------
        Extensively tested and heavily used. Should be bug-free.

        Note
        ----
        This function is intended to be 100% consistent with its reverse
        function `datetime_to_seconds`. If you find exceptions to this, the
        author would be very thankful for a note!

        See also
        --------
        datetime_to_seconds : reverse function
    """
    ref = reference_date

    # convert date to datetime:
    if ref.__class__ == dt.date:
        ref = dt.datetime(ref.year, ref.month, ref.day)

    if isinstance(seconds, Iterable):
        return [seconds_to_datetime(s, ref) for s in seconds]

    if np.isnan(seconds):
        return None
    else:
        return ref + dt.timedelta(seconds=float(seconds))


###################################################
# JULIAN DAYS                                     #
###################################################
def julian_days_to_datetime(days):
    """Return a dt.datetime.

        Parameters
        ----------
        days : float or an Iterable of such

        Returns
        -------
        dt.datetime or list of such

        Tested
        ------
        Moderately. Seems to be working but no guarantee.
    """
    if isinstance(days, Iterable):
        return([julian_days_to_datetime(d) for d in days])
    JD0 = 1721425.5  # Julian date of 1st Jan 1, 00:00
    return dt.datetime(1, 1, 1) + dt.timedelta(days=days - JD0)

def datetime_to_julian_days(time):
    """Return a float.

        Parameters
        ----------
        time : datetime.datetime object or list of such

        Returns
        -------
        float or array of such
    """
    if isinstance(time, Iterable):
        return np.array([datetime_to_julian_days(t) for t in time])
    JD0 = 1721425.5  # Julian date of 1st Jan 1, 00:00
    diff = (time - dt.datetime(1, 1, 1))
    D = diff.days
    S = diff.seconds / 86400.
    U = diff.microseconds / (86400 * 1e6)
    return D + S + U + JD0


###################################################
# DAY OF YEAR                                     #
###################################################
def day_of_year(*args, **kwargs):
    """Alias to doy()."""
    return doy(*args, **kwargs)

def doy(date):
    """Return day of year and year as pair of ints.

        Deprecated. Use sod instead!

        Parameters
        ----------
        date : datetime.date or datetime.datetime or list of such

        Returns
        -------
        doy : int or array of such
            day of year. Between 1 and 366 (inclusive, in leap years)
        year : int or array of such

        Note
        ----
        1st January is DOY 1 (not 0).
    """
    ###################################################
    # RECURSIVE FUNCTION CALL FOR LISTS               #
    ###################################################
    if isinstance(date, Iterable):
        S = np.shape(date)
        assert len(S) == 1
        N = S[0]

        days = np.zeros(N, dtype=int)
        years = np.zeros(N, dtype=int)
        for n in range(N):
            days[n], years[n] = doy(date[n])
        return days, years

    ###################################################
    # TYPE CONVERSION                                 #
    ###################################################
    if isinstance(date, dt.datetime):
        date = date.date()

    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    assert isinstance(date, dt.date)

    year = date.year
    ref = dt.date(year, 1, 1)
    day = (date - ref).days + 1
    return day, year

def date_from_doy(doy, year):
    """Return a datetime.date.

        Parameters
        ----------
        doy : int
            day of year. Can be negative or positve. Can be larger than 366.
        year : int

        Returns
        -------
        date : datetime.date or datetime.datetime

        Note
        ----
        1st January is DOY 1 (not 0).
    """
    return dt.date(year, 1, 1) + dt.timedelta(days=doy-1)

def datetime_from_doy(doy, year, hour=0, minute=0, second=0, microsecond=0):
    """Return a datetime.date.

        Parameters
        ----------
        doy : int
            day of year. Can be negative or positve. Can be larger than 366.
        year : int

        Returns
        -------
        date : datetime.date or datetime.datetime

        Note
        ----
        1st January is DOY 1 (not 0).
    """
    base = dt.datetime(year, 1, 1, hour, minute, second, microsecond)
    increment = dt.timedelta(days=doy-1)
    return base + increment

def days_in_month(date):
    """Return the number of days of that month.

        Parameters
        ----------
        date : dt.date

        Returns
        -------
        int
            number of days of that months.
    """
    if not isinstance(date, int):
        raise TypeError('`date` must be dt.date, got %s' % type(date))

    year = date.year
    month = date.month
    if month in [1, 3, 5, 7, 8, 10, 12]:
        return 31
    if month in [4, 6, 9, 11]:
        return 30
    if calendar.isleap(year):
        return 29
    return 28

def days_in_year(date):
    """Return the number of days of that year.

        Parameters
        ----------
        date : dt.date or int
            if int: the year

        Returns
        -------
        int
            number of days of that months.
    """
    if isinstance(date, dt.date):
        return days_in_year(date.year)

    if not isinstance(date, int):
        raise TypeError('`date` must be dt.date or int, got %s' % type(date))

    year = date
    if calendar.isleap(year):
        return 366
    return 365


###################################################
# SECONDS OF DAY                                  #
###################################################
def sod_doy(time):
    """Return second of day, day of year and year as triple of ints.

        Parameters
        ----------
        time : datetime.datetime or list of such

        Returns
        -------
        sod : int or array of such
            second of day. Between 1 and 86400 (inclusive)
        doy : int or array of such
            day of year. Between 1 and 366 (inclusive, in leap years)
        year : int or array of such

        Notes
        -----
        1st January is DOY 1 (not 0).
        00:00:00 is SOD 1 (not 0).
    """
    ###################################################
    # RECURSIVE FUNCTION CALL FOR LISTS               #
    ###################################################
    if isinstance(time, Iterable):
        S = np.shape(time)
        assert len(S) == 1
        N = S[0]

        secs = np.zeros(N, dtype=int)
        days = np.zeros(N, dtype=int)
        years = np.zeros(N, dtype=int)
        for n in range(N):
            secs[n], days[n], years[n] = sod_doy(time[n])
        return secs, days, years

    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    assert isinstance(time, dt.datetime)

    year = time.year
    month = time.month
    day = time.day

    ref_doy = dt.datetime(year, 1, 1)
    ref_sod = dt.datetime(year, month, day)
    day_of_year = (time - ref_doy).days + 1
    sec_of_day = (time - ref_sod).seconds + 1
    return sec_of_day, day_of_year, year


################################################################
# ROUNDING                                                     #
################################################################
def floor(time, comp='day', step=1):
    """Round to the chosen component.

        Parameters
        ----------
        time : datetime.datetime
        comp : str, optional
            (default: 'day'). Allowed values are {'year', 'month', 'day',
            'hour', 'minute', 'second', microsecond'}
        step : int, optional
            (default: 1). If larger than one, `time` is rounded the highest
            multiple of `step` that is <= `time`.

        Returns
        -------
        datetime.datetime
            rounded to the highest multiple of the chosen component that is not
            larger than the input

        History
        -------
        2019-02-13: (AA) created.
    """
    # constants --------------------------------------
    time_components = ('microsecond', 'second', 'minute', 'hour')
    date_components = ('day', 'month', 'year')
    components = time_components + date_components

    # perform input checks
    # ------------------------------------------------
    # time
    if not isinstance(time, dt.datetime):
        raise TypeError('time must be datetime.datetime.')

    # comp
    if not isinstance(comp, str):
        raise TypeError('comp must be str.')
    if comp not in components:
        raise ValueError('Unknown datetime component: %s' % comp)

    # step
    if not isinstance(step, int):
        raise TypeError('step must be int.')
    step = int(step)
    if step < 1:
        raise ValueError('Step must be a positive int.')
    # ------------------------------------------------

    # set rounded component values
    # ------------------------------------------------
    # - Start from least significant component.
    # - Set them to 0 (or 1 for date components) until the requested component
    #   is reached.
    # - Set the requested component to the adequant multiple of step.
    # - Exit the loop (leaving all more significant components as they are).
    kwargs = {}
    for component in components:
        if component == comp:
            # requested component has been reached
            value = getattr(time, component)
            new_value = value - value % step
            kwargs[component] = new_value

            # all more significant componants remain untouched.
            # --> exit the loop
            break

        # requested component has not yet been reached
        elif component in time_components:
            kwargs[component] = 0
        elif component in date_components:
            kwargs[component] = 1

    # make sure that day and month are >= 1
    for component in date_components:
        if component not in kwargs:
            continue
        if kwargs[component] < 1:
            kwargs[component] = 1

    return time.replace(**kwargs)

###################################################
# STRINGS                                         #
###################################################
def name_of_month(month, kind='full'):
    """Return a str which represents the name of the month.

        Parameters
        ----------
        month : int, between 1 and 12.
            number of the month
        kind : {'full', 'short', 'abbr'}, optional
            'full' : return full month name
            'short' : return the first 3-characters of the month name with no
                      trailing '.'
            'abbr' : return a max. 4-char str, Either the month name, if it is
                     no longer than 4, otherwise its first 3 characters plus
                     trailing '.'
            Default : 'full'

        Returns
        -------
        str
            The month name.
    """
    assert isinstance(month, int)
    assert 0 < month < 13

    names = ['January', 'February', 'March', 'April', 'May', 'June',
            'July', 'August', 'September', 'October', 'November',
            'December']

    if kind == 'full':
        pass
    elif kind == 'short':
        # shorten names
        names = [name[:3] for name in names]
    elif kind == 'abbr':
        # abbreviate names
        for n, name in enumerate(names):
            if len(name) > 4:
                names[n] = name[:3] + '.'
    else:
        raise ValueError('Unrecognized kind: ' + str(kind))

    return names[month-1]

def name_of_weekday(weekday, Nchar=None):
    """Return a str which represents the name of the month.

        Parameters
        ----------
        month : int, between 1 and 12.
            number of the month
        kind : {'full', 'short', 'abbr'}, optional
            'full' : return full month name
            'short' : return the first 3-characters of the month name with no
                      trailing '.'
            'abbr' : return a max. 4-char str, Either the month name, if it is
                     no longer than 4, otherwise its first 3 characters plus
                     trailing '.'
            Default : 'full'

        Returns
        -------
        str
            The month name.
    """
    assert isinstance(weekday, int)
    assert 0 <= weekday < 7
    if Nchar is not None:
        assert isinstance(Nchar, int)
        assert Nchar > 0

    names = [
            'Monday',
            'Tuesday',
            'Wednesday',
            'Thursday',
            'Friday',
            'Saturday',
            'Sunday',
            ]
    name = names[weekday]

    # shorten name
    if Nchar is not None:
        name = name[:Nchar]

    return name

def str_to_datetime(s, fmt=None):
    """Return a dt.datetime.

        If input is insufficient, the function tries to guess smartly to
        extend.

        It understands 'now' and 'today'.

        Parameters
        ----------
        s : str or iterable of such
        fmt : str or None
            The formatter

        Returns
        -------
        datetime.datetime or list of such
    """
    if s is None:
        return dt.datetime.now()

    # check type
    if isinstance(s, str):
        pass
    elif isinstance(s, Iterable):
        return [str_to_datetime(ss, fmt) for ss in s]
    else:
        raise TypeError('Cannot convert type %s' % type(s))

    # now
    if s.lower() == 'now':
        return dt.datetime.now()

    # today
    if s.lower() == 'today':
        date = dt.date.today()
        time = dt.time()
        return dt.datetime.combine(date, time)

    # try to guess formatter
    if fmt is not None:
        pass
    elif len(s) == 8:
        # yyyymmdd
        fmt = '%Y%m%d'
    elif len(s) == 10 and s[4] in ' -_':
        # yyyy?mm?dd
        gap = s[4]
        fmt = '%Y' + gap + '%m' + gap + '%d'
    elif len(s) == 11 and s[8] == '_':
        # yyyymmdd_HH
        fmt = '%Y%m%d_%H'
    elif len(s) == 14:
        # yyyymmddHHMMSS
        fmt = '%Y%m%d%H%M%S'
    elif len(s) == 15 and s[8] == '_':
        # yyyymmdd_HHMMSS
        fmt = '%Y%m%d_%H%M%S'
    elif len(s) == 19:
        # yyyy?mm?dd?HH?MM?SS
        fmt = (''
                + '%Y' + s[4]
                + '%m' + s[7]
                + '%d' + s[10]
                + '%H' + s[13]
                + '%M' + s[16]
                + '%S'
                )
    else:
        message = (
                "I don't know how to convert this to a datetime: '%s'" % s
                + "\nYou can teach me how to by editing the file that"
                + " raises this error."
                )
        raise ValueError(message)

    return dt.datetime.strptime(s, fmt)

def str_to_date(s, fmt=None):
    """Return a dt.datetime.

        If input is insufficient, the function tries to guess smartly to
        extend.

        It understands 'now' and 'today'.

        Parameters
        ----------
        s : str or iterable of such
        fmt : str or None
            The formatter

        Returns
        -------
        datetime.date or list of such
    """
    # check input
    if not isinstance(s, Iterable):
        message = 'Input must be str or iterable of such, got %s' % type(s)
        raise TypeError(message)

    # Iterable
    if not isinstance(s, str):
        return [str_to_date(ss, fmt) for ss in s]

    datetime = str_to_datetime(s, fmt)
    return datetime.date()

def str_to_timedelta(words, enhanced=False):
    if isinstance(words, str):
        separators = ('_', '-')
        for separator in separators:
            words = words.replace(separator, ' ')
        words = words.strip().split()

    if len(words) == 1:
        number = words[0]
        unit = 'd'
        string = number + ' ' + unit
        return str_to_timedelta(string)

    if len(words) != 2:
        message = "Don't know how to interpret this as timedelta: %s" % string
        raise ValueError(message)

    try:
        number = int(words[0])
    except ValueError:
        raise ValueError('Need int in first place, got %s' % number)

    unit = words[1]

    aliases = {
            'us' : 'microseconds',
            'microsecond' : 'microseconds',
            's' : 'seconds',
            'sec' : 'seconds',
            'secs' : 'seconds',
            'min' : 'minutes',
            'mins' : 'minutes',
            'minute' : 'minutes',
            'h' : 'hours',
            'hour' : 'hours',
            'd' : 'days',
            'day' : 'days',
            'w' : 'weeks',
            'week' : 'weeks',
            'mon' : 'months',
            'month' : 'months',
            'a' : 'years',
            'year' : 'years',
            'yr' : 'years',
            'yrs' : 'years',
            }

    if unit in aliases:
        unit = aliases[unit]

    if unit == 'weeks':
        unit = 'days'
        number *= 7

    if enhanced and unit == 'months':
        return MonthTimeDelta(number)

    if enhanced and unit == 'years':
        return YearTimeDelta(number)

    return dt.timedelta(**{unit : number})

def timedelta_to_str(timedelta):
    """Convert a timedelta to number, unit and return as str, str.

        This function is aimes to be fully inverse to `str_to_timedelta`

        Parameters
        ----------
        timedelta : dt.timedelta or MonthTimeDelta or YearTimeDelta

        Returns
        -------
        number : str
        unit : str
    """
    if isinstance(timedelta, MonthTimeDelta):
        num = str(timedelta.number)
        unit = 'month'
        return num, unit

    if isinstance(timedelta, YearTimeDelta):
        num = str(timedelta.number)
        unit = 'a'
        return num, unit

    if not isinstance(timedelta, dt.timedelta):
        raise NotImplementedError(type(timedelta))

    seconds = timedelta.total_seconds()
    if seconds % 1 != 0:
        raise NotImplementedError('Fractions of seconds not implemented.')

    if seconds == 0:
        return 0, 's'

    units_and_factors = (
            ('s', 1),
            ('min', 60),
            ('h', 3600),
            ('d', 86400),
            )

    int_unit, int_factor = units_and_factors[0]

    for unit, factor in units_and_factors:
        if seconds % factor != 0:
            break
        int_unit = unit
        int_factor = factor

    int_number = int(seconds / int_factor)

    return str(int_number), int_unit


################################################################
# MISC                                                         #
################################################################
def decompose(times):
    """Return integer components as dict.

        Parameters
        ----------
        times : None or dt.date or dt.time or iterable of such

        Returns
        -------
        components : dict of int or None
            keys: 'year', 'month', 'day', 'hour', 'minute', 'second',
                    'microsecond'

        """
    ############################################################
    # just one time                                            #
    ############################################################
    if not isinstance(times, Iterable):
        times = [times]
        components = decompose(times)
        for key in components:
            components[key] = components[key][0]
        return components

    ############################################################
    # iterable of times                                        #
    ############################################################
    Nt = len(times)
    keys = ('year', 'month', 'day', 'hour', 'minute', 'second', 'microsecond')

    components = {}
    for key in keys:
        components[key] = [None] * Nt

    for nt, time in enumerate(times):
        if time is None:
            continue

        if isinstance(time, dt.time):
            components['hour'][nt] = time.hour
            components['minute'][nt] = time.minute
            components['second'][nt] = time.second
            components['microsecond'][nt] = time.microsecond
            continue

        if isinstance(time, dt.date):
            components['year'][nt] = time.year
            components['month'][nt] = time.month
            components['day'][nt] = time.day
            continue

        if isinstance(time, dt.datetime):
            components['year'][nt] = time.year
            components['month'][nt] = time.month
            components['day'][nt] = time.day
            components['microsecond'][nt] = time.microsecond
            components['hour'][nt] = time.hour
            components['minute'][nt] = time.minute
            components['second'][nt] = time.second
            components['microsecond'][nt] = time.microsecond
            continue

        raise TypeError('Cannot handle type %s' % type(time))

    return components
