# standard modules
import datetime as dt

# PyPI modules
import numpy as np

# misc
from misc.text import string_utils

################################################################
# chrono arithmetics                                           #
################################################################
def fraction_done(chrono):
    total_count = chrono.total_count
    if total_count > 0:
        return 1. * chrono.count / chrono.total_count
    else:
        return 1.

def fraction_todo(chrono):
    return 1. - fraction_done(chrono)

def time_done(chrono):
    """Return a dt.timedelta."""
    return chrono.global_timer.get('t')

def time_todo(chrono):
    """Return a dt.timedelta."""
    if chrono.count == chrono.total_count:
        return dt.timedelta()

    time = time_total(chrono)
    if time is None:
        return None
    return time - time_done(chrono)

def time_total(chrono):
    """Return a dt.timedelta."""
    seconds_done = time_done(chrono).total_seconds()
    fraction = fraction_done(chrono)
    if fraction == 0:
        return None
    seconds_total = seconds_done / fraction
    try:
        return dt.timedelta(seconds=seconds_total)
    except OverflowError:
        return None

def speed(chrono):
    seconds_done = time_done(chrono).total_seconds()
    return chrono.count / seconds_done


################################################################
# srings                                                       #
################################################################
def time_string(x):
    """Convert timedelta to str.

        Parameters
        ----------
        x : dt.timedelta or None

        Returns
        -------
        str

        Author
        ------
        Andreas Anhaeuser (AA) <anhaeus@meteo.uni-koeln.de>
        Institute for Geophysics and Meteorology
        University of Cologne, Germany

        History
        -------
        2017-04-06 (AA): Created
    """
    ###################################################
    # SPECIAL CASES                                   #
    ###################################################
    # None
    if x is None:
        return 'unknown'

    # negative
    if x < dt.timedelta():
        return '-' + time_string(-x)

    ###################################################
    # REGULAR CASE                                    #
    ###################################################
    secs = int(x.total_seconds())
    return str(dt.timedelta(seconds=secs))

def short_time_string(seconds, sep=''):
    """Convert to minutes, hours, days if neccessary.

        Parameters
        ----------
        seconds : float
        sep : str, optional
            separator between number and unit. Default: ''

        Returns
        -------
        str
            something like '3.0s', '51min', etc

        Author
        ------
        Andreas Anhaeuser (AA) <anhaeus@meteo.uni-koeln.de>
        Institute for Geophysics and Meteorology
        University of Cologne, Germany

        History
        -------
        2017-03-08 (AA): Created
    """
    x = seconds
    units = 's'
    # minutes
    if abs(x) > 60:
        x /= 60
        units = 'min'
    # hours
    if abs(x) > 60:
        x /= 60
        units = 'h'
        # days
        if abs(x) > 24:
            x /= 24
            units = 'd'
    return '%s%s' % (string_utils.human_format(x, sep=sep), units)

def nice_time_string(seconds):
    """Return N.Mu or NNu or NNNu."""
    value = seconds
    units = 's'
    if value > 300:
        value /= 60
        units = 'm'
        if value > 300:
            value /= 60
            units = 'h'
            if value > 100:
                value /= 24
                units = 'd'

    if value >= 10:
        fmt = '%1.0f'
    else:
        fmt = '%1.1f'

    return fmt % value + units

def short_time_string_no_frac(seconds, digits=3, sep=''):
    x = int(np.round(seconds))
    units = 's'
    # minutes
    if len(str(x)) > digits:
        x //= 60
        units = 'm'
    # hours
    if len(str(x)) > digits:
        x //= 60
        units = 'h'
    # days
    if len(str(x)) > digits:
        x //= 24
        units = 'd'
    # years
    if len(str(x)) > digits:
        x //= 365
        units = 'a'
    return '%s%s%s' % (str(x), sep, units)

def count_string(count):
    if count < 0:
        return '-' + count_string(-count)

    plain = '%1.0f' % count
    len_plain = len(plain)
    nice = plain[-3:]
    pos = -3
    while pos > - len_plain:
        nice = plain[pos-3:pos] + ',' + nice
        pos -= 3
    return nice

################################################################
# chronometer sring functions                                  #
################################################################
def speed_string(chrono):
    x = speed(chrono)
    if x == 0:
        return '---'
    units = 's'

    # minutes
    if x < 1:
        x *= 60
        units = 'min'

    # hours
    if x < 1:
        x *= 60
        units = 'h'

    # days
    if x < 1:
        x *= 24
        units = 'd'

    return '%s %s/%s' % (
            string_utils.human_format(x), chrono.item_plural, units,
            )

def inverse_speed_string(chrono):
    s = speed(chrono)
    if s == 0:
        return '---'
    inverse_speed = 1 / s
    time_str = short_time_string(inverse_speed, sep=' ')
    return time_str + '/' + chrono.item_name

def end_string(chrono, fmt):
    now = dt.datetime.now()
    time = time_todo(chrono)
    if time is None:
        return '---'
    try:
        end = now + time
        return end.strftime(fmt)
    except OverflowError:
        return '---'

################################################################
# grammar                                                      #
################################################################
def plural(word):
    """Return `word` + '(e)s'."""
    if word[:1] in ('s', 'x', 'z'):
        return word + 'es'
    else:
        return word + 's'
