"""Utility functions for BaseInterval and sub-classes."""

# standard modules
import datetime as dt

def daytimediff(minuend, subtrahend, mode=None):
    """Return the difference in daytime regardless of the absolute date.

        Parameters
        ----------
        minuend, subtrahend : datetime.datime or datetime.time

        Returns
        -------
        datetime.timedelta

        mode: (None, 'abs', 'pos', 'neg')
        * None: a value between -12h and + 12h is returned
        * 'abs': a value between 0 and 12h is returned. The absolute difference
        * 'pos': a value between 0 and 24h is returned.
        * 'neg': a value between -24h and 0 is returned.

        Example
        -------
        >>> minuend    = datetime.datetime(2015, 1, 1, 12, 0, 0)
        >>> subtrahend = datetime.datetime(1970,12,16,  8, 0, 0)
        >>> diff = daytimediff(minuend, subtrahend)
        >>> print(repr(diff))
        datetime.timedelta(0, 14400)
        >>> print diff
        4:00:00

        Reliability
        -----------
        Moderately tested.
    """
    ###################################
    # RECURSIVELY CALL FUNCTION       #
    ###################################
    if isinstance(minuend, list):
        return [daytimediff(m, subtrahend) for m in minuend]
    if isinstance(subtrahend, list):
        return [daytimediff(minuend, s) for s in subtrahend]

    ###################################
    # INPUT CHECK                     #
    ###################################
    for d in (minuend, subtrahend):
        if not isinstance(d, (dt.time, dt.datetime)):
            raise TypeError(
                    'Arguments must be datetime.datime or datetime.time.'
                    )

    if isinstance(minuend, dt.datetime):
        m = minuend.replace(1, 1, 1)
    else:
        hour = minuend.hour
        minute = minuend.minute
        second = minuend.second
        microsec = minuend.microsecond
        m = dt.datetime(1, 1, 1, hour, minute, second, microsec)

    if isinstance(subtrahend, dt.datetime):
        s = subtrahend.replace(1, 1, 1)
    else:
        hour = subtrahend.hour
        minute = subtrahend.minute
        second = subtrahend.second
        microsec = subtrahend.microsecond
        s = dt.datetime(1, 1, 1, hour, minute, second, microsec)

    diff = m - s
    if mode is None:
        while diff.total_seconds() > 86400/2:
            diff -= dt.timedelta(days=1)
        while diff.total_seconds() < -86400/2:
            diff += dt.timedelta(days=1)

    if mode == 'abs':
        diff = abs(diff, -diff)

    if mode == 'pos':
        while diff.total_seconds() < 0:
            diff += dt.timedelta(days=1)

    if mode == 'neg':
        while diff.total_seconds() > 0:
            diff -= dt.timedelta(days=1)

    return diff
