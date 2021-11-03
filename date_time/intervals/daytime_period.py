# standard modules
import datetime as dt
from collections import Iterable

# PyPI
import numpy as np

# local
from .base_interval import CyclicInterval
from .interval import Interval


class DaytimePeriod(CyclicInterval):
    """A portion of the day cycle.

        This type is unaware of its absolute (calendar) date.
        It can extend beyond midnight, e. g. 23:00 to 01:00.

        Examples
        --------
        >>> import datetime as dt
        >>> lon = 6.9
        >>> lat = 50.9
        >>> now = dt.datetime.now()
        >>> SR, SS = utc_sunrise_sunset(now, lon, lat)
        >>> day_period = DaytimePeriod(SR, SS)
        >>> if day_period.contains(now):
        >>>     print 'It is day.'
        >>> else:
        >>>     print 'It is night'

        # absolute dates do not matter:
        >>> SR_new = SR + dt.timedelta(days=100)
        >>> SS_new = SS + dt.timedelta(day=-2)
        >>> day_period = DaytimePeriod(SR_new, SS_new)
        >>> if day_period.contains(now):
        >>>     print 'It is day.'
        >>> else:
        >>>     print 'It is night'
    """
    def __init__(self, dt_start=None, dt_end=None, allow_whole_day=True):
        """dt_start and dt_end must be dt.datetime or dt.time objects.

            The absolute calendar dates of the input arguments are ignored.
            This means that calling this function with dt_start 1st Jan 1970
            23:00 and dt_end 30th Mar 1970 01:00 is equivalent to calling it
            with dt_start 1st Jan 2014 23:00 and dt_end 30th Mar 1789 01:00.
            Both will result in a time period between 23:00 and 01:00.

            Boundaries: The lower end dt_start is inclusive, the upper end
            dt_end is exclusive.

            Parameters
            ----------
            dt_start : dt.datetime or dt.time, optional
                default: midnight
            dt_end : dt.datetime or dt.time, optional
                default: midnight
            allow_whole_day : bool, optional
               applies only if start and end are equal (or equivalent) see Note
               below. Default: True.

            Note
            ----
            In case, dt_start == dt_end
            * if allow_whole_day is True, the time period will contain the
              whole day.
            * if allow_whole_day is False, the time period will not contain
              anything.

            In all other cases allow_whole_day is ignored.
        """
        # default ------------------------------------
        if dt_start is None:
            dt_start = dt.datetime(1, 1, 1)
        if dt_end is None:
            dt_end = dt.datetime(1, 1, 1)

        # input check --------------------------------
        for d in (dt_start, dt_end):
            if not isinstance(d, (dt.time, dt.datetime)):
                raise TypeError(
                        'dt_start and dt_end must be instances of '
                        + 'datetime. datetime or datetime.time.'
                        )
        if not isinstance(allow_whole_day, bool):
            raise TypeError('allow_whole_day must be a boolean.')

        # convert to dt.datetime
        # --------------------------------------------
        if isinstance(dt_start, dt.time):
            start = dt.datetime.combine(dt.date(1, 1, 1), dt_start)
        else:
            start = dt_start
        if isinstance(dt_end, dt.time):
            end = dt.datetime.combine(dt.date(1, 1, 1), dt_end)
        else:
            end = dt_end

        # shift to year 1
        # --------------------------------------------
        # self.start will be on Jan 1st 0001 CE.
        # self.end will be between 0 and 1 day later than self.end, all will
        # thus be on Jan 1st or Jan 2nd in year 0001 CE.
        start = start.replace(1, 1, 1)
        end = end.replace(1, 1, 1)

        # check sequence
        # --------------------------------------------
        # make sure that end is not earlier that start:
        if end < start:
            end = end.replace(day=2)
        if end == start and allow_whole_day:
            end = end.replace(day=2)

        # create fields ------------------------------
        self.start = start
        self.end = end

    @classmethod
    def from_interval(cls, interval, allow_whole_day=None):
        """Construct a DaytimePeriod from an Interval.

            interval : Interval
            allow_whole_day : bool or None
                if None, True is assumed for intervals of finite length

        """
        if not isinstance(interval, Interval):
            raise TypeError('Expected Interval, got %s' % interval.__class__)

        if allow_whole_day is None:
            allow_whole_day = interval.length() != 0

        bounds = interval.get_bounds()
        return cls(*bounds, allow_whole_day=allow_whole_day)

    def __repr__(self):
        """Return a str."""
        return 'DaytimePeriod ' + str(self)

    def __str__(self):
        """Return a str."""
        s = self.start.time()
        e = self.end.time()
        if s != e:
            return '[%s, %s)' % (s, e)

        if self.start == self.end:
            return '(empty)'

        return '(full day)'

    def __contains__(self, other):
        """Return a bool or a list of such."""
        # other is iterable
        # -----------------------------------------------
        if isinstance(other, Iterable):
            return self.contains_iterable(other)
        # -----------------------------------------------

        # select comparator
        # -----------------------------------------------
        if isinstance(other, dt.time):
            return self._contains_time(other)
        if isinstance(other, dt.datetime):
            return self._contains_time(other)

        if isinstance(other, Interval):
            return self._contains_interval(other)
        if isinstance(other, DaytimePeriod):
            return self._contains_daytime_period(other)
        # -----------------------------------------------

        # try inverse
        # -----------------------------------------------
        method = 'is_contained_in'
        if hasattr(other, method) and callable(getattr(other, method)):
            return other.is_contained_in(self)
        # -----------------------------------------------

        raise TypeError('Cannot compare to %s' % type(other))

    def is_full_day(self):
        """Return a bool."""
        return self.length() == dt.timedelta(days=1)

    ############################################################
    # helpers                                                  #
    ############################################################
    def _contains_time(self, d):
        """Check whether d is within DaytimePeriod or not.

            Parameters
            ----------
            d : dt.datetime or dt.time or Iterable of such

            Returns
            -------
            bool or list of such
        """

        ####################################
        # CONVERT TO dt.datetime           #
        ####################################
        if d.__class__ is dt.time:
            dd = dt.datetime.combine(dt.date(1, 1, 1), d)
        else:
            dd = d.replace(1, 1, 1)

        ####################################
        # CHECK RELATIVE POSITIONS         #
        ####################################
        # make sure that dd is later than self.start:
        if dd < self.start:
            dd = dd.replace(day=2)
        # check whether dd is earlier than self.end:
        if dd < self.end:
            return True

        return False

    def _contains_interval(self, interval):
        """Return a bool."""
        if interval.length() > self.length():
            return False
        if interval.end_inclusive and interval.end == self.end:
            return False

        period = DaytimePeriod.from_interval(interval)
        return self.__contains__(period)

    def _contains_daytime_period(self, other):
        """Return a bool."""
        if self.start > other.start:
            return False
        if self.end < other.end:
            return False
        return True

    ############################################################
    # helpers: is_contained_in                                 #
    ############################################################
    def is_contained_in(self, other):
        """Return a bool."""
        if isinstance(other, Interval):
            return self.is_contained_in_interval(other)
        if isinstance(other, DaytimePeriod):
            return __other__.contains(self)
        raise TypeError(type(other))

    def is_contained_in_interval(self, interval):
        """Return a bool."""
        if interval.length() >= dt.timedelta(days=1):
            return True
        period = DaytimePeriod.from_interval(interval)
        return period.__contains__(self)

