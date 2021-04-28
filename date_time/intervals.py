#!/usr/bin/python
"""Implentations of time intervals, daytime periods and a seosons."""

import datetime as dt
from collections import Iterable

class Interval(object):
    """A time interval.

        Each bound may be inclusive or exclusive. Default is inclusive for
        lower bound and exclusive for upper bound. Intervals in reverse
        direction are not allowed.

        History
        -------
        2021-04-17 (AA):  Allow int for `start` and `end`
        2021-04-17 (AA):  Add overlap_length()
    """
    def __init__(self, start, end, start_inclusive=True, end_inclusive=False):
        """A time interval.

            Parameters
            ----------
            start : datetime.datetime or int or Iterable of int
                int is interpreted as year
                Iterable of int is interpreted as args to dt.datetime
            end : datetime.datetime or int or Iterable of int
                must not be smaller than `start`
            start_inclusive : bool, optional
                Default: True
            end_inclusive : bool, optional
                Default: False

            `start` must not be larger than `end`.
            if `start` == `end`, the interval cannot be left-open and
            right-closed.
        """
        # convert Iterable -> datetime
        if isinstance(start, Iterable):
            start = dt.datetime(*start)
        if isinstance(end, Iterable):
            end = dt.datetime(*end)

        # convert int -> datetime
        if isinstance(start, int):
            start = dt.datetime(start, 1, 1,)
        if isinstance(end, int):
            end = dt.datetime(end, 1, 1,)

        # convert date -> datetime
        if type(start) == dt.date:
            start = dt.datetime.combine(start, dt.time())
        if type(end) == dt.date:
            end = dt.datetime.combine(end, dt.time())

        # input check
        if not isinstance(start, dt.datetime):
            raise TypeError('`start` must be datetime.datetime.')
        if not isinstance(end, dt.datetime):
            raise TypeError('`end` must be datetime.datetime.')

        # Do not allow intervals in reverse direction
        if end < start:
            raise ValueError('Reverse intervals not implemented')
        elif end == start and end_inclusive and not start_inclusive:
            raise ValueError('Reverse intervals not implemented')

        self.start = start
        self.end = end

        self.start_inclusive = start_inclusive
        self.end_inclusive = end_inclusive

    def __eq__(self, other):
        """Return a bool."""
        if not isinstance(other, Interval):
            return False
        if other.start != self.start:
            return False
        if other.end != self.end:
            return False
        if other.start_inclusive != self.start_inclusive:
            return False
        if other.end_inclusive != self.end_inclusive:
            return False
        return True

    def __lt__(self, other):
        """Return True if strictly before, False otherwise.

            Parameters
            ----------
            other : Interval or datetime.datetime

            Returns
            -------
            bool
                True if self is strictly before other, False otherwise.
        """
        if isinstance(other, Interval):
            return self._completely_earlier(other)
        elif isinstance(other, dt.datetime):
            return self._ends_before(other)
        else:
            raise TypeError('Cannot compare to %s.' % other.__class__)

    def __gt__(self, other):
        """Return True if strictly after, False otherwise.

            Parameters
            ----------
            other : Interval or datetime.datetime

            Returns
            -------
            bool
                True if self is strictly after other, False otherwise.
        """
        if isinstance(other, Interval):
            return other < self
        elif isinstance(other, dt.datetime):
            return self._starts_after(other)
        else:
            raise TypeError('Cannot compare to %s.' % other.__class__)

    def __str__(self):
        """Return a str."""
        if self.start_inclusive:
            lower_par = '['
        else:
            lower_par = '('
        if self.end_inclusive:
            upper_par = ']'
        else:
            upper_par = ')'

        lower_bound = str(self.start)
        upper_bound = str(self.end)

        s = '%s%s, %s%s' % (lower_par, lower_bound, upper_bound, upper_par)
        return s

    def __repr__(self):
        """Return a str."""
        return 'Interval %s' % str(self)

    def __contains__(self, other):
        """Return a bool or a list of such."""
        if isinstance(other, Interval):
            return self._contains_interval(other)
        if isinstance(other, dt.datetime):
            return self._contains_datetime(other)
        if isinstance(other, DaytimePeriod):
            return self._contains_daytime_period(other)
        if isinstance(other, Season):
            return self._contains_season(other)
        raise TypeError('Cannot compare to %s' % type(other))

    def center(self):
        """Return the center of the interval as datetime.time."""
        return self.start + 0.5 * self.length()

    def get_bounds(self):
        """Return a pair of datetime.datetime."""
        return (self.start, self.end)

    ############################################################
    # timedelta returning methods                              #
    ############################################################
    def length(self):
        """Return a datetime.timedelta."""
        return self.end - self.start

    def overlap_length(self, other):
        """Return a datetime.timedelta."""
        if not self.overlaps(other):
            return dt.timedelta()

        overlap = self.intersect(other)
        return overlap.length()

    ############################################################
    # boolean methods                                          #
    ############################################################
    def contains(self, other):
        """Alias for backward compatibility."""
        if isinstance(other, Iterable):
            return [self.contains(item) for item in other]
        return self.__contains__(other)

    def overlaps(self, other):
        """Return True if intervals overlap, False otherwise."""
        if not isinstance(other, Interval):
            raise TypeError('other must be Interval.')
        if self < other:
            return False
        if self > other:
            return False
        return True

    def is_continued_by(self, other):
        """Return True if other is a seamless continuation of self.

            Returns True if and only if both these are True:
            - `other` starts exactly where `self` ends
            - the upper bound of `self` or the lower bound of `other` or both
              are inclusive

            Parameters
            ----------
            other : Interval

            Returns
            -------
            bool
        """
        if self.end != other.start:
            return False
        if self.end_inclusive:
            return True
        if other.start_inclusive:
            return True
        return False

    def is_addable_to(self, other):
        """Return True if the intervals overlap or touch, False otherwise."""
        if not isinstance(other, Interval):
            raise TypeError('other must be Interval.')
        if self.overlaps(other):
            return True
        if self.is_continued_by(other):
            return True
        if other.is_continued_by(self):
            return True
        return False

    def is_subtractable(self, other):
        """Return True if the intervals are subtractable, False, otherwise."""
        if self._overhangs_left(other):
            return True
        if other._overhangs_left(self):
            return True
        return False

    ############################################################
    # Interval-returning methods                               #
    ############################################################
    def intersect(self, other):
        """Return the intersection interval.

            The intersect exists if the intervals overlap.
        """
        if not isinstance(other, Interval):
            raise TypeError('other must be Interval.')
        if not self.overlaps(other):
            message = 'Intervals do not overlap: %s, %s' % (self, other)
            raise ValueError(message)

        # start
        if self.start == other.start:
            start = self.start
            start_inclusive = self.start_inclusive and other.start_inclusive
        elif self.start > other.start:
            start = self.start
            start_inclusive = self.start_inclusive
        else:
            start = other.start
            start_inclusive = other.start_inclusive

        # end
        if self.end == other.end:
            end = self.end
            end_inclusive = self.end_inclusive and other.end_inclusive
        elif self.end < other.end:
            end = self.end
            end_inclusive = self.end_inclusive
        else:
            end = other.end
            end_inclusive = other.end_inclusive

        return Interval(
                start=start, end=end,
                start_inclusive=start_inclusive, end_inclusive=end_inclusive,
                )

    def add(self, other):
        """Return the union interval."""
        if not isinstance(other, Interval):
            raise TypeError('other must be Interval.')
        if not self.is_addable_to(other):
            raise ValueError('Intervals not unitable: %s, %s' % (self, other))

        # start
        if self.start == other.start:
            start = self.start
            start_inclusive = self.start_inclusive or other.start_inclusive
        elif self.start < other.start:
            start = self.start
            start_inclusive = self.start_inclusive
        else:
            start = other.start
            start_inclusive = other.start_inclusive

        # end
        if self.end == other.end:
            end = self.end
            end_inclusive = self.end_inclusive or other.end_inclusive
        elif self.end > other.end:
            end = self.end
            end_inclusive = self.end_inclusive
        else:
            end = other.end
            end_inclusive = other.end_inclusive

        return Interval(
                start=start, end=end,
                start_inclusive=start_inclusive, end_inclusive=end_inclusive,
                )

    def subtract(self, other):
        """Return the difference interval."""
        if not isinstance(other, Interval):
            raise TypeError('other must be Interval.')
        if not self.overlaps(other):
            raise ValueError(
                'Intervals do not overlap: %s, %s' % (self, other)
                )

        if other._overhangs_left(self):
            return other.subtract(self)

        if not self._overhangs_left(other):
            raise ValueError('Result would not be and Interval.')

        # when this line is reached, self overhangs other on the left
        start = self.start
        start_inclusive = self.start_inclusive
        end = other.start
        end_inclusive = not other.start_inclusive
        return Interval(start, end, start_inclusive, end_inclusive)

    ############################################################
    # helpers                                                  #
    ############################################################
    def _overhangs_left(self, other):
        """Return True if `self` overhangs to the left of `other`."""
        if not isinstance(other, Interval):
            raise TypeError('other must be Interval.')

        if self == other:
            return False

        # start
        if self.start > other.start:
            return False
        if self.start == other.start:
            if (not self.start_inclusive) and other.start_inclusive:
                return False

        # end
        if self.end > other.end:
            return False
        if self.end == other.end:
            if self.end_inclusive and not other.end_inclusive:
                return False

        return True

    def _contains_interval(self, other):
        """Return a bool."""
        if self.start_inclusive or (not other.start_inclusive):
            lower_cond = self.start <= other.start
        else:
            lower_cond = self.start < other.start

        if self.end_inclusive or (not other.end_inclusive):
            upper_cond = other.end <= self.end
        else:
            upper_cond = other.end < self.end

        return (lower_cond and upper_cond)

    def _contains_datetime(self, time):
        """Return a bool.

            Check whether `time` is in the interval.

            Parameters
            ----------
            time : datetime.datetime

            Returns
            -------
            bool
        """
        if self.start_inclusive:
            lower_cond = self.start <= time
        else:
            lower_cond = self.start < time

        if self.end_inclusive:
            upper_cond = time <= self.end
        else:
            upper_cond = time < self.end

        return (lower_cond and upper_cond)

    def _contains_daytime_period(self, period):
        """Return a bool."""
        if self.length() >= dt.timedelta(days=1):
            return True
        speriod = DaytimePeriod.from_interval(self)
        return speriod.__contains__(period)

    def _contains_season(self, season):
        """Return a bool."""
        length = self.length()
        if length >= dt.timedelta(days=366):
            return True
        elif length > dt.timedelta(days=365):
            raise NotImplementedError('Not straightforward due to leap years.')

        s = Season.from_interval(self)
        return s.__contains__(season)

    def _ends_before(self, time):
        """Return a bool.

            Parameters
            ----------
            time : datetime.datetime

            Returns
            -------
            bool
                True if `self` ends before `time`, False otherwise.
        """
        if not isinstance(time, dt.datetime):
            raise TypeError('time must be datetime.datetime.')
        if self.end_inclusive:
            return self.end < time
        else:
            return self.end <= time

    def _starts_after(self, time):
        """Return a bool.

            Parameters
            ----------
            time : datetime.datetime

            Returns
            -------
            bool
                True if `self` start after `time`, False otherwise.
        """
        if not isinstance(time, dt.datetime):
            raise TypeError('time must be datetime.datetime.')
        if self.start_inclusive:
            return self.start > time
        else:
            return self.start >= time

    def _completely_earlier(self, other):
        """Return a bool."""
        if self.end_inclusive and other.start_inclusive:
            return self.end < other.start
        else:
            return self.end <= other.start


class DaytimePeriod(object):
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
        """Construct a Season from an Interval.

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
        if isinstance(other, dt.time):
          return self._contains_time(other)
        if isinstance(other, dt.datetime):
          return self._contains_time(other)
        if isinstance(other, DaytimePeriod):
          return self._contains_daytime_period(other)
        if isinstance(other, Interval):
            return self._contains_interval(other)
        raise TypeError('Cannot compare to %s' % other.__class__)

    def __eq__(self, other):
        """return a bool."""
        if not isinstance(other, DaytimePeriod):
            return False
        if other.start != self.start:
            return False
        if other.end != self.end:
            return False
        return True

    def contains(self, other):
        """Alias for backward compatibility."""
        return self.__contains__(other)

    def is_full_day(self):
        """Return a bool."""
        return self.length() == dt.timedelta(days=1)

    def length(self):
        """Return an instance of datetime.timedelta."""
        return self.end - self.start

    def overlaps(self, other):
        """Return a bool."""
        if not isinstance(other, DaytimePeriod):
            raise TypeError(
                    'Expected DaytimePeriod, got %s' % other.__class__
                    )
        if self.start in other:
            return True
        if self.end in other:
          return True
        if other.start in self:
          return True
        if other.end in self:
            return True
        return False

    ############################################################
    # helpers                                                  #
    ############################################################
    def _contains_daytime_period(self, other):
        """Return a bool."""
        if self.start > other.start:
            return False
        if self.end < other.end:
            return False
        return True

    def _contains_interval(self, interval):
        """Return a bool."""
        if interval.length() > self.length():
            return False
        if interval.end_inclusive and interval.end == self.end:
            return False

        period = DaytimePeriod.from_interval(interval)
        return self.__contains__(period)

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


class Season(object):
    """A section of the year cycle.

        This type is unaware of its absolute year number. It can extend beyond
        New Year, e. g. Nov 2nd to Jan 10th.

        Examples
        --------
        >>> import datetime as dt
        >>> lon = 6.9
        >>> lat = 50.9
        >>> now = dt.datetime.now()
        >>> beg = dt.datetime(1, 6, 21)
        >>> end = dt.datetime(1, 9, 23)

        # instantiate with beg and end:
        >>> season = Season(beg, end)
        >>> if season.contains(now):
        >>>     print 'It is summer!'
        >>> else:
        >>>     print 'It is not summer.'

        # alternatively, instantiate with months:
        >>> season = Season(months='DJF')
        >>> if season.contains(now):
        >>>     print 'It is meteorological winter.'
        >>> else:
        >>>     print 'It is not meteorological winter'

        # show length of season:
        >>> print season.length()

        # special case: beg == end
        >>> season = Season(beg, beg)
        >>> print season.length()
        365 days, 0:00:00
        >>> season = Season(beg, beg, allow_whole_year=False)
        >>> print season.length()
        0:00:00

        # the function is unaware of the absolute year, so pay attention in
        # leap years!
        >>> beg = dt.datetime(2016, 2, 1)
        >>> end = dt.datetime(2016, 3, 1)
        >>> season = Season(beg, end)
        >>> print season.length()
        28 days, 0:00:00

        # the function is microsecond-precise:
        >>> beg = dt.datetime(1, 6, 1)
        >>> end = dt.datetime(1, 6, 30, 12)
        >>> instant1 = end - dt.timedelta(microseconds=2)
        >>> instant2 = end + dt.timedelta(microseconds=2)
        >>> season = Season(beg, end)
        >>> print season.contains(instant1)
        True
        >>> print season.contains(instant2)
        False
       """
    def __init__(
            self,
            dt_start=dt.datetime(1, 1, 1),
            dt_end=dt.datetime(1, 1, 1),
            months='',
            allow_whole_year=True,
        ):
        """dt_start and dt_end must be dt.datetime or dt.date objects.

            The absolute year numbers of the input arguments are ignored. This
            means that calling this function with dt_start 1st Jan 1970 and
            dt_end

            30th Mar 1970 is equivalent to calling it with dt_start 1st Jan
            2014 and dt_end 30th Mar 1789.

            Parameters
            ----------
            dt_start : dt.datetime, optional
                lower boundary, inclusive
            dt_end : dt.datetime, optional
                upper boundary, exclusive
            months : str, optional
                This overrides `dt_start` and `dt_end`. Valid values are:
                1. three-char abbreviation of one month, such as
                    {'jan', 'feb', ...} or
                2. sequence of chars of consecutive month initials, such as
                    {'djf', 'mam', ...}, but also {'fm', 'jjaso', ...}
                    (first and last months are inclusive) or
                3. 'year' for the whole year or
                4. '' to disable and use `dt_start` and `dt_end` instead
                    (default)
                `months` is not case-sensitive
            allow_whole_year : bool, optional
                only applies if `dt_start` and `dt_end` are identical and, in
                this case, determines whether Season contains the whole year or
                nothing.

            Special cases
            -------------
            1. dt_start == dt_end
                * if allow_whole_year : Season contains the whole year.
                * else                : Season does not contain anything.

            2. February 29th
                If the date of dt_start and/or dt_end is Feb 29th, they are
                treated as Mar 1st 00:00.
        """
        ###################################################
        # INPUT CHECK                                     #
        ###################################################
        for d in (dt_start, dt_end):
            assert isinstance(d, dt.date)
        assert isinstance(allow_whole_year, bool)

        # months:
        if months.lower() == 'year':
            mon = 'jfmamjjasond'
        else:
            mon = months.lower()
        allowed_month_seqs = 'jfmamjjasondjfmamjjasond'
        allowed_months = (
                'jan', 'feb', 'mar', 'apr', 'may', 'jun',
                'jul', 'aug', 'sep', 'oct', 'nov', 'dec',
                )
        if not isinstance(mon, str):
            raise TypeError('months must be str.')
        if len(mon) == 1:
            msg = 'months must be a str of at least two month initial letters.'
            raise ValueError(msg)
        if len(mon) > 12:
            raise ValueError('months must not contain more than 12 letters.')
        if mon not in allowed_months:
            if mon[:3].lower() not in allowed_month_seqs:
                msg = 'Not a sequence of month initial letters: %s' % mon
                raise ValueError(msg)

        ###################################################
        # INITIALIZE                                      #
        ###################################################
        # self.start will be in year 1
        # self.end will be between 0 and 1 year later than self.end, and will
        # thus be in year 1 or 2
        if mon == '':
            start = year_one(dt_start)
            end   = year_one(dt_end)
        elif mon in allowed_month_seqs:
            first_month = allowed_month_seqs.find(mon) + 1
            last_month = (first_month + len(mon)) % 12
            if last_month == 0:
                last_month = 12
            start = dt.datetime(1, first_month, 1)
            end   = dt.datetime(1, last_month, 1)
        elif mon[:3].lower() in allowed_months:
            first_month = allowed_months.index(mon) + 1
            last_month  = (first_month + 1) % 12
            if last_month == 0:
                last_month = 12
            start = dt.datetime(1, first_month, 1)
            end   = dt.datetime(1, last_month, 1)


        # make sure that end is not earlier than start:
        if end < start:
            end = end.replace(year=2)
        if end == start and allow_whole_year:
            end = end.replace(year=2)

        self.start = start
        self.end   = end

    @classmethod
    def from_interval(cls, interval, allow_whole_year=True):
        """Construct a Season from an Interval."""
        bounds = interval.get_bounds()
        return cls(*bounds, allow_whole_year=allow_whole_year)

    @classmethod
    def from_tuple(cls, beg, end, allow_whole_year=True):
        """Shorthand constructer using datetime args as tuples.

            beg : tuple, length >= 2
                interpreted as month, day, [hour, minute, second, microsecond]
            end : tuple, length >= 2
                interpreted as month, day, [hour, minute, second, microsecond]
            allow_whole_day : bool or None
                if None, True is assumed for intervals of finite length
        """
        dt_beg = dt.datetime(1, *beg)
        dt_end = dt.datetime(1, *end)
        return cls(dt_beg, dt_end, allow_whole_year=allow_whole_year)

    def __repr__(self):
        """Return a str."""
        return 'Season ' + str(self)

    def __str__(self):
        """Return a str."""
        s = self.start.replace(year=1900)
        e = self.end.replace(year=1900)

        fmt = '%d %b'
        add = ''
        if s.hour != 0 or e.hour != 0:
            add = ' %H:%M'
        if s.second != 0 or e.second != 0:
            add = ' %H:%M:%S'
        if s.microsecond != 0 or e.microsecond != 0:
            add = ' %H:%M:%S.%f'
        fmt = fmt + add
        beg_str = s.strftime(fmt)
        end_str = e.strftime(fmt)

        if s != e:
            text = '[%s, %s)' % (beg_str, end_str)
        elif self.start == self.end:
            text = '(empty)'
        else:
            text = '(full year)'

        return text

    def __eq__(self, other):
        """return a bool."""
        if not isinstance(other, Season):
            return False
        if other.start != self.start:
            return False
        if other.end != self.end:
            return False
        return True

    def __contains__(self, other):
        """Return a bool or a list of such."""
        if isinstance(other, Season):
            return self._contains_season(other)
        if isinstance(other, Interval):
            return self._contains_interval(other)
        if isinstance(other, dt.date):
            return self._contains_datetime(other)
        raise TypeError('Cannot compare to %s' % other.__class__)

    def contains(self, other):
        """Alias for backward compatibility."""
        return self.__contains__(other)

    def months(self):
        """Return a string in the form of 'JJA' of 'Jun'.

            A month is considered to be 'covered' if at least 15 of its days
            are within the season.

            If the season covers more than one month, a sequence of capitals,
            such as 'JJA' (June-August) or 'NDJFMA' (November-April) is
            returned.

            If the season covers one month, than the first three characters of
            its name are returned.

            If the season covers less than 15 days in any month, '' is
            returned.
        """
        raise NotImplementedError('')

    def length(self):
        """Return an instance of datetime.timedelta."""
        return self.end - self.start

    ############################################################
    # helpers                                                  #
    ############################################################
    def _contains_season(self, other):
        """Return a bool."""
        if self.start > other.start:
            return False
        if self.end < other.end:
            return False
        return True

    def _contains_interval(self, interval):
        """Return a bool."""
        if interval.length() > self.length():
            return False
        if interval.end_inclusive and interval.end == self.end:
            return False

        season = Season.from_interval(interval)
        return self.__contains__(season)

    def _contains_datetime(self, time):
        """Check whether time is within Season or not.

            Parameters
            ----------
            time : dt.datetime or dt.date or Iterable of such

            Returns
            -------
            bool or list of such
        """
        ###################################################
        # INPUT CHECK                                     #
        ###################################################
        if not isinstance(time, dt.date):
            raise TypeError(
                    'Argument must be datetime.datetime or datetime.date.'
                    )

        # make time_yo be after self.start
        time_yo = year_one(time)
        if time_yo < self.start:
            time_yo = time_yo.replace(year=2)

        # check whether time_yo is earlier than self.end:
        if time_yo < self.end:
            return True
        else:
            return False


################################################################
# helper functions                                             #
################################################################
def year_one(d):
    """Return a datetime.datetime in year 1 CE.

        This is designed as a helper function for the Season class.

        Parameters
        ----------
        d : an instance of datetime.datetime or datetime.date

        Returns
        -------
        A datetime.datetime object. Date and time is be the same as of the
        input object, only the year is set to 1.

        Note
        ----
        Dates on 29th Feb will are converted to 1st Mar 00:00.

        Tested
        ------
        Intensively. Bug-free.
    """
    # input check
    if not isinstance(d, dt.date):
        raise TypeError('Expected dt.date or dt.datetime, got %s' % type(d))

    # Feb 29th
    if d.month == 2 and d.day == 29:
        return dt.datetime(1, 3, 1)

    # dt.datetime
    if isinstance(d, dt.datetime):
        return d.replace(year=1)

    # dt.date
    return dt.datetime.combine(d.replace(year=1), dt.time())

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
