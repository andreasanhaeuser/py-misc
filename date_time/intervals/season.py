"""A season."""
# standard modules
import datetime as dt
from collections import Iterable

# PyPI
import numpy as np

# local
from .base_interval import CyclicInterval
from .interval import Interval
from .daytime_period import DaytimePeriod


class Season(CyclicInterval):
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
            interval=None,
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
        # ---------------------------------------
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
        # ---------------------------------------

        ###################################################
        # INITIALIZE                                      #
        ###################################################
        # self.start will be in year 1
        # self.end will be between 0 and 1 year later than self.end, and will
        # thus be in year 1 or 2

        # init from months
        # ------------------------------------------------
        if len(mon) >= 2 and mon in allowed_month_seqs:
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

        # init from start, end
        # ------------------------------------------------
        else:
            start = year_one(dt_start)
            end   = year_one(dt_end)
        # ------------------------------------------------

        # make sure that end is not earlier than start:
        # ------------------------------------------------
        if end < start:
            end = end.replace(year=2)
        if end == start and allow_whole_year:
            end = end.replace(year=2)
        # ------------------------------------------------

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
        return self.nice_string()

    def nice_string(self, month_fmt='%B'):
        """Return a nice human-readable string if possible."""
        start = self.start
        end = self.end

        # empty
        # ==================================================
        if start == end:
            return 'empty season'
        # ==================================================

        # full year
        # ==================================================
        if start.replace(year=end.year) == end:
            return 'full year'
        # ==================================================

        # one full month
        # ==================================================
        if self.is_one_calendar_month():
            return self.start.strftime('%B')
        # ==================================================

        # several full months
        # ==================================================
        if self.is_calendar_months():
            end_mod = end - dt.timedelta(1)
            start_str = start.strftime(month_fmt)
            end_str = end_mod.strftime(month_fmt)
            return start_str + ' - ' + end_str
        # ==================================================

        # Default
        # ==================================================
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
        return '[%s, %s)' % (beg_str, end_str)
        # ==================================================

    def __contains__(self, other):
        """Return a bool or a list of such."""
        # other is iterable
        # -----------------------------------------------
        if isinstance(other, Iterable):
            return self.contains_iterable(other)
        # -----------------------------------------------

        # select comparator
        # -----------------------------------------------
        if isinstance(other, dt.date):
            return self._contains_datetime(other)
        if isinstance(other, Interval):
            return self._contains_interval(other)
        if isinstance(other, Season):
            return self._contains_season(other)
        if isinstance(other, DaytimePeriod):
            return self._contains_daytime_period(other)
        # -----------------------------------------------

        # try inverse
        # -----------------------------------------------
        method = 'is_contained_in'
        if hasattr(other, method) and callable(getattr(other, method)):
            return other.is_contained_in(self)
        # -----------------------------------------------

        raise TypeError('Cannot compare to %s' % other.__class__)

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

    ############################################################
    # special intervals -- check                               #
    ############################################################
    def is_one_calendar_month(self):
        """Return True if its exactly one calendar month, False otherwise."""
        interval = Interval(self.start, self.end)
        return interval.is_one_calendar_month()

    def is_calendar_months(self):
        """Return True if its exactly calendar months, False otherwise."""
        interval = Interval(self.start, self.end)
        return interval.is_calendar_months()

    ############################################################
    # helpers: contains                                        #
    ############################################################
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

    def _contains_interval(self, interval):
        """Return a bool."""
        if interval.length() > self.length():
            return False
        if interval.end_inclusive and interval.end == self.end:
            return False

        season = Season.from_interval(interval)
        return self.__contains__(season)

    def _contains_season(self, other):
        """Return a bool."""
        if self.start > other.start:
            return False
        if self.end < other.end:
            return False
        return True

    def _contains_daytime_period(self, daytime_period):
        raise NotImplementedError()

    ############################################################
    # helpers: is_contained_in                                 #
    ############################################################
    def is_contained_in(self, other):
        """Return a bool."""
        if isinstance(other, Interval):
            return self.is_contained_in_interval(other)
        if isinstance(other, DaytimePeriod):
            return self.is_contained_in_daytime_period(other)
        if isinstance(other, Season):
            return other.__contains__(self)
        raise TypeError(type(other))

    def is_contained_in_interval(self, interval):
        length = interval.length()
        if length >= dt.timedelta(days=366):
            return True
        elif length > dt.timedelta(days=365):
            raise NotImplementedError('Not straightforward due to leap years.')

        season = Season.from_interval(interval)
        return season.__contains__(self)

    def is_contained_in_daytime_period(self, daytime_period):
        raise NotImplementedError()


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
