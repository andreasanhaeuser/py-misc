"""A time interval."""

# standard modules
import datetime as dt
from collections import Iterable

# PyPI
import numpy as np

# local
from .base_interval import BaseInterval


class Interval(BaseInterval):
    """A time interval.

        Each bound may be inclusive or exclusive. Default is inclusive for
        lower bound and exclusive for upper bound. Intervals in reverse
        direction are not allowed.

        History
        -------
        2021-05-19 (AA):  Allow dt.date in __contains__()
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

        # other is iterable
        # -----------------------------------------------
        if isinstance(other, Iterable):
            return self.contains_iterable(other)
        # -----------------------------------------------

        # select comparator
        # -----------------------------------------------
        if isinstance(other, Interval):
            return self._contains_interval(other)

        if isinstance(other, dt.datetime):
            return self._contains_datetime(other)

        if isinstance(other, dt.date):
            # date -> datetime
            other = dt.datetime.combine(other, dt.time())
            return self._contains_datetime(other)
        # -----------------------------------------------

        # try inverse
        # -----------------------------------------------
        method = 'is_contained_in'
        if hasattr(other, method) and callable(getattr(other, method)):
            return other.is_contained_in(self)
        # -----------------------------------------------

        raise TypeError('Cannot compare to %s' % type(other))

    def center(self):
        """Return the center of the interval as datetime.time."""
        return self.start + 0.5 * self.length()

    def get_bounds(self, fmt=None):
        """Return a pair of datetime.datetime or str.

            Parameters
            ----------
            fmt : str or None
                string conversion format
                None: return datedate.datetime objects

            Returns
            -------
            (datetime.datetime, datetime.datetime) or (str, str)
        """
        if fmt is None:
            return (self.start, self.end)

        if isinstance(fmt, str):
            start = self.start.strftime(fmt)
            end = self.end.strftime(fmt)
            return start, end

        raise TypeError('fmt must be None or str, got %s' % type(fmt))

    ############################################################
    # timedelta returning methods                              #
    ############################################################
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

    def is_one_calendar_month(self):
        """Return True if its exactly one calendar month, False otherwise."""
        start = self.start
        end = self.end

        if not self.is_calendar_months():
            return False

        if start.year == end.year and start.month + 1 == end.month:
            return True

        if start.year+1 == end.year and start.month == 12 and end.month == 1:
            return True

        return False

    def is_calendar_months(self):
        """Return True if its exactly calendar months, False otherwise."""
        start = self.start
        end = self.end

        # Is `start` exactly at the beginning of a month?
        if start != dt.datetime(start.year, start.month, 1):
            return False

        # Is `end` exactly at the beginning of a month?
        if end != dt.datetime(end.year, end.month, 1):
            return False

        return True

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
