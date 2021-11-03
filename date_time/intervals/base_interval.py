"""An abstract parent class for Interval, Season and DaytimePeriod."""

# standard
from collections import Iterable

# PyPI
import numpy as np

class BaseInterval(object):
    def length(self):
        """Return an instance of datetime.timedelta."""
        return self.end - self.start

    ############################################################
    # contains                                                 #
    ############################################################
    def contains_iterable(self, other):
        """Return a ndarray of bool."""
        assert isinstance(other, Iterable)

        idx = np.zeros_like(other, dtype=bool)
        for n, element in enumerate(other):
            idx[n] = self.__contains__(element)
        return idx

    ############################################################
    # length in human readable units                           #
    ############################################################
    def total_seconds(self):
        """Return interval duration in seconds as float."""
        return self.length().total_seconds()

    def total_minutes(self):
        """Return interval duration in minutes as float."""
        return self.total_seconds() / 60.

    def total_hours(self):
        """Return interval duration in hours as float."""
        return self.total_minutes() / 60.

    def total_days(self):
        """Return interval duration in days as float."""
        return self.total_hours() / 24.

    def total_weeks(self):
        """Return interval duration in weeks as float."""
        return self.total_days() / 7.


class CyclicInterval(BaseInterval):
    """An abstract parent class for Season and DaytimePeriod."""
    def __eq__(self, other):
        """return a bool."""
        self_class = type(self)
        if not isinstance(other, self_class):
            return False
        if other.start != self.start:
            return False
        if other.end != self.end:
            return False
        return True

    def contains(self, other):
        """Alias for backward compatibility."""
        return self.__contains__(other)

    def overlaps(self, other):
        """Return a bool."""
        self_class = type(self)

        if not isinstance(other, self_class):
            message = 'Expected %s, got %s' % (self_class, type(other))
            raise TypeError(message)
        if self.start in other:
            return True
        if self.end in other:
          return True
        if other.start in self:
          return True
        if other.end in self:
            return True
        return False
