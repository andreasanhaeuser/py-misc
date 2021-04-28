# standard modules
"""timedelta classes."""
import calendar
import datetime as dt

import numpy as np

class MonthTimeDelta():
    def __init__(self, number):
        assert isinstance(number, int)
        self.number = number
        self.days = 0

    def __add__(self, other):
        if not isinstance(other, dt.date):
            raise TypeError('Cannot handle %s' % type(other))

        return add_months(other, self.number)

    def __radd__(self, other):
        return self.__add__(other)

    def __neg__(self):
        return MonthTimeDelta(-self.number)

    def __sub__(self, other):
        neg = - self
        return neg.__add__(other)

    def __rsub__(self, other):
        return self.__sub__(other)

    def __eq__(self, other):
        if isinstance(other, MonthTimeDelta):
            if self.number != other.number:
                return False
            if self.days != other.days:
                return False
            return True

        if isinstance(other, YearTimeDelta):
            if self.number != other.number * 12:
                return False
            if self.days != other.days:
                return False
            return True

        if isinstance(other, dt.timedelta):
            other_days = get_total_days(other)
            if self.max_days() < other_days:
                return False
            if self.min_days() > other_days:
                return False
            raise ValueError(
                    'Comparison between %s and %s is amiguous.' 
                    % (str(self), str(other))
                    )

        return False

    def __lt__(self, other):
        if isinstance(other, MonthTimeDelta):
            if self.number < other.number:
                return True
            if self.number > other.number:
                return False
            return self.days < other.days

        if isinstance(other, YearTimeDelta):
            if self.number < other.number * 12:
                return True
            if self.number > other.number * 12:
                return False
            return self.days < other.days

        if isinstance(other, dt.timedelta):
            other_days = get_total_days(other)
            if self.max_days() < other_days:
                return True
            if self.min_days() > other_days:
                return False
            raise ValueError(
                    'Comparison between %s and %s is amiguous.' 
                    % (str(self), str(other))
                    )

        return False

    def __gt__(self, other):
        if isinstance(other, MonthTimeDelta):
            if self.number > other.number:
                return True
            if self.number < other.number:
                return False
            return self.days > other.days

        if isinstance(other, YearTimeDelta):
            if self.number > other.number * 12:
                return True
            if self.number < other.number * 12:
                return False
            return self.days > other.days

        if isinstance(other, dt.timedelta):
            other_days = get_total_days(other)
            if self.min_days() > other_days:
                return True
            if self.max_days() < other_days:
                return False
            raise ValueError(
                    'Comparison between %s and %s is amiguous.' 
                    % (str(self), str(other))
                    )

        return False

    def __le__(self, other):
        if isinstance(other, MonthTimeDelta):
            if self.number < other.number:
                return True
            if self.number > other.number:
                return False
            return self.days <= other.days

        if isinstance(other, YearTimeDelta):
            if self.number < other.number * 12:
                return True
            if self.number > other.number * 12:
                return False
            return self.days <= other.days

        if isinstance(other, dt.timedelta):
            other_days = get_total_days(other)
            if self.max_days() < other_days:
                return True
            if self.min_days() > other_days:
                return False
            raise ValueError(
                    'Comparison between %s and %s is amiguous.' 
                    % (str(self), str(other))
                    )

        return False

    def __ge__(self, other):
        if isinstance(other, MonthTimeDelta):
            if self.number > other.number:
                return True
            if self.number < other.number:
                return False
            return self.days >= other.days

        if isinstance(other, YearTimeDelta):
            if self.number > other.number * 12:
                return True
            if self.number < other.number * 12:
                return False
            return self.days >= other.days

        if isinstance(other, dt.timedelta):
            other_days = get_total_days(other)
            if self.min_days() > other_days:
                return True
            if self.max_days() < other_days:
                return False
            raise ValueError(
                    'Comparison between %s and %s is amiguous.' 
                    % (str(self), str(other))
                    )

        return False

    def min_days(self):
        """Return the minimum number of days that this can represent."""
        month_lengths = [28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31]
        days_years = 365.24 * (self.number // 12)
        days_months = np.sum(month_lengths[:self.number % 12])
        days_frac = days_years + days_months + self.days
        return int(np.floor(days_frac))

    def max_days(self):
        """Return the maximum number of days that this can represent."""
        month_lengths = [31, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 29]
        days_years = 365.24 * (self.number // 12)
        days_months = np.sum(month_lengths[:self.number % 12])
        days_frac = days_years + days_months + self.days
        return int(np.ceil(days_frac))

class YearTimeDelta():
    def __init__(self, number):
        assert isinstance(number, int)
        self.number = number
        self.days = 0

    def __add__(self, other):
        if not isinstance(other, dt.date):
            raise TypeError('Cannot handle %s' % type(other))

        return add_years(other, self.number)

    def __radd__(self, other):
        return self.__add__(other)

    def __neg__(self):
        return YearTimeDelta(-self.number)

    def __sub__(self, other):
        neg = - self
        return neg.__add__(other)

    def __rsub__(self, other):
        return self.__sub__(other)

    def __eq__(self, other):
        if isinstance(other, MonthTimeDelta):
            return other == self

        if isinstance(other, YearTimeDelta):
            if self.number != other.number:
                return False
            if self.days != other.days:
                return False
            return True

        if isinstance(other, dt.timedelta):
            other_days = get_total_days(other)
            if self.max_days() < other_days:
                return False
            if self.min_days() > other_days:
                return False
            raise ValueError(
                    'Comparison between %s and %s is amiguous.' 
                    % (str(self), str(other))
                    )

        return False

    def __lt__(self, other):
        if isinstance(other, MonthTimeDelta):
            return other > self

        if isinstance(other, YearTimeDelta):
            if self.number < other.number:
                return True
            if self.number > other.number:
                return False
            return self.days < other.days

        if isinstance(other, dt.timedelta):
            other_days = get_total_days(other)
            if self.max_days() < other_days:
                return True
            if self.min_days() > other_days:
                return False
            raise ValueError(
                    'Comparison between %s and %s is amiguous.' 
                    % (str(self), str(other))
                    )

        return False

    def __gt__(self, other):
        if isinstance(other, MonthTimeDelta):
            return other < self

        if isinstance(other, YearTimeDelta):
            if self.number > other.number:
                return True
            if self.number < other.number:
                return False
            return self.days > other.days

        if isinstance(other, dt.timedelta):
            other_days = get_total_days(other)
            if self.min_days() > other_days:
                return True
            if self.max_days() < other_days:
                return False
            raise ValueError(
                    'Comparison between %s and %s is amiguous.' 
                    % (str(self), str(other))
                    )

        return False

    def __le__(self, other):
        if isinstance(other, MonthTimeDelta):
            return other >= self

        if isinstance(other, YearTimeDelta):
            if self.number < other.number:
                return True
            if self.number > other.number:
                return False
            return self.days <= other.days

        if isinstance(other, dt.timedelta):
            other_days = get_total_days(other)
            if self.max_days() < other_days:
                return True
            if self.min_days() > other_days:
                return False
            raise ValueError(
                    'Comparison between %s and %s is amiguous.' 
                    % (str(self), str(other))
                    )

        return False

    def __ge__(self, other):
        if isinstance(other, MonthTimeDelta):
            return other <= self

        if isinstance(other, YearTimeDelta):
            if self.number > other.number:
                return True
            if self.number < other.number:
                return False
            return self.days >= other.days

        if isinstance(other, dt.timedelta):
            other_days = get_total_days(other)
            if self.min_days() > other_days:
                return True
            if self.max_days() < other_days:
                return False
            raise ValueError(
                    'Comparison between %s and %s is amiguous.' 
                    % (str(self), str(other))
                    )

        return False

    def min_days(self):
        """Return the minimum number of days that this can represent."""
        # This can be refined, as there are no more than 1 consecutive months
        # with 28 days
        days_fractional = self.number * 365.24 + self.days
        return int(np.floor(days_fractional))

    def max_days(self):
        """Return the maximum number of days that this can represent."""
        # This can be refined, as there are at max 2 consecutive months with 31
        # days
        days_fractional = self.number * 365.24 + self.days
        return int(np.ceil(days_fractional))

class EnhancedTimeDelta(object):
    def __init__(
            self, days=0, seconds=0, microseconds=0, minutes=0, hours=0,
            weeks=0, years=0, months=0,
            ):
        raise NotImplementedError()
        if years == 0 and months == 0:
            return super().__init__(
                    days=days, seconds=seconds, microseconds=microseconds,
                    minutes=minutes, hours=hours, weeks=weeks,
                    )

###################################################
# ITERATION                                       #
###################################################
def next_month(year, month):
    """Return (year, month) of the following month.
    
        Parameters
        ----------
        year : int
        month : int
        
        Returns
        -------
        year_next : int
        month_next : int
    """
    if month < 12:
        month += 1
    else:
        month = 1
        year += 1
    return year, month

def previous_month(year, month):
    """Return (year, month) of the previous month.
    
        Parameters
        ----------
        year : int
        month : int
        
        Returns
        -------
        year_prev : int
        month_prev : int
    """
    if month >= 1:
        month -= 1
    else:
        month = 12
        year -= 1
    return year, month

def next_day(year, month, day):
    """Return (year, month, day) of the following day.
    
        Parameters
        ----------
        year : int
        month : int
        day : int
        
        Returns
        -------
        year_next : int
        month_next : int
        day_next : int
    """
    thisday = dt.datetime(year, month, day)
    nextday = thisday + dt.timedelta(days=1)
    y = nextday.year
    m = nextday.month
    d = nextday.day
    return y, m, d

def add_months(time, number=1):
    """Increase by number of months.
        
        Parameters
        ----------
        time : datetime.date
        number : int (pos or neg)
            the number of months to add (subtract if negative)

        Returns
        -------
        datetime.datetime

        Raises
        ------
        ValueError
            If the day of the input does not exist in the output (29th...31st)
    """
    assert isinstance(time, dt.date)
    assert isinstance(number, int)
    number = int(number)

    year = time.year
    month = time.month
    day = time.day
    if number >= 0:
        for n in range(number):
            year, month = next_month(year, month)
    else:
        for n in range(abs(number)):
            year, month = previous_month(year, month)

    days_in_month = calendar.mdays[month]
    if calendar.isleap(year) and month == 2:
        days_in_month = 29

    if day > days_in_month:
        whatis = 'Invalid date: %04i-%02i-%02i.' % (year, month, day)
        day = days_in_month
        willdo = 'Setting to:   %04i-%02i-%02i.' % (year, month, day)
        message = ' '.join(whatis, willdo)
        warnings.warn(message)

    return time.replace(year=year, month=month, day=day)

def add_years(time, number=1):
    """Increase by number of years."""
    assert isinstance(time, dt.date)
    assert isinstance(number, int)
    number = int(number)

    year = time.year + number
    month = time.month
    day = time.day

    if calendar.isleap(year) and month == 2 and day == 29:
        whatis = '29th Feb in non-leap year %i.' % year
        willdo = 'Setting to 28th Feb.'
        message = ' '.join(whatis, willdo)
        warnings.warn(message)

        day = 28

    return time.replace(year=year, day=day)


################################################################
# helpers                                                      #
################################################################
def get_total_days(x):
    """Return total days of timedelta."""
    if not isinstance(x, dt.timedelta):
        raise TypeError(str(type(x)))

    return x.total_seconds() / 86400.
