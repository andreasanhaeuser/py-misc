#!/usr/bin/python3
"""Testing suite for datetime_utils."""

# standard
import datetime as dt
import unittest

# to be tested
import misc.date_time.utils as du
from misc.date_time.utils import Interval, DaytimePeriod, Season

class MonthsInInterval(unittest.TestCase):
    def test_months_in_interval(self):
        time_beg = dt.datetime(2000, 7, 2, 0, 3, 1, 323)
        time_end = dt.datetime(2001, 4, 5, 0, 0, 1, 11)
        time_end_flat = dt.datetime(2001, 4, 1)
        interval = Interval(time_beg, time_end)
        interval_flat = Interval(time_beg, time_end_flat)
        interval_flat_incl = Interval(
                time_beg, time_end_flat, end_inclusive=True)

        
        last_beg = dt.datetime(2001, 4, 1)
        last_end = dt.datetime(2001, 5, 1)
        last_month = Interval(last_beg, last_end)

        months = du.months_in_interval(interval)
        self.assertEqual(months[-1], last_month)

        months_flat = du.months_in_interval(interval_flat)
        self.assertNotEqual(months_flat[-1], last_month)

        months_flat_incl = du.months_in_interval(interval_flat_incl)
        self.assertEqual(months_flat_incl[-1], last_month)

class TimeDelta(unittest.TestCase):
    def test_string_to_timedelta(self):
        print('*' * 40)
        print('STRING TO TIMDELTA')
        in_out = (
                ('1 d', dt.timedelta(1)),
                ('_1 d_', dt.timedelta(1)),
                (' 1 d ', dt.timedelta(1)),
                (' 1_d ', dt.timedelta(1)),
                ('1 day', dt.timedelta(1)),
                ('2 days', dt.timedelta(2)),
                ('2_days', dt.timedelta(2)),
                ('2', dt.timedelta(2)),
                (['2', 'days'], dt.timedelta(2)),
                ('2 min', dt.timedelta(minutes=2)),
                ('4 s', dt.timedelta(seconds=4)),
                ('3_hours', dt.timedelta(hours=3)),
                ('1_hour', dt.timedelta(hours=1)),
                ('-1_hour', dt.timedelta(hours=-1)),
                )

        for pair in in_out:
            string, expected = pair
            result = du.str_to_timedelta(string)
            self.assertEqual(result, expected)
            print(str(string).ljust(15), '->', str(result).rjust(16))
        print('*' * 40)

class DatetimeRange(unittest.TestCase):
    def test_monthly(self):
        beg = dt.datetime(2020, 1, 3)
        end = dt.datetime(2021, 3, 5)
        inc = du.MonthTimeDelta(2)
        times = du.datetime_range(beg, end, inc)
        expected = [
                dt.datetime(2020, 1, 3),
                dt.datetime(2020, 3, 3),
                dt.datetime(2020, 5, 3),
                dt.datetime(2020, 7, 3),
                dt.datetime(2020, 9, 3),
                dt.datetime(2020, 11, 3),
                dt.datetime(2021, 1, 3),
                dt.datetime(2021, 3, 3),
                ]
        N = len(times)
        self.assertEqual(N, len(expected))
        for n in range(N):
            self.assertEqual(times[n], expected[n])

if __name__ == '__main__':
    unittest.main()
