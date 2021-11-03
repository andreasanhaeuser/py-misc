#!/usr/bin/python3
"""Testing suite for datetime_intervals."""

# standard
import datetime as dt
import unittest

# modules to be tested
import misc.date_time.utils as du
# from misc.date_time.utils import Interval, DaytimePeriod, Season

Interval = du.Interval
DaytimePeriod = du.DaytimePeriod
Season = du.Season

class Comparisons(unittest.TestCase):
    def setUp(self):
        """Prepare some intervals, daytime periods and seasons."""
        tbeg = dt.datetime(2018, 6, 1)
        tend = dt.datetime(2018, 6, 2)
        tinc = dt.timedelta(hours=1)
        times = du.datetime_range(tbeg, tend, tinc)
        self.times = times

        # intervals
        self.i_inner = Interval(times[4], times[5])
        self.i_outer = Interval(times[2], times[7])
        self.i_early = Interval(times[2], times[5])
        self.i_late = Interval(times[4], times[7])
        self.i_before = Interval(times[0], times[2])

        # daytime periods
        self.p_inner = DaytimePeriod.from_interval(self.i_inner)
        self.p_outer = DaytimePeriod.from_interval(self.i_outer)
        self.p_early = DaytimePeriod.from_interval(self.i_early)
        self.p_late = DaytimePeriod.from_interval(self.i_late)
        self.p_before = DaytimePeriod.from_interval(self.i_before)

        # seasons
        self.s_inner = Season.from_interval(self.i_inner)
        self.s_outer = Season.from_interval(self.i_outer)
        self.s_early = Season.from_interval(self.i_early)
        self.s_late = Season.from_interval(self.i_late)
        self.s_before = Season.from_interval(self.i_before)

    def test_interval_in_interval(self):
        self.assertIn(self.i_inner, self.i_outer)
        self.assertIn(self.i_inner, self.i_early)
        self.assertIn(self.i_inner, self.i_late)
        self.assertNotIn(self.i_outer, self.i_inner)
        self.assertNotIn(self.i_early, self.i_late)

    def test_interval_overlaps_interval(self):
        self.assertTrue(self.i_early.overlaps(self.i_late))
        self.assertTrue(self.i_late.overlaps(self.i_early))
        self.assertTrue(self.i_inner.overlaps(self.i_outer))
        self.assertTrue(self.i_outer.overlaps(self.i_inner))
        self.assertFalse(self.i_before.overlaps(self.i_early))
        self.assertFalse(self.i_before.overlaps(self.i_late))

        times = self.times

        early_open = Interval(times[0], times[1], False, False)
        early_closed = Interval(times[0], times[1], True, True)
        early_left_open = Interval(times[0], times[1], False, True)
        early_right_open = Interval(times[0], times[1], True, False)

        late_open = Interval(times[1], times[2], False, False)
        late_closed = Interval(times[1], times[2], True, True)
        late_left_open = Interval(times[1], times[2], False, True)
        late_right_open = Interval(times[1], times[2], True, False)

        self.assertTrue(early_closed.overlaps(late_closed))
        self.assertTrue(early_closed.overlaps(late_right_open))
        self.assertTrue(early_left_open.overlaps(late_closed))
        self.assertTrue(early_left_open.overlaps(late_right_open))

        self.assertFalse(early_open.overlaps(late_open))
        self.assertFalse(early_open.overlaps(late_closed))
        self.assertFalse(early_open.overlaps(late_left_open))
        self.assertFalse(early_open.overlaps(late_right_open))
        self.assertFalse(early_closed.overlaps(late_open))
        self.assertFalse(early_closed.overlaps(late_left_open))
        self.assertFalse(early_left_open.overlaps(late_open))
        self.assertFalse(early_left_open.overlaps(late_left_open))
        self.assertFalse(early_right_open.overlaps(late_open))
        self.assertFalse(early_right_open.overlaps(late_closed))
        self.assertFalse(early_right_open.overlaps(late_left_open))
        self.assertFalse(early_right_open.overlaps(late_right_open))

    def test_interval_in_daytime_period(self):
        self.assertIn(self.i_inner, self.p_inner)
        self.assertIn(self.i_inner, self.p_outer)
        self.assertIn(self.i_inner, self.p_early)
        self.assertIn(self.i_inner, self.p_late)
        self.assertNotIn(self.i_early, self.p_late)
        self.assertNotIn(self.i_outer, self.p_inner)
        self.assertIn(self.i_early, self.p_outer)
        self.assertIn(self.i_late, self.p_outer)

    def test_daytime_period_in_interval(self):
        self.assertIn(self.p_inner, self.i_inner)
        self.assertIn(self.p_inner, self.i_outer)
        self.assertIn(self.p_early, self.i_outer)
        self.assertIn(self.p_late, self.i_outer)
        self.assertNotIn(self.p_outer, self.i_inner)

    def test_daytime_period_in_daytime_period(self):
        self.assertIn(self.p_inner, self.p_inner)
        self.assertIn(self.p_inner, self.p_outer)
        self.assertIn(self.p_early, self.p_outer)
        self.assertIn(self.p_late, self.p_outer)
        self.assertNotIn(self.p_outer, self.p_inner)

    def test_interval_in_season(self):
        self.assertIn(self.i_inner, self.s_inner)
        self.assertIn(self.i_inner, self.s_outer)
        self.assertIn(self.i_inner, self.s_early)
        self.assertIn(self.i_inner, self.s_late)
        self.assertNotIn(self.i_early, self.s_late)
        self.assertNotIn(self.i_outer, self.s_inner)
        self.assertIn(self.i_early, self.s_outer)
        self.assertIn(self.i_late, self.s_outer)

    def test_season_in_season(self):
        self.assertIn(self.s_inner, self.s_inner)
        self.assertIn(self.s_inner, self.s_outer)
        self.assertIn(self.s_early, self.s_outer)
        self.assertIn(self.s_late, self.s_outer)
        self.assertNotIn(self.s_outer, self.s_inner)

    def test_season_in_interval(self):
        self.assertIn(self.s_inner, self.i_inner)
        self.assertIn(self.s_inner, self.i_outer)
        self.assertIn(self.s_early, self.i_outer)
        self.assertIn(self.s_late, self.i_outer)
        self.assertNotIn(self.s_outer, self.i_inner)

    def test_compare_intervals(self):
        self.assertTrue(self.i_before < self.i_early)
        self.assertFalse(self.i_early < self.i_late)
        self.assertFalse(self.i_early > self.i_late)
        self.assertFalse(self.i_inner < self.i_outer)

class SetOperations(unittest.TestCase):
    def setUp(self):
        """Prepare some intervals, daytime periods and seasons."""
        tbeg = dt.datetime(2018, 6, 1)
        tend = dt.datetime(2018, 6, 2)
        tinc = dt.timedelta(hours=1)
        times = du.datetime_range(tbeg, tend, tinc)
        self.times = times

        # intervals
        self.i_inner = Interval(times[4], times[5])
        self.i_outer = Interval(times[2], times[7])
        self.i_early = Interval(times[2], times[5])
        self.i_late = Interval(times[4], times[7])
        self.i_before = Interval(times[0], times[2])

        # daytime periods
        self.p_inner = DaytimePeriod.from_interval(self.i_inner)
        self.p_outer = DaytimePeriod.from_interval(self.i_outer)
        self.p_early = DaytimePeriod.from_interval(self.i_early)
        self.p_late = DaytimePeriod.from_interval(self.i_late)
        self.p_before = DaytimePeriod.from_interval(self.i_before)

        # seasons
        self.s_inner = Season.from_interval(self.i_inner)
        self.s_outer = Season.from_interval(self.i_outer)
        self.s_early = Season.from_interval(self.i_early)
        self.s_late = Season.from_interval(self.i_late)
        self.s_before = Season.from_interval(self.i_before)


if __name__ == '__main__':
    unittest.main()
