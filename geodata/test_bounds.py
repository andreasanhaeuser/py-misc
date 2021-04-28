#!/usr/bin/python3
"""A test suite to the bounds module."""

# standard modules
import unittest

# PyPI
import geopandas as gp

# the module to be tested
from misc.geodata import bounds


class Basic(unittest.TestCase):
    def setUp(self):
        self.bounds = (10, -60, 20, -40)

    def test_gdf(self):
        gdf = bounds.convert_bounds_to_gdf(self.bounds)
        self.assertIsInstance(gdf, gp.GeoDataFrame)
        tpl = bounds.convert_bounds_to_tuple(gdf)
        self.assertIsInstance(tpl, tuple)
        self.assertEqual(tpl, self.bounds)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
