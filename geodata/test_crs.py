#!/usr/bin/python3
"""A test suite to the crs module."""

# standard modules
import unittest

# the module to be tested
from misc.geodata import crs

class UTMTransform(unittest.TestCase):
    def setUp(self):
        self.units = 'km'
        self.lon = 3
        self.lat = 10
        self.x = 500
        self.y = self.lat/90 * 10000
        self.crs = crs.get_utm_proj4string(self.lon, self.lat, self.units)
        self.delta_x = 20
        self.delta_y = 20
        self.delta_lon = 0.1
        self.delta_lat = 0.1

    def test_lonlat_to_utm(self):
        x, y = crs.lonlat_to_utm(self.lon, self.lat, crs=self.crs)
        self.assertAlmostEqual(x, self.x, delta=self.delta_x)
        self.assertAlmostEqual(y, self.y, delta=self.delta_y)

    def test_utm_to_lonlat(self):
        lon, lat = crs.utm_to_lonlat(self.x, self.y, crs=self.crs)
        self.assertAlmostEqual(lon, self.lon, delta=self.delta_lon)
        self.assertAlmostEqual(lat, self.lat, delta=self.delta_lat)

    def tearDown(self):
        pass


class UTMProj4String(unittest.TestCase):
    def setUp(self):
        self.fractional = False
        self.units = 'm'

    def test_NE(self):
        self.lon = 0.2
        self.lat = 0.2
        self.units = 'km'
        self.exp = '+proj=utm +datum=WGS84 +no_defs +zone=31 +north +units=km'

    def test_SE(self):
        self.lon = 11.8
        self.lat = -88
        self.units = 'km'
        self.exp = '+proj=utm +datum=WGS84 +no_defs +zone=32 +south +units=km'

    def test_NW(self):
        self.lon = -0.1
        self.lat = 88
        self.units = 'km'
        self.exp = '+proj=utm +datum=WGS84 +no_defs +zone=30 +north +units=km'

    def test_SW(self):
        self.lon = -179
        self.lat = -88
        self.units = 'km'
        self.exp = '+proj=utm +datum=WGS84 +no_defs +zone=1 +south +units=km'

    def test_beyond_180_east(self):
        self.lon = 184
        self.lat = -4
        self.units = 'm'
        self.exp = '+proj=utm +datum=WGS84 +no_defs +zone=1 +south +units=m'

    def test_fractional(self):
        self.fractional = True
        self.lon = -6
        self.lat = -4
        self.exp = '+proj=utm +datum=WGS84 +no_defs +zone=29.5 +south +units=m'

    def tearDown(self):
        args = (self.lon, self.lat, self.units, self.fractional)
        result = crs.get_utm_proj4string(*args)
        print(result)
        self.assertEqual(result, self.exp)


if __name__ == '__main__':
    unittest.main()
    from misc.geodata import crs as crs_utils

    units = 'km'
    lon = 3
    lat = 10
    x = 500
    y = lat/90 * 10000
    crs = crs_utils.get_utm_proj4string(lon, lat, units)
    delta_x = 20
    delta_y = 20
    delta_lon = 0.1
    delta_lat = 0.1

    x, y = crs_utils.lonlat_to_utm(lon, lat, crs=crs)

    lon, lat = crs_utils.utm_to_lonlat(x, y, crs=crs)
