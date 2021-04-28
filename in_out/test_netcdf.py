#!/usr/bin/python3
"""Testing suite for datetime_utils."""

# standard
import unittest
import warnings

# PyPI
import numpy as np
import matplotlib.pyplot as plt

# misc
from misc.date_time.timer import Timer

# to be tested
import misc.in_out.netcdf as nc


_filename_in = 'test_files/file_in.nc'
_filename_out = 'test_files/file_out.nc'

class SelectVariables(unittest.TestCase):
    def test_select(self):
        filename_in = _filename_in
        filename_out = _filename_out
        varnames = ('time', 'qa_value')
        nc.select_variables(filename_in, filename_out, varnames)

        print(
                'Selected %s from\n%s -> %s.\nCheck manually.'
                % (varnames, filename_in, filename_out)
                )


class Bounds(unittest.TestCase):
    def test_bounds(self):
        filename_in = _filename_in
        filename_out = 'test_bounds.png'
        cmap = plt.get_cmap('gist_ncar')
        bounds = (124, 35, 128, 38)
        xlo, ylo, xhi, yhi = bounds
        xcorners = (xlo, xhi, xhi, xlo, xlo)
        ycorners = (ylo, ylo, yhi, yhi, ylo)

        # load
        timer = Timer().start()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            section, meta = nc.read_file(filename_in, bounds=bounds)
        timer.show().reset()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            full, meta = nc.read_file(filename_in)
        timer.show().stop()

        plt.close()
        plt.figure(figsize=(16, 9))

        x = full['longitude']
        y = full['latitude']
        z = full['nitrogendioxide_tropospheric_column']
        vmin = 0
        vmax = np.nanmax(z)
        plt.subplot(1, 2, 1)
        plt.pcolormesh(x, y, z, cmap=cmap, vmin=vmin, vmax=vmax)
        plt.plot(xcorners, ycorners, 'k-')
        plt.plot(xcorners, ycorners, 'w--')
        plt.axis('equal')

        plt.subplot(1, 2, 2)
        x = section['longitude']
        y = section['latitude']
        z = section['nitrogendioxide_tropospheric_column']
        plt.pcolormesh(x, y, z, cmap=cmap, vmin=vmin, vmax=vmax)
        plt.plot(xcorners, ycorners, 'k-')
        plt.plot(xcorners, ycorners, 'w--')
        plt.axis('equal')

        plt.savefig(filename_out)
        plt.close()
        print('Check manually:', filename_out)


if __name__ == '__main__':
    unittest.main()
