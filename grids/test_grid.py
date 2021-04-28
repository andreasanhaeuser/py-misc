#!/usr/bin/python3
"""Testing suite for grid module."""

# standard modules
import unittest

# PyPI modules
import numpy as np
import matplotlib.pyplot as plt

# module to be tested
from misc.grids import grid

class Grid(unittest.TestCase):
    def setUp(self):
        self.x0 = -4.3
        self.y0 = 1.7
        self.size_x = 8
        self.size_y = 6
        self.dx = 2
        self.dy = 3

    def test_coarse(self):
        self.title = 'coarse'

    def test_fine(self):
        self.dx = 0.2
        self.dy = 0.3
        self.title = 'fine'

    def test_exclude(self):
        self.size_x = 8 - 0.0001
        self.title = 'exlude x'

    def tearDown(self):
        x, y = grid.get_regular_2d_grid(
                self.x0, 
                self.y0, 
                self.size_x, 
                self.size_y, 
                self.dx, 
                self.dy,
                )
        xm, ym = np.meshgrid(x, y)
        x0 = self.x0
        y0 = self.y0
        size_x = self.size_x
        size_y = self.size_y
        xlo = x0 - size_x / 2
        xhi = x0 + size_x / 2
        ylo = y0 - size_y / 2
        yhi = y0 + size_y / 2
        plt.plot(
                [xlo, xhi, xhi, xlo, xlo], 
                [ylo, ylo, yhi, yhi, ylo],
                'b-',
                )
        plt.plot([x0], [y0], 'bo')
        plt.plot(xm, ym, 'r.')
        plt.title(self.title)
        plt.grid()
        ax = plt.gca()
        ax.axis('equal')
        plt.show()


if __name__ == '__main__':
    unittest.main()
