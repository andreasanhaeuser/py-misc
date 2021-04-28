#!/usr/bin/python3
"""Testing suite for chronometer.Chronometer."""

# standard
import unittest

# PyPI
import numpy as np
import matplotlib.pyplot as plt

# main
from misc.math.linreg import LinearRegression

class Graphically(unittest.TestCase):
    def test_graphically(self):
        #======== SETUP ========#
        x_center = 1
        x_spread = 1.0
        m_true = 1.3
        t_true = -0.1
        x_spread_noise = 0
        y_spread_noise = 0.25
        N = 10

        xmin = -1
        xmax = 10
        xinc = 1e-3
        #=======================#

        # create ideal data set
        x_ideal = np.random.normal(x_center, x_spread, N)
        y_ideal = m_true * x_ideal + t_true

        # add noise:
        x_noise = np.random.normal(0, x_spread_noise, N)
        y_noise = np.random.normal(0, y_spread_noise, N)
        x = x_ideal + x_noise
        y = y_ideal + y_noise

        linreg = LinearRegression(x, y)

        x_range = np.arange(xmin, xmax, xinc)
        y_model = linreg.model(x_range)
        y_pred = linreg.prediction_interval(x_range)
        y_conf = linreg.model_confidence_interval(x_range)

        color_model = 'b'
        color_conf = (0.5, 0.5, 1)
        color_pred = (1, 0.5, 0.5)
        color_sample = 'k'

        plt.fill_between(x_range, y_pred[1], y_pred[0], color=color_pred)
        plt.fill_between(x_range, y_conf[1], y_conf[0], color=color_conf)
        plt.plot(x, y, '.', color=color_sample)
        plt.plot(x_range, y_model, '-', color=color_model)
        plt.xlim(xmin, xmax)
        plt.show()


if __name__ == '__main__':
    unittest.main()
