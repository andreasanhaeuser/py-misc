import numpy as np
import matplotlib.pyplot as plt

class Location(object):
    """Abstract class for locations."""
    def __init__(self, name=None):
        self.name = name

    def plot(
            self, ax=None, lw=None, facecolor=None, edgecolor=None,
            zorder=None,
            ):
        pass

def PolygonLocation(Location):
    def __init__(self, polygon, name=None):
        self.polygon = polygon
        self.name = name

    def plot(
            self, ax=None, lw=None, facecolor=None, edgecolor=None,
            zorder=None,
            ):
        pass

def CircleLocation(Location):
    def __init__(self, x, y, radius, name=None, Nsegments=1024):
        self.name = name

        phi_inc = 2 * np.pi / Nsamples
        phi_min = 0.
        phi_max = 2 * np.pi
        phis = np.arange(phi_min, phi_max, phi_inc)

        x_edge = x + radius * np.cos(phis)
        y_edge = y + radius * np.sin(phis)

        ax_orig = plt.gca()
        if ax is not None:
            plt.sca(ax)
