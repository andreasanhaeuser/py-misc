"""Bar plot functions."""

# PyPI
import numpy as np
import matplotlib.pyplot as plt

def grouped_bars(x, height, width=0.8, bottom=0., color=None, labels=None, **kwargs):
    """Plots N groups of M bars.

        Parameters
        ----------
        x : array of shape (N,)
            bar positions.
        height : array of shape (M, N)
            bar heights
        width : float
            width of the whole group

        Returns
        -------
        list of plt.BarContainer
    """

    N = len(x)
    shape = np.shape(height)
    assert len(shape) == 2
    if shape[1] != N and shape[0] == N:
        height = height.T
        shape = np.shape(height)
    M = shape[0]
    assert np.shape(height) == (M, N)

    barwidth = width / M

    bottom = cast_2d(bottom, M, N)
    color = cast_1d(color, M)
    labels = cast_1d(labels, M)

    for m in range(M):
        xoffset = barwidth * (m - (M-1)/2)

        x_m = x + xoffset
        height_m = height[m]
        bottom_m = bottom[m]
        color_m = color[m]
        label_m = labels[m]

        plt.bar(
                x_m, height_m, barwidth, bottom_m, color=color_m, label=label_m,
                **kwargs
                )

################################################################
# helpers                                                      #
################################################################
def cast_1d(x, M):
    if x is None:
        return [None for m in range(M)]

    if len(x) == M:
        return x

    if np.shape(x)[:1] != (1,):
        x = x * np.ones(M)

    assert len(x) == M
    return x

def cast_2d(x, M, N):
    """Make sure output is of shape (M, N)."""
    if np.isscalar(x):
        x = x * np.ones(M)

    is_flat = len(np.shape(x)) == 1

    if is_flat and len(x) == M:
        x = np.outer(x, np.ones(N))
    elif is_flat and len(x) == N:
        x = np.outer(np.ones(M), x)

    assert np.shape(x) == (M, N)
    return x
