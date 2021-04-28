# PyPI
import numpy as np

def get_regular_2d_grid(x0, y0, size_x, size_y, dx, dy):
    """Return x- and y-coordinates pair of 1d-arrays of odd length.

        (x0, y0) is the center point of the output array.

        Parameters
        ----------
        x0 : float
            x-coordinate of grid center
        y0 : float
            y-coordinate of grid center
        size_x : float
            total grid size in x-direction
        size_y : float
            total grid size in y-direction
        dx : float
            size of a grid cell in x-direction
        dy : float
            size of a grid cell in y-direction
            

        Returns
        -------
        x, y : 1d-arrays of odd length
            x- and y-coordindates, with (x0, y0) at the center
    """
    x = get_regular_sequence(x0, size_x, dx)
    y = get_regular_sequence(y0, size_y, dy)
    return x, y


################################################################
# helpers                                                      #
################################################################
def get_regular_sequence(center, size, delta):
    """Center will be a grid point."""
    assert size > 0
    assert delta > 0

    N = size // delta + 1

    # make N an odd number
    if N % 2 == 0:
        N -= 1

    N_half = N // 2

    low = center - delta * N_half
    high = low + delta * (N - 0.5)
    return np.arange(low, high, delta)
