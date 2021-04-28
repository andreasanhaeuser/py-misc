# PyPI
import numpy as np
import pyproj


def compute_coordinates_from_transform(meta):
    """Return x and y as two 1d-arrays."""
    Nx = meta['width']
    Ny = meta['height']
    transform = meta['transform']

    # i : indices in x-direction
    # j : inideces in y-direction
    i = np.arange(Nx)
    j = np.arange(Ny)

    # apply affine transformation from index to coordinate
    # `x0` and `y0` remain unused
    # (the zeros can be replace by any value without altering the output)
    x, y0 = transform * (i, 0)
    x0, y = transform * (0, j)
    return x, y

def get_coordinate_order():
    """Return 'lonlat' or 'latlon'."""
    # setup
    # =====================================================
    module = pyproj

    # starting at this version, use the new order
    main_version_change = 2
    sub_version_change = 6      # not sure about this value;
                                # possibly some sub-versions earlier
    old_order = 'lonlat'
    new_order = 'latlon'
    # =====================================================

    # determine version
    # =====================================================
    count = 0
    pivot = 0
    version_str = module.__version__
    version = []
    while '.' in version_str[pivot:]:
        pivot_end = version_str.index('.', pivot)
        value_str = version_str[pivot:pivot_end]

        try:
            value = int(value_str)
        except ValueError:
            value = value_str
        version.append(value)

        pivot = pivot_end + 1

    main_version = version[0]
    sub_version = version[1]
    # =====================================================


    # determine coordinate order
    # =====================================================

    # check main version
    # ------------------------------------------------
    if main_version < main_version_change:
        return old_order

    if main_version > main_version_change:
        return new_order

    # sanity check
    assert main_version == main_version_change
    # ------------------------------------------------

    # check sub-version
    # ------------------------------------------------
    if sub_version < sub_version_change:
        return old_order

    return new_order
