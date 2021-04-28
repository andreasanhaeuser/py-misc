# PyPI
import numpy as np

# misc
from . import tif

def read_tiled_tif(filenames, bounds, aggregation='last'):
    """Load a tif that is stored in several tiles

        Parameters
        ----------
        filenames : iterable of str
        bounds : tuple or GeoDataFrame
        aggregation : {'last', 'mean', 'max', 'min', 'sum'}

        Returns
        -------
        data : array
        meta : dict
    """
    if aggregation != 'last':
        raise NotImplementedError(aggregation)

    init = False
    for filename in filenames:
        # load
        try:
            fdata, fmeta = tif.load_file_raw(filename, bounds)
        except ValueError as e:
            # This happens if bounds does not overlap tile
            continue

        # initialize
        # ==============================
        if not init:
            data = fdata
            meta = fmeta
            init = True
            continue
        # ==============================

        data, meta = extend_coordinates(data, meta, fdata, fmeta)
        data, meta = include_tile(data, meta, fdata, fmeta)

    assert init

    return data, meta

################################################################
# helpers - combine                                            #
################################################################
def include_tile(data, meta, fdata, fmeta):
    nxs = get_tile_index(meta, fmeta, 'x_native')
    nys = get_tile_index(meta, fmeta, 'y_native')
    for fny, ny in enumerate(nys):
        data[ny, nxs] = fdata[fny]

    return data, meta

def get_tile_index(meta, fmeta, key):
    coords = meta[key]
    fcoords = fmeta[key]
    N = len(coords)
    idx = np.zeros_like(coords, dtype=bool)
    for n, coord in enumerate(coords):
        if coord in fcoords:
            idx[n] = True

    ns = np.arange(N, dtype=int)[idx]
    return ns

def extend_coordinates(data, meta, fdata, fmeta):
    f = extend_coordinates_one_axis
    data, meta = f(data, meta, fdata, fmeta, 'x_native', 0)
    data, meta = f(data, meta, fdata, fmeta, 'y_native', 1)
    return data, meta

def extend_coordinates_one_axis(data, meta, fdata, fmeta, key, axis):
    data = data.T
    fdata = fdata.T

    coords = meta[key]
    fcoords = fmeta[key]

    if is_sorted_ascending(coords):
        assert is_sorted_ascending(fcoords)
        sorting = 'a'
    elif is_sorted_descending(coords):
        assert is_sorted_descending(fcoords)
        sorting = 'd'
    else:
        print(coords)
        raise Exception()

    # create empty data slice to be inserted
    shape = np.shape(data)
    shape_empty = shape[:axis] + (1,) + shape[axis+1:]
    empty = np.nan * np.zeros(shape_empty)

    for coord in fcoords:
        if coord in coords:
            continue

        if sorting == 'a':
            if coord < coords[0]:
                i = 0
            else:
                i = np.where(coord > coords)[0][-1] + 1

        if sorting == 'd':
            if coord > coords[0]:
                i = 0
            else:
                i = np.where(coord < coords)[0][-1] + 1

        coords = np.insert(coords, i, coord)
        data = np.insert(data, i, empty, axis=axis)
        assert data.shape[axis] == len(coords)

    meta[key] = coords
    data = data.T

    if sorting == 'a':
        assert is_sorted_ascending(coords)
    else:
        assert is_sorted_descending(coords)
    return data, meta

def is_sorted_ascending(a, strict=False):
    N = len(a)
    for n in range(N-1):
        if a[n] < a[n+1]:
            continue
        if a[n] == a[n+1] and not strict:
            continue

        return False

    return True

def is_sorted_descending(a, strict=False):
    N = len(a)
    for n in range(N-1):
        if a[n] > a[n+1]:
            continue
        if a[n] == a[n+1] and not strict:
            continue

        return False

    return True
