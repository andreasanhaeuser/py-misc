"""NetCDF helpers."""

# misc
from misc.geodata.bounds import convert_bounds_to_tuple

def get_cropper(fid, bounds):
    """Return a dict that contains cropping indices."""
    if bounds is None:
        return None

    dimensions = ('lat', 'lon')
    bounds = convert_bounds_to_tuple(bounds)
    cropper = {}

    for dimension in dimensions:
        limits = get_limits(bounds, dimension)
        name = get_dimension_name(fid, dimension)
        idx = get_dimension_idx(fid, name, limits)

        cropper[dimension] = {}
        cropper[dimension]['name'] = name
        cropper[dimension]['idx'] = idx

    return cropper

def get_crop_indices(vid, cropper):
    """Return a tuple of slices."""
    if cropper is None:
        return slice(None)

    dimnames = vid.dimensions
    idx = [slice(None) for dimname in dimnames]

    for dimension in cropper:
        cropper1d = cropper[dimension]
        dimname = cropper1d['name']
        idx1d = cropper1d['idx']
        if dimname is None:
            continue

        if dimname not in dimnames:
            continue

        if idx1d is None:
            continue

        axis = dimnames.index(dimname)
        idx[axis] = idx1d

    return tuple(idx)

################################################################
# helpers                                                      #
################################################################
def get_limits(bounds, dimension):
    """Return lower and upper bounds as pair."""
    if bounds is None:
        return None

    if dimension == 'lon':
        values = bounds[0], bounds[2]

    elif dimension == 'lat':
        values = bounds[1], bounds[3]

    else:
        raise NotImplementedError(dimension)

    low = min(values)
    high = max(values)
    return low, high

def get_dimension_name(fid, dimension):
    """Return the dimension name as it appears in the nc file.

        Parameters
        ----------
        fid : open nc file handle
        dimension : str

        Returns
        -------
        str or None
            None if no match found
    """
    candidates = get_name_candidates(dimension)
    for dimname in fid.dimensions:
        # check whether it is a dimname
        if dimname.lower() not in candidates:
            continue

        # check whether file contains coordinates to the dimension
        if dimname not in fid.variables:
            continue

        # it's a hit
        return dimname

    # fall-back value
    return None

def get_name_candidates(dimension):
    """Return dimension name candidates as tuple of str."""
    if dimension == 'lon':
        return ('lon', 'long', 'longitude', 'xlon')
    if dimension == 'lat':
        return ('lat', 'latitude', 'xlat')
    raise ValueError(dimension)

def get_dimension_idx(fid, name, limits):
    """Return a bool array or None."""
    if name is None:
        return None

    if limits is None:
        return None

    values = fid.variables[name][:]
    assert len(values.shape) == 1
    low, high = limits

    idx_low = values >= low
    idx_high = values <= high
    idx = idx_low & idx_high

    N = len(values)
    ns = [n for n in range(N) if idx[n]]
    nlo = min(ns)
    nhi = max(ns) + 1
    return slice(nlo, nhi)
