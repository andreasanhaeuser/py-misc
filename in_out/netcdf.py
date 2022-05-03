#!/usr/bin/python
"""A collection of functions related to netcdf reading and writing.

    References
    ----------
    http://gfesuite.noaa.gov/developer/netCDFPythonInterface.html

    Authors
    -------
    Andreas Anhaeuser (AA) <anhaeus@meteo.uni-koeln.de>
    & Christopher Frank (CF) <cfrank@meteo.uni-koeln.de>
    Institute for Geophysics and Meteorology
    University of Cologne, Germany

    History
    -------
    2014-2016 (AA): Creation
    2016      (CF): Extending cast()-function, debugging.
"""

# standard modules
import os
from collections import Iterable
import warnings

# PyPI modules
import numpy as np
from netCDF4 import Dataset

# misc
from . import netcdf_helpers as helpers

################################################################
# main                                                         #
################################################################
def read_file(
        filename, varnames=None, ignore_varnames=None, str_method='regular',
        bounds=None,
        ):
    """Return content of netcdf file in two dicts.

        Parameters
        ----------
        filename : str
            path to netcdf file
        varnames : list of str, optional
            if given, only these variables are loaded
        ignore_varnames : list of str, optional
            these are ignored, even if explicitly listed in `varnames`
        str_method : str
            (e. g. 'regular', 'utf8', ...)

        Returns
        -------
        data : dict
        meta : dict

        History
        -------
        2022-01-18 (AA): expanduser() on filename
        2020-10-26 (AA): Add global dimension loader
        2019-09-20 (AA): Add per-variable dimension loader
        2018-01-03 (AA): Created
    """
    fn_long = os.path.expanduser(filename)
    if not os.path.isfile(fn_long):
        return None, None

    if ignore_varnames is None:
        ignore_varnames = []

    # initialize
    data = {}
    meta = {}

    varnames_load = varnames

    with Dataset(fn_long, 'r') as fid:
        # load dimensions
        dims = fid.dimensions
        meta['dimensions'] = [dim for dim in dims]
        meta['dimension_sizes'] = [dims[dim].size for dim in dims]
        
        cropper = helpers.get_cropper(fid, bounds)

        # create loading list if necessary
        if varnames_load is None:
            varnames_load = fid.variables.keys()

        # load variables
        # =================================================
        for varname in varnames_load:
            # ignore
            if varname in ignore_varnames:
                continue

            # load
            meta[varname] = read_var_atts(fid, varname)
            data[varname] = read_var_value(
                    fid, varname, method=str_method, cropper=cropper,
                    )
        # =================================================

    return data, meta

def write_file(*args, **kwargs):
    """Alias."""
    return write_file_v2(*args, **kwargs)

def write_file_v2(data, meta, filename):
    """Write netcdf file."""
    dimnames = meta['dimensions']

    with Dataset(filename, 'w') as fid:
        fid.set_auto_mask(False)

        # dimensions
        for dimname in dimnames:
            N = len(data[dimname])
            fid.createDimension(dimname, N)

        # variables
        for varname in sorted(data):
            if varname == 'char_pos':
                continue

            dimnames = meta[varname]['dimensions']
            dtype = meta[varname]['dtype']
            values = data[varname]
            if dtype == 'c':
                values = cast_string(values)
            if np.issubdtype(dtype, np.dtype('S')):
                values = cast_string(values)
                dtype = 'c'

            vid = fid.createVariable(varname, dtype, dimnames)
            vid[:] = values

            # attributes
            attributes = meta[varname]
            ignore_keys = ('dimensions', 'dtype')
            for attname in attributes:
                if attname in ignore_keys:
                    continue
                value = attributes[attname]
                vid.setncattr(attname, value)

    return filename

def select_variables(file_in, file_out, varnames):
    """Create a copy with only selected variables.

        Parameters
        ----------
        file_in : str or netCDF4.Dataset
            name of input file or an opened netCDF4.Dataset of such
        file_out : str or netCDF4.Dataset
            name of ouput file or an opened netCDF4.Dataset of such
        varname : list of str
            variables to be copied

        Returns
        -------
        None
    """
    # Open file_in
    # ------------------------------------------------------
    if not isinstance(file_in, Dataset):
        # check type
        if not isinstance(file_in, str):
            message = (
                'file_in must be a filename or netCDF4.Dataset'
                + ', got %s' % file_in
                )
            raise TypeError(message)

        # check if file exists
        if not os.path.isfile(file_in):
            raise OSError('Input file does not exist: %s' % file_in)

        # call function recursively
        with Dataset(file_in, 'r') as fi:
            return select_variables(fi, file_out, varnames)

    # If this line is reached, file_in is a Dataset
    fi = file_in
    # ------------------------------------------------------

    # Open file_out
    # ------------------------------------------------------
    if not isinstance(file_out, Dataset):
        # check type
        if not isinstance(file_out, str):
            message = (
                'file_out must be a filename or netCDF4.Dataset'
                + ', got %s' % file_out
                )
            raise TypeError(message)

        # build directory
        dirname = os.path.dirname(file_out)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)

        # call function recursively
        data_model = fi.data_model
        with Dataset(file_out, 'w', format=data_model) as fo:
            return select_variables(fi, fo, varnames)

    # If this line is reached, file_out is a Dataset
    fo = file_out
    # ------------------------------------------------------

    # global attributes
    # ------------------------------------------------------
    attributes = fi.ncattrs()
    for attribute in attributes:
        value = fi.getncattr(attribute)
        fo.setncattr(attribute, value)
    # ------------------------------------------------------

    for varname in varnames:
        if varname not in fi.variables:
            continue

        vi = fi.variables[varname]
        dimensions = vi.dimensions
        dtype= vi.dtype
        attributes = vi.ncattrs()

        # create dimensions
        # ------------------------------------------------------
        for dimension in dimensions:
            if dimension in fo.dimensions:
                continue

            size = fi.dimensions[dimension].size
            fo.createDimension(dimension, size)
        # ------------------------------------------------------

        # create variable
        # ------------------------------------------------------
        vo = fo.createVariable(varname, dtype, dimensions)
        # ------------------------------------------------------

        # create attributes
        # ------------------------------------------------------
        for attribute in attributes:
            value = vi.getncattr(attribute)
            vo.setncattr(attribute, value)
        # ------------------------------------------------------

        # copy data
        # ------------------------------------------------------
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            vo[:] = vi[:]
        # ------------------------------------------------------

    return None


################################################################
# helpers - read                                               #
################################################################
def read_var_value(fid, varname, method='regular', cropper=None):
    """Return one nc variable, strings already processed."""
    # variable in file?
    if varname not in fid.variables.keys():
        filename = fid.name
        raise IOError('Variable "%s" not in %s' % (varname, filename))
    vid = fid.variables[varname]

    # string
    dtypes_str = ('c', 'S1')
    if vid.dtype in dtypes_str:
        return read_string_variable(fid, varname, method)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        if cropper is None:
            values = vid[:]
        else:
            idx = helpers.get_crop_indices(vid, cropper)

            values = vid[idx]
    return values

def read_var_atts(fid, varname):
    """Return variable attributes as dict."""
    vid = fid[varname]
    attnames = vid.ncattrs()
    atts = {}
    atts['dimensions'] = vid.dimensions
    atts['dtype'] = vid.dtype
    for attname in attnames:
        atts[attname] = vid.getncattr(attname)
    return atts

def read_string_variable(fid, varname, method='regular'):
    """Convert netcdf char-array into a python string."""
    raw = fid.variables[varname][:]
    if method == 'regular':
        return [x.tostring().decode().strip('\0') for x in raw]
    elif method == 'utf8':
        return [''.join([b.decode('UTF-8') for b in line]) for line in raw]
    else:
        raise ValueError('Unknown string type: %s' % method)

################################################################
# helpers - write                                              #
################################################################
def cast_string(data, Nchar=64, fill_char='\0'):
    """Cast to a predefined length."""
    assert isinstance(data, Iterable)

    if not isinstance(data, str):
        return [cast_string(word) for word in data]

    return data[:Nchar].ljust(Nchar, fill_char)

def cast(data, nctype, Nchar=64):
    """Cast value(s) into a format suitable for netcdf writing."""
    if nctype == 'c':
        return cast_string(data, Nchar)

    # scalar handling:
    if np.isscalar(data):
        x = cast(np.array([data]), nctype)
        return x[0]

    if nctype == 'f':
        return np.array(data, dtype=np.float32)

    ###################################################
    # INT TYPES                                       #
    ###################################################
    inttypes = {
            'i8' : [np.int64, -2**64, 2**64 - 1],
            'i'  : [np.int32, -2**31, 2**31 - 1],
            'h'  : [np.uint16,     0, 2**16 - 1],
            's'  : [np.int16, -2**15, 2**15 - 1],
            }

    found = False
    for inttype in inttypes:
        if nctype == inttype:
            dtype = inttypes[inttype][0]
            lo = inttypes[inttype][1]
            hi = inttypes[inttype][2]
            if lo < 0:
                filler = lo
            else:
                filler = hi
            x = np.round(data)
            x[np.isnan(x)] = filler
            assert np.min(x) >= lo,   'Value(s) too small for' + inttype
            assert np.max(x) <= hi,   'Value(s) too large for' + inttype
            found = True
    assert found, 'could not find nctype ' + nctype
    return np.array(x, dtype=dtype)
