# local modules
from misc import string_utils as aa_str

################################################################
# legacy                                                       #
################################################################
def write_file_v1(
        data, meta, filename, vartable,
        ignore_key_errors=False, add_var_meta=True, add_glob_meta=True,
        string_length=None,
        ):
    """Create a netcdf file.

    Parameters
    ----------
    data : dict
        contains the data to be written (see Usage section)
    meta : dict
        additional information on data (see Usage section)
    filename : str
        name of the output netcdf file
    vartable : str
        path to a vartable (see Usage section)
    ignore_key_errors : bool, optional
        if True, then keys that are listed in the vartable but not in data
        are ignored, otherwise a KeyError is raised, default: False
    add_var_meta : bool, optional
        if True, then all information in meta is written as attributes to the
        corresponding variables. Keys that are already listed in the variable
        table are ignored. Default: True
    add_glob_meta : bool, optional
        if True, then the entry 'glob' of meta (if present) is written as
        global attributes.
    string_length : int, optional
        if this is set, variables that have nc-type 'c' are automatically
        adjusted to this length (either filled by '\0'-characters of cut).


    Returns
    -------
    None


    Usage
    -----
    data and meta must contain all keys that appear in the py_name column of
    the vartable.

    The vartable file must be of the following shape:

    py_name  | nc_name     | is_dim | is_var | dims    |  nc_type |  nc_units |
    ---------+-------------+--------+--------+---------+----------+-----------+
    level    | z_level     |      2 |        |         |          |           |
    pres     | pressure    |        |      1 |       2 |        f |        Pa |
    time     | time        |      1 |      1 |       1 |        f |         s |
    lat      | latitude    |        |      1 |         |        f |       deg |
    lon      | longitude   |        |      1 |         |        f |       deg |
    alt      | hsurf       |        |      1 |         |        i |         m |
    temp     | temperature |        |      1 |     2,1 |        f |         K |

    Syntax:
     * columns / separators :
        columns must be separated by exactly one '|' character. There is no
        separator before the first column, but one after the last column. This
        is essential for the function to work properly.
    * headers :
        The first line contains the column headers.
    * ignored lines / comments / decorators:
        - Lines that have not exactly the same numbers of columns as the header
          line are ignored.
        - Lines starting with '#' are also ignored.
        - The second line (horizontal rule) is just a decorator for the human
          reader but not neccessary. The function has the same behaviour
          without this line.

    Compulsory headers:
    * py_name : name of the variable as it appears in data
    * nc_name : name of the variable as it shall be written to the nc-file
    * is_dim :  leave blank, if that variable is not a nc-dimension, otherwise
                label each dimension variable with a unique string (or number)
                that does not contain the '|' character
    * is_var :  leave blank or '0' if it is only a dimension but not a
                variable, otherwise use any other character. '1' is
                recommended.
    * dims :    dimensions labels of the variable as they appear in the is_dim
                column of the appropriate dimension. Separate multiple
                dimensions by ','. It is important that the dimensions appear
                in the right order, i. e. that they are consistent with the
                dimensions of data[py_name]. For dimensionless variables, leave
                this column blank.

    Recommended headers:
    * nc_type : The netcdf type of the variable, such as 'f', 'i', 'c', ...
                Default: 'f'
    * gain, offset :
                the python variable is converted to the netcdf variable by
                     nc_value = (py_value - offset) / gain
                This is constistent with read_file(), i. e. if you read and
                then write a dataset with the same variable table, nothing will
                have changed. Default: gain = 1, offset = 0.
    * nc_units :
                Units for the nc_variable. Note: The variable attribute name
                will be 'units', not 'nc_units'.
    * standard_name :
                CF standard name as in
                http://cfconventions.org/standard-names.html

    Optional headers:
      Any other headers will be added as variable attributes.


    Note
    ----
    The nc-type 'c' (character/string) has not been tested and will probably
    cause problems in most cases.


    Bugs
    ----
    Has not been tested extensively. Probably contains bugs. Better check your
    output for correctness. Please report bugs!


    Author
    ------
    Andreas Anhaeuser (AA) <anhaeus@meteo.uni-koeln.de>
    Institute for Geophysics and Meteorology
    University of Cologne, Germany

    History
    -------
    2015 (AA): Created
    2016-09-13 (AA): `string_length` extention
    """
    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    assert string_length is None or isinstance(string_length, int)
    assert isinstance(data, dict)
    assert isinstance(meta, dict)
    assert isinstance(filename, str)
    assert isinstance(vartable, str)
    assert os.path.isfile(vartable)

    ###################################
    # OPEN NC-FILE                    #
    ###################################
    if os.path.isfile(filename):
        os.remove(filename)
    i = filename.rfind('/')
    path = filename[:i]
    if not os.path.isdir(path):
        os.makedirs(path)
    nc = ncfile(filename, 'a')

    ###################################
    # READ VARTABLE                   #
    ###################################
    init = False
    vt = open(vartable, 'r')
    lines = vt.readlines()
    vt.close()
    Nlines = len(lines)

    for mode in ['dims', 'vars']:
        for nline in range(Nlines):
            line = lines[nline]
            # the last word is removed because it is just the new line char:
            words = line.split('|')[:-1]

            # delete leading and trailing white spaces:
            for i in range(len(words)):
                word = words[i]
                while len(word) > 0 and word[0] == ' ':
                    word = word[1:]
                while len(word) > 0 and word[-1] == ' ':
                    word = word[:-1]
                words[i] = word

            ###############################
            # HEADER LINE                 #
            ###############################
            if nline == 0:
                if init:
                    continue
                keys = []
                keypos = []
                pypos = None      # py_name position
                ncpos = None      # nc_name position
                gainpos = None    # gain position
                offsetpos = None  # offset position
                typos = None      # nc_type position
                isdimpos = None   # is_dim position
                isvarpos = None   # is_var position
                dimpos = None     # dims position
                otherpos = []     # position of remaining headers

                i = 0
                for word in words:
                    if word == 'py_name':
                        pypos = i
                    elif word == 'nc_name':
                        ncpos = i
                    elif word == 'gain':
                        gainpos = i
                    elif word == 'offset':
                        offsetpos = i
                    elif word == 'nc_type':
                        typos = i
                    elif word == 'is_dim':
                        isdimpos = i
                    elif word == 'is_var':
                        isvarpos = i
                    elif word == 'dims':
                        dimpos = i
                    else:
                        otherpos.append(i)
                    # change header 'nc_units' to key 'units':
                    if word == 'nc_units':
                        word = 'units'
                    keys.append(word)
                    i += 1

                # initialize:
                dimnames = {}
                Ncols = len(words)
                init = True

                # next line:
                continue

            ###############################
            # SKIP SEPARATOR LINES        #
            ###############################
            if len(words) != Ncols:
                continue

            ###############################
            # SKIP COMMENT LINES          #
            ###############################
            if any(words[0]) and words[0][0] == '#':
                continue

            ###############################
            # CREATE DIMENSIONS           #
            ###############################
            if mode == 'dims':
                is_dim = words[isdimpos]
                if not is_dim:
                    continue
                nckey = words[ncpos]
                pykey = words[pypos]
                if np.isscalar(data[pykey]):
                    L = 1
                else:
                    L = len(data[pykey])
                nc.createDimension(nckey, L)
                dimnames[is_dim] = nckey
                continue

            ###########################
            # VARIABLE SETTINGS       #
            ###########################
            # variable names
            nckey      = words[ncpos]
            pykey      = words[pypos]

            # dimnumbers
            dimnumbers = words[dimpos].split(',')
            if dimnumbers == ['']:
                dimnumbers = []

            # isvar
            isvar = words[isvarpos]
            if isvar.lower() in ['', 'n', 'no', '0']:
                continue

            # nctype
            if typos is None:
                nctype = 'f'
            else:
                nctype = words[typos]
            if nctype == '':
                nctype = 'f'

            # gain
            if gainpos is None:
                gain = 1
            elif words[gainpos] == '':
                gain = 1
            else:
                gain = float(words[gainpos])

            # offset
            if offsetpos is None:
                offset = 0
            elif words[offsetpos] == '':
                offset = 0
            else:
                offset = float(words[offsetpos])

            # attributes
            atts = []

            ##############################
            # CREATE VARIABLE            #
            ##############################
            try:
                if nctype != 'c':
                    raw = (data[pykey] - offset) / gain
                else:
                    raw = data[pykey]
            except KeyError as e:
                if ignore_key_errors:
                    continue
                else:
                    raise e
            dims = tuple([dimnames[dn] for dn in dimnumbers])
            varid = nc.createVariable(nckey, nctype, dims)

            ##############################
            # ADJUST STRING LENGTH       #
            ##############################
            if nctype == 'c' and string_length:
                raw = aa_str.adjust_string_length(raw, string_length)

            ##############################
            # WRITE DATA TO VARIABLE     #
            ##############################
            value = cast(raw, nctype)

            try:
                varid.assignValue(value)
            except Exception as e:
                message = ''.join([
                    'pykey: %s' % pykey,
                    '\nnckey: %s' % nckey,
                    '\nshape: %s' % str(np.shape(value)),
                    ])
                print(message)
                raise e

            ###############################
            # WRITE VARIABLE ATTRIBUTES   #
            ###############################
            var = nc.variables[nckey]
            for i in range(Ncols):
                if i not in otherpos:
                    continue
                if keys[i] in ['py_units']:
                    continue
                setattr(varid, keys[i], words[i])

            ###############################
            # ADD VARIABLE META DATA      #
            ###############################
            if add_var_meta:
                if pykey not in meta.keys():
                    continue
                metakeys = meta[pykey].keys()
                for att in metakeys:
                    if att not in keys:
                        setattr(varid, att, meta[pykey][att])
        vt.close()

    ###################################
    # ADD GLOBAL META                 #
    ###################################
    if add_glob_meta and 'glob' in meta.keys():
        gatts = sorted(meta['glob'].keys())
        for gatt in gatts:
            setattr(nc, gatt, meta['glob'][gatt])

    nc.close()


#################################################
# HELPER FUNCTIONS                              #
#################################################

################################################################
# deprecated                                                   #
################################################################
def read_file_without_vartable_old(filename, read_atts=True):
    """*Return two dictionaries.

        Parameters
        ----------
        filename : str
            path to the nc-file
        read_atts : bool, optional
            if True, then all variable attributes are read. Default: True

        Returns
        -------
        (dict, dict)
            The first dictionary contains the variable values. The second
            dictionary consist itself of dictionaries which contains meta
            information.

        Note
        ----
        This is designed as a helper function to read_file()
    """
    # open file
    fid = ncfile(filename, 'r')
    data = {}
    meta = {}

    # get variable list
    vnames = fid.variables.keys()
    for vname in vnames:
        varid = fid.variables[vname]
        typechar = varid.typecode()
        meta[vname] = {}
        offset  = 0
        scale   = 1
        filler  = np.nan
        missing = np.nan
        val_min = np.nan
        val_max = np.nan

        # read variable attributes
        if read_atts:
            # get attribute list
            attnames = dir(varid)
            for attname in attnames:
                if attname == 'assignValue':
                    continue
                if attname == 'getValue':
                    continue
                if attname == 'typecode':
                    continue
                if attname == '_FillValue':
                    filler = getattr(varid, attname)[0]
                if attname == 'missing_value':
                    missing = getattr(varid, attname)[0]
                if attname == 'valid_min':
                    val_min = getattr(varid, attname)[0]
                if attname == 'valid_max':
                    val_max = getattr(varid, attname)[0]
                if attname == 'valid_range':
                    val_min, val_max = getattr(varid, attname)
                if attname == 'add_offset':
                    offset = getattr(varid, attname)[0]
                    continue
                if attname == 'scale_factor':
                    scale = getattr(varid, attname)[0]
                    continue
                meta[vname][attname] = getattr(varid, attname)

        # exclude invalid data:
        raw = fid.variables[vname].getValue()
        shape = np.shape(raw)

        if typechar in 'df' and read_atts:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                # CASE: array
                if np.shape(raw) != ():
                    if not np.isnan(filler):
                        raw[raw==filler] = np.nan
                    if not np.isnan(missing):
                        raw[raw==filler] = np.nan
                    if not np.isnan(val_min):
                        raw[raw < val_min] = np.nan
                    if not np.isnan(val_max):
                        raw[raw > val_max] = np.nan
                # CASE: scalar
                else:
                    if not np.isnan(filler):
                        if raw == filler:
                            raw = np.nan
                    if not np.isnan(missing):
                        if raw == filler:
                            raw = np.nan
                    if not np.isnan(val_min):
                        if raw < val_min:
                            raw = np.nan
                    if not np.isnan(val_max):
                        if raw > val_max:
                            raw = np.nan
            data[vname] = scale * raw + offset
        else:
            data[vname] = raw

    # get global attributes
    if read_atts:
        meta['glob'] = {}
        gattnames = dir(fid)
        for gattname in gattnames:
            meta['glob'][gattname] = getattr(fid, gattname)

    fid.close()
    return data, meta

def read_file_old(
        filename,
        vartable=None,
        read_atts=True,
        ignore_key_errors=False,
        read_dims=True,
        messages=True,
        merge_strings=True,
        ):
    """Return two dictionaries.

    Parameters
    ----------
    filename : str
        path to the nc-file
    vartable : str, optional
        path to the vartable file
    read_atts : bool, optional
        if True, then all variable attributes are read. Attributes that have
        same name as a header in the variable table are ignored. Default is
        True.
    ignore_key_errors : bool, optional
    messages : bool, optional
        Issue warning and info messages on screen or not.
    merge_strings : bool, optional
        If True, 2d-character matrices are converted to arrays of strings


    Returns
    -------
    (dict, dict)
        The first dictionary contains the variable values.
        The seconds dictionary consists itself of dictionaries which contains
        meta information.


    Note
    ----
    If no vartable is specified, then all variables are read.

    The vartable file must be of the following shape:

    py_name     | nc_name | gain | offset | other_optional_headers
    ------------+---------+------+--------+-----------------------
    varname1    | ncvn1   |    1 |      0 | further info
    varname2    | ncvn2   | 1000 |      0 | further info
    ...

    All headers except py_name and nc_name are optional.

    py_units :
    The header py_units is converted to the key 'units'.

    gain, offset :
    If gain and/or offset are given, they must be numbers. If there are no
    headers called 'gain' and/or 'offset', the values will be set to 1 and/or
    0, respectively for all variables. They are used for conversions:

        py_value = nc_value * gain + offset


    Example
    -------
    For example, if height is given in km above ground in the nc-file and you
    want it in meters above sea level in your python variable, your table will
    look something like this, assuming the ground to be 205m ASL:

    py_name  | nc_name    | gain | offset | long_name              | py_units
    ---------+------------+------+--------+------------------------+---------
    height   | elevation  | 1000 |    205 | height above sea level |        m
    varname2 | ncvn2      |    1 |      0 | long variable name     |    kg/mK

    General netcdf write information
    --------------------------------
    For general information see
    http://gfesuite.noaa.gov/developer/netCDFPythonInterface.html

    Author
    ------
    Written in 2015
    by Andreas Anhaeuser
    Insitute for Geophysics and Meteorology
    University of Cologne
    Germany
    <anhaues@meteo.uni-koeln.de>
    """
    if not os.path.isfile(filename):
        raise IOError('Cannot find file ' + filename)

    if vartable is None:
        return read_file_without_vartable(
                filename=filename, read_atts=read_atts)

    ###################################
    # READ VARTABLE                   #
    ###################################
    vt = open(vartable, 'r')
    nc = ncfile(filename, 'r')

    init = False
    for line in vt:
        # the last word is removed because it is just the new line character:
        words = line.split('|')[:-1]

        # delete leading and trailing white spaces:
        for i in range(len(words)):
            word = words[i]
            while len(word) > 0 and word[0] == ' ':
                word = word[1:]
            while len(word) > 0 and word[-1] == ' ':
                word = word[:-1]
            words[i] = word

        ###############################
        # HEADER LINE                 #
        ###############################
        if not init:
            keys   = []
            keypos = []
            pypos     = None  # py_name position
            ncpos     = None  # nc_name position
            isvarpos  = None  # nc_name position
            gainpos   = None  # gain position
            offsetpos = None  # offset position
            otherpos = []     # position of remaining headers

            i = 0
            for word in words:
                if word == 'py_name':
                    pypos = i
                elif word == 'nc_name':
                    ncpos = i
                elif word == 'is_var':
                    isvarpos = i
                elif word == 'gain':
                    gainpos = i
                elif word == 'offset':
                    offsetpos = i
                else:
                    otherpos.append(i)
                # change header 'py_units' to key 'units':
                if word == 'py_units':
                    word = 'units'
                keys.append(word)
                i += 1

            # initialize:
            data = {}
            meta = {}
            Ncols = len(words)
            init = True

            # next line:
            continue

        ###############################
        # SKIP SEPERATOR LINES        #
        ###############################
        if len(words) != Ncols:
            continue

        ###############################
        # SKIP COMMENT LINES          #
        ###############################
        if any(words[0]) and words[0][0] == '#':
            continue

        pykey = words[pypos]
        nckey = words[ncpos]

        ###############################
        # VARIABLE ?                  #
        ###############################
        if isvarpos is None:
            isvar = True
        elif words[isvarpos] in ['', '0']:
            isvar = False
        else:
            isvar = True

        if not isvar and not read_dims:
            continue

        ###############################
        # READ DIMENSION              #
        ###############################
        if not isvar and read_dims:
            try:
                data[pykey] = np.arange(nc.dimensions[nckey])
                meta[pykey] = {'is_dim': '1'}
                continue
            except KeyError as error:
                if messages:
                    print('WARNING: No dimension named "' + nckey +
                            '" in nc-file ' + filename)
                # next line:
                continue
            else:
                raise error

        ###############################
        # READ DATA                   #
        ###############################
        try:
            raw = nc.variables[nckey].getValue()
        except KeyError as error:
            if ignore_key_errors:
                if messages:
                    print('WARNING: No variable named "' + nckey +
                            '" in nc-file ' + filename)
                # next line:
                continue
            else:
                raise error


        #################################
        # GAIN AND OFFSET FROM FILE     #
        #################################
        meta[pykey] = {}
        varid = nc.variables[nckey]
        typechar = varid.typecode()
        offset  = 0
        scale   = 1
        filler  = np.nan
        missing = np.nan
        val_min = np.nan
        val_max = np.nan

        # read variable attributes
        # get attribute list
        attnames = dir(varid)
        for attname in attnames:
            if attname == 'assignValue':
                continue
            if attname == 'getValue':
                continue
            if attname == 'typecode':
                continue
            if attname == '_FillValue':
                filler = getattr(varid, attname)[0]
            if attname == 'missing_value':
                missing = getattr(varid, attname)[0]
            if attname == 'valid_min':
                val_min = getattr(varid, attname)[0]
            if attname == 'valid_max':
                val_max = getattr(varid, attname)[0]
            if attname == 'valid_range':
                val_min, val_max = getattr(varid, attname)
            if attname == 'add_offset':
                offset = getattr(varid, attname)[0]
                continue
            if attname == 'scale_factor':
                scale = getattr(varid, attname)[0]
                continue
            if read_atts and attname not in keys:
                meta[pykey][attname] = getattr(varid, attname)


        # exclude invalid data:
        shape = np.shape(raw)
        if typechar in 'df' and read_atts:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                if isinstance(raw, Iterable):
                    if not np.isnan(filler):
                        raw[raw==filler] = np.nan
                    if not np.isnan(missing):
                        raw[raw==filler] = np.nan
                    if not np.isnan(val_min):
                        raw[raw < val_min] = np.nan
                    if not np.isnan(val_max):
                        raw[raw > val_max] = np.nan
                else:
                    if not np.isnan(filler) and raw==filler:
                        raw = np.nan
                    if not np.isnan(missing) and raw==filler:
                        raw = np.nan
                    if not np.isnan(val_min) and raw < val_min:
                        raw = np.nan
                    if not np.isnan(val_max) and raw > val_max:
                        raw = np.nan

        if typechar not in 'c' and read_atts:
            scaled = scale * raw + offset
        else:
            scaled = raw

        #################################
        # GAIN AND OFFSET FROM VARTABLE #
        #################################
        # gain:
        if gainpos is None:
            gain = 1
        elif words[gainpos] == '':
            gain = 1
        else:
            gain = float(words[gainpos])
        # avoid cast from int to float if not necessary
        if int(gain) == gain:
            gain = int(gain)

        # offset:
        if offsetpos is None:
            offset = 0
        elif words[offsetpos] == '':
            offset = 0
        else:
            offset = float(words[offsetpos])
        # avoid cast from int to float if not necessary
        if int(offset) == offset:
            offset = int(offset)

        # scale:
        if raw.dtype != 'S1':
            data[pykey] = scaled * gain + offset
        else:
            data[pykey] = scaled

        ###################################################
        # MERGE STRINGS                                   #
        ###################################################
        if raw.dtype == 'S1' and len(np.shape(scaled)) == 2:
            data[pykey] = map(''.join, scaled)

        ###############################
        # ADD META DATA               #
        ###############################
        for i in range(Ncols):
            if i not in otherpos:
                continue
            meta[pykey][keys[i]] = words[i]

    ###############################
    # READ GLOBAL ATTRIBUTES      #
    ###############################
    if read_atts:
        if 'glob' not in meta.keys():
            meta['glob'] = {}
        gatts = dir(nc)
        ignore = ['close', 'createDimension', 'createVariable', 'flush','sync']
        for gatt in gatts:
            if gatt in ignore:
                continue
            meta['glob'][gatt] = getattr(nc, gatt)

    nc.close()
    vt.close()
    return data, meta
