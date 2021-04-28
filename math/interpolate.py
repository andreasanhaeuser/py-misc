#!/usr/bin/python3

# standard
from copy import deepcopy as copy
import collections
import datetime as dt
import itertools
import numbers

# PyPI
import numpy as np

# misc
import misc.datetime_utils as du
import misc.string_utils as str_utils
from misc.math.nearest_neighbour import nearest_neighbour
from misc.chronometer import Chronometer
    

#####################################
# INTERPOLATION                     #
#####################################
def interpolate_1d_linear(
        x_in,
        val_in,
        x_out,
        x_tolerance=0,
        out_of_bounds='nearest',
        assume_sorted=True,
        many_nans=False,
        ):
    """Linear 1D-interpolation for numbers or datetime.

        x_in and x_out can be floats or datetime.datetime instances.
       
        Parameters
        ----------
        x_in : array of floats or list of datetime.datetime
        x_out : array of floats or list of datetime.datetime
        val_in : array of floats
            same length as x_in
        x_tolerance : float or datetime.timedelta, optional
            only nearest neighbours that are close than x_tolerance are
            considered. 0 causes infinite tolerance! Default: 0.
        out_of_bounds : {'zero', 'nan', 'nearest', 'ext', float}
            out of bounds value are replaced by this
            'ext' : extrapolate
        assume_sorted : bool
            if True, the function assumes that x_in be sorted in rising order
        many_nans : bool
            no effect (kept for backwards compatibility)
        
        Returns
        -------
        val_out : array
            The output is a linear interpolation between the sample values
            given in values_sample at the times given in val_in.     

        Notes
        -----
        float vs datetime.datetime:
            If some of x_in, x_out, or x_tolerance are given in numerical
            values and others as instances of datetime.datetime (or
            datetime.timedelta in the case of x_tolerance), the numericals are
            considered as seconds (since 1970).

        x_tolercance:
            if 0, this is interpreted as infinite tolerance!
       
        Author
        ------
        Written in 2014-2016
        by Andreas Anhaeuser
        Insitute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """
    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    assert len(x_in) == len(val_in)

    ###################################################
    # CONVERSIONS                                     #
    ###################################################
    if isinstance(x_in[0], dt.datetime):
        x_in = du.datetime_to_seconds(x_in)
    if isinstance(x_out[0], dt.datetime):
        x_out = du.datetime_to_seconds(x_out)
    if isinstance(x_tolerance, dt.timedelta):
        x_tolerance = x_tolerance.total_seconds()

    xis = np.array(x_in)
    xos = np.array(x_out)
    vis = np.array(val_in)
    tol = x_tolerance

    # tol: 0 --> inf
    if tol == 0:
        tol = np.inf

    ###################################################
    # SORT                                            #
    ###################################################
    # (recursively call function if unsorted)
    if not assume_sorted:
        idx = np.argsort(xis)
        xis = xis[idx]
        vis  = vis[idx]
        
        return interpolate_1d_linear(xis=xis,
            vis=vis, xos=xos, x_tolerance=tol, 
            out_of_bounds=out_of_bounds, assume_sorted=True)
        
    ###################################################
    # OUT OF BOUNDS                                   #
    ###################################################
    if out_of_bounds in ['nan', 'NaN']:
        oob  = 'val'
        oobv = np.nan
    elif out_of_bounds in ['zero', 'zeros', '0']:
        oob  = 'val'
        oobv = 0.
    elif isinstance(out_of_bounds, numbers.Number):
        oob  = 'val'
        oobv = out_of_bounds
    elif out_of_bounds == 'nearest':
        oob  = 'nn'
        oobv = None
    elif out_of_bounds == 'ext':
        raise NotImplementedError(
                'Current implementation seems buggy (AA 2016-09-25.)')
        oob  = 'ex'
        oobv = None
    else:
        raise ValueError('Uncrecognized out_of_bounds value: ' + 
                str(out_of_bounds))
    
    ###################################################
    # SPECIAL CASES                                   #
    ###################################################
    # special case: empty input array:
    if len(x_in) == 0:
        return np.array([oobv] * len(x_out))
    
    # special case: empty output array:
    if len(xos) == 0:
        return np.array([])
    
    ###################################################
    # MANY                                            #
    ###################################################
    neg = np.isnan(xis) | np.isnan(vis)
    use = ~neg
    xis = xis[use]
    vis = vis[use]

    ###################################################
    # SPECIAL CASES                                   #
    ###################################################
    # special case: input array length 0
    if len(xis) == 0:
        return np.nan * xos

    # special case: input array length 1
    if len(xis) == 1:
        if oobv is not None:
            # idea: all entries which are close enough (indices `pos`), are
            # assigned vis[0], all others are oobv

            # default
            vos = oobv * np.ones(len(xos))

            # find positions close enough
            diff = np.abs(xos - xis)
            pos = diff <= tol

            # assign them the only existing value
            vos[pos] = vis[0]
        else:
            # all out-values are equal to the single in-value
            vos = vis[0] * np.ones(len(xos))
        return vos
            
    ###################################################
    # INTERPOLATION FUNCTION                          #
    ###################################################
    def f(x1, x2, y1, y2, x):
        """Linear interpolation function."""
        # special case
        if x1 == x2:
            return y1

        # regular case
        m = (y2 - y1) / (x2 - x1)
        return y1 + m * (x - x1)        

    ###################################################
    # INITIALIZATION                                  #
    ###################################################
    N = len(xos)
    I = len(xis)
    vos = np.nan * np.empty(N)
    xlos = np.nan * np.empty(N)
    xhis = np.nan * np.empty(N)

    chrono = Chronometer(N, header='Interpolation')
    for n in range(N):
        chrono.loop_and_show()
        ###################################################
        # FIND NEAREST NEIGHBOUR                          #
        ###################################################
        xo = xos[n]
        xlo, ilo, vlo = nearest_neighbour(
            xis, xo,
            direction='down', values=vis, skip_nans=False, assume_sorted=True)
            
        xhi, ihi, vhi = nearest_neighbour(
                xis, xo, direction='up', values=vis,
                skip_nans=True, assume_sorted=True)

        # out of bounds:
        if xlo is None:
            xlo = -np.inf
        if xhi is None:
            xhi = +np.inf

        # out of tolerance:
        if xo - xlo > tol:
            xlo = -np.inf
        if xhi - xo > tol:
            xhi = +np.inf
        
        # distances:
        dlo = xo - xlo
        dhi = xhi - xo

        ###################################################
        # REGULAR CASE                                    #
        ###################################################
        if np.isfinite(dlo) and np.isfinite(dhi):
            vos[n] = f(xlo, xhi, vlo, vhi, xo)
            continue
      
        ###################################################
        # OUT OF TOLERANCE ON BOTH SIDES                  #
        ###################################################
        if ~np.isfinite(dlo) and ~np.isfinite(dhi):
            if oob == 'val':
                vos[n] = oobv
            continue

        ###################################################
        # OUT OF BOUNDS ON LOWER SIDE                     #
        ###################################################
        if ~np.isfinite(dlo) and np.isfinite(dhi):
            # constant value
            if oob == 'val':
                vos[n] = oobv
                continue

            # nearest neighbour
            elif oob == 'nn':
                vos[n] = vhi
                continue

            # extrapolate
            elif oob  == 'ex':
                # upper bound
                if ihi == I - 1:
                    vos[n] = vhi
                    continue

                # get next higher x
                xhi2 = xis[ihi+1]
                vhi2 = vis[ihi+1]
                dhi2 = xhi2 - xhi

                # xhi2 out of tolerance
                if dhi2 > tol:
                    vos[n] = vhi
                    continue

                # extrapolate
                vos[n] = f(xhi, xhi2, vhi, vhi2, xo)
                continue

        ###################################################
        # OUT OF BOUNDS ON UPPER SIDE                     #
        ###################################################
        if np.isfinite(dlo) and ~np.isfinite(dhi):
            # constant value
            if oob == 'val':
                vos[n] = oobv
                continue

            # nearest neighbour
            elif oob == 'nn':
                vos[n] = vlo
                continue

            # extrapolate
            elif oob  == 'ex':
                # lower bound
                if ilo == 0:
                    vos[n] = vlo
                    continue

                # get next higher x
                xlo2 = xis[ilo-1]
                vlo2 = vis[ilo-1]
                dlo2 = xlo - xlo2

                # xhi2 out of tolerance
                if dlo2 > tol:
                    vos[n] = vlo
                    continue

                # extrapolate
                vos[n] = f(xlo, xlo2, vlo, vlo2, xo)
                continue
    chrono.resumee()

    return vos
        
def interpolate_1d_linear_old(
        x_in,
        val_in,
        x_out,
        x_tolerance=0,
        out_of_bounds='nan',
        assume_sorted=True,
        many_nans=False,
        ):
    """Linear 1D-interpolation for numbers or datetime.

        x_in and x_out can be floats or datetime.datetime instances.
       
        Parameters
        ----------
        x_in : array of floats or list of datetime.datetime
        x_out : array of floats or list of datetime.datetime
        val_in : array of floats
            same length as x_in
        x_tolerance : float or datetime.timedelta, optional
            only nearest neighbours that are close than x_tolerance are
            considered. 0 causes infinite tolerance! Default: 0.
        out_of_bounds : {'zero', 'nan', 'nearest', float}
            out of bounds value are replaced by this
        assume_sorted : bool
            if True, the function assumes that x_in be sorted in rising order
        many_nans : bool
            Switch this to True, if either x_in or val_in contain a lot of
            nan's.  The function yields exactly the same result, regardless on
            whether many_nans is True or False. It is just a matter of
            performance.  If x_in or val_in contain a lot of nan's, switching
            this to True will speed up the function. If they don't, it will
            slow down the function.
        
        Returns
        -------
        val_out : array
            The output is a linear interpolation between the sample values
            given in values_sample at the times given in val_in.     

        Notes
        -----
        float vs datetime.datetime:
            If some of x_in, x_out, or x_tolerance are given in numerical
            values and others as instances of datetime.datetime (or
            datetime.timedelta in the case of x_tolerance), the numericals are
            considered as seconds (since 1970).

        x_tolercance:
            if 0, this is interpreted as infinite tolerance!
       
        Author
        ------
        Written in 2014-2016
        by Andreas Anhaeuser
        Insitute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """
    
    #### CHECK INPUT ###
    if not len(x_in) == len(val_in):
        raise LookupError('x_in and val_in must be of the same length.')

    #### Unsorted
    # recursively call function if unsorted:
    if not assume_sorted:
        idx_old = np.argsort(x_in)
        xi_sor  = sorted(x_in)
        vi_sor  = [val_in[i] for i in idx_old]
        
        o       = out_of_bounds
        return interpolate_1d_linear(x_in = xi_sor,
            val_in = vi_sor, x_out = x_out, x_tolerance = x_tolerance, 
            out_of_bounds = o, assume_sorted = True)
        
    # this is for dealing with out of bounds values lateron:
    if out_of_bounds in ['nan', 'NaN']:
        oob  = 'val'
        oobv = np.nan
    elif out_of_bounds in ['zero', 'zeros', '0']:
        oob  = 'val'
        oobv = 0.
    elif isinstance(out_of_bounds, numbers.Number):
        oob  = 'val'
        oobv = out_of_bounds
    else:
        oob  = 'nn'
        oobv = None
    
    # special case: empty input array:
    if len(x_in) == 0:
        return np.array([oobv] * len(x_out))
    
    # special case: empty output array:
    if len(x_out) == 0:
        return np.array([])
    
    # copy arrays and convert to np.array:
    xi  = x_in[:]
    xo  = x_out[:]
    vi  = val_in[:]
    tol = x_tolerance
    
    # convert datetime_list into a numerical list:
    for c in [xi, xo, vi]:
        if isinstance(c[0], dt.datetime):
            c[:] = du.datetime_to_seconds(c)

    # convert numpy arrays to lists:
    def convert(x):
        if isinstance(x, np.ndarray):
            return x.tolist()
        else:
            return x
    xi = convert(xi)
    xo = convert(xo)
    vi = convert(vi)

    # convert tol from timedelta to numerical:
    if isinstance(tol, dt.timedelta):
        tol = tol.seconds

    # convert tol==0 to inf:
    if tol == 0:
        tol = np.inf

    #### NaN's
    # If vi contains a lot of nan's, this will slow down nearest_neighbour.
    # For this reason, delete these elements from xi and vi:
    if many_nans:
        n = 0
        numbers = []
        for n in range(len(xi)):
            if not (np.isnan(vi[n]) or np.isnan(xi[n])):
                numbers.append(n)
        xi = [xi[n] for n in numbers]
        vi = [vi[n] for n in numbers]
            
    #### INTERPOLATION    
    # interpolation function:
    def f(xlo, xhi, vlo, vhi, x):
        """Linear interpolation function."""
        return ((x-xlo) * vhi + (xhi-x) * vlo) / (xhi-xlo)
        
    # initialization:
    valout = []
    N = len(xo)
    n = 0
    for x in xo:
        n += 1
        ###
        # FIND NEAREST NEIGHBOURS AND DETERMINE DISTANCES
        # nearest neighbours:
        xlo, ilo, vlo = nearest_neighbour(
            xi, x,
            direction='down', values=vi, skip_nans=True, assume_sorted=True)
            
        xhi, ihi, vhi = nearest_neighbour(
                xi, x, direction='up', values=vi,
                skip_nans=True, assume_sorted=True)


        # if no nearest neighbour is found:
        if xlo==None: xlo = -np.inf
        if xhi==None: xhi = +np.inf
        if vlo==None: vlo =  np.nan
        if vhi==None: vhi =  np.nan
        
        # distances:
        distlo = x-xlo
        disthi = xhi-x
        dist   = [distlo, disthi]
      
        ###
        # DETERMINE THE INTERPOLATED VALUE v CORRESPONDING TO x
        # special case: x is exactly on a non-nan input point:
        if max(dist) == 0:
            v = vlo   # (vlo==vhi in this case)
        
        # finite tolerance:
        elif tol < np.inf:
            inbounds = (distlo <= tol, disthi <= tol)
            
            if inbounds == (True, True):
                v = f(xlo, xhi, vlo, vhi, x)
            elif inbounds == (True, False):
                v = vlo
            elif inbounds == (False, True):
                v = vhi
            elif inbounds == (False, False) and oob == 'val':
                v = oobv
            elif inbounds == (False, False) and oob == 'nn':
                if distlo <= disthi:
                    v = vlo
                else:
                    v = vhi
        
        # infinite tolerance:
        elif tol == np.inf:
            neighbours = (distlo<np.inf,disthi<np.inf)
            
            if neighbours==(True, True):
                v = f(xlo,xhi,vlo,vhi,x)
            elif False in neighbours and oob=='val':
                v = oobv
            elif False in neighbours and oob=='nn':
                if distlo<=disthi: v = vlo
                else:              v = vhi
                

        # EXTEND OUTPUT ARRAY BY v:  
        valout.append(v)
        
    return np.array(valout)

def interpolate_2d_linear(*args, **kwargs):
    return interpolate_2d_linear_v2(*args, **kwargs)
        
def interpolate_2d_linear_v1(
        x_in,
        y_in,
        samples,
        x_out,
        y_out,
        x_tolerance=0.,
        y_tolerance=0.,
        out_of_bounds='ext',
        info_on_screen=True,
        prefix='interpolate_2D_linear: ',
        ):
    """Linear 2D-interpolation for numericals or datetime.
    
        Parameters
        ----------
        x_in : array (1d) of floats or datetime, length N
        y_in : array (1d) of floats or datetime, length M
        samples : array (shape N, M) of floats
        x_out, y_out : arrays (1d)
            the grid on which samples will be interpolated. If the elements of
            x_in are instances of datetime, the those of x_out should be, too.
            The same holds for y_in and y_out.
        x_tolerance : float or datetime.timedelta, optional
            only nearest neighbours that are close than x_tolerance are
            considered.  0 causes infinite tolerance! Default: 0.
        y_tolerance : float or datetime.timedelta, optional
            only nearest neighbours that are close than y_tolerance are
            considered.  0 causes infinite tolerance! Default: 0.
        out_of_bounds : {'zero', 'nan', 'nearest', 'ext', float}
                out of bounds value are replaced by this 'ext' : extrapolate
   
        Returns
        -------
        val_out : array
            The output is a linear interpolation between the sample values
            given in values_sample at the times given in val_in.     

        Notes
        -----
        float vs datetime.datetime:
            If some of x_in, x_out, or x_tolerance are given in numerical
            values and others as instances of datetime.datetime (or
            datetime.timedelta in the case of x_tolerance), the numericals are
            considered as seconds (since 1970).

        x_tolercance:
            if 0, this is interpreted as infinite tolerance!
       
        Author
        ------
        Written in 2016
        by Andreas Anhaeuser
        Insitute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """
    N = len(x_in)
    M = len(y_in)

    I = len(x_out)
    J = len(y_out)

    intermediate = np.nan * np.empty((I, M))
    out = np.nan * np.empty((I, J))

    # interpolate in x direction:
    for m in range(M):
        intermediate[:, m] = interpolate_1d_linear(
                x_in=x_in,
                val_in=samples[:, m],
                x_out=x_out,
                x_tolerance=x_tolerance,
                out_of_bounds=out_of_bounds,
                assume_sorted=True,
                )

    # interpolate in y direction:
    for i in range(I):
        out[i] = interpolate_1d_linear(
                x_in=y_in,
                val_in=intermediate[i, :],
                x_out=y_out,
                x_tolerance=y_tolerance,
                out_of_bounds=out_of_bounds,
                assume_sorted=True,
                )
    return out
    
def interpolate_2d_linear_v2(
        x_in,
        y_in,
        samples,
        x_out,
        y_out,
        x_tolerance=0.,
        y_tolerance=0.,
        out_of_bounds = 'nan',
        info_on_screen=True,
        prefix='interpolate_2D_linear: ',
        assume_sorted=False,
        dtype=float,
        ):
    """Linear 2D-interpolation for numericals or datetime.
    
        Parameters
        ----------
        x_in, y_in
            the sample grid. They can both either be 1D-arrays of numercials or
            of instances of datetime
        samples: 2d-array
            the function values on the (x_in,y_in)-grid. Must be given in the
            form samples[x_in, y_in]
        x_out, y_out
            the grid on which samples will be interpolated. If the elements of
            x_in are instances of datetime, the those of x_out should be, too.
            The same holds for y_in and y_out.
        out_of_bounds : {'nan', 'zero'}
            if x_out or y_out are outside the input grid, then the output will
            be assigned this value
        dtype
            datatype of the output array.
       
        Returns
        -------
        The function returns a 2D-array which is a linear interpolation
        between the sample values in the form samples_out[x_out,y_out]
        
        To do
        -----
        implement `assume_sorted`
       
        Written in 2014
        by Andreas Anhaeuser
        Insitute for Geophysics and Meteorology
        University of Cologne
        Germany
        <andreas.anhaeuser.data_analyst@posteo.net>
    """
    #### INPUT CHECK ###
    if not (len(x_in), len(y_in)) == np.shape(samples):
        raise AssertionError('Dimensions of x_in, y_in and samples do ' + \
            'not match.')
    xi = np.array(x_in)
    yi = np.array(y_in)
    xo = np.array(x_out)
    yo = np.array(y_out)
    zi = np.array(samples)
    
    # convert datetime_list into a numerical list:
    for c in [xi, yi, xo, yo]:
        if c[0].__class__ == dt.datetime:
            c[:] = np.array(du.datetime_to_seconds(c))

    # 2D linear interpolation function:
    f = interpolate_2d_one_point

    # out of bounds value:
    if out_of_bounds in ['zero', 'Zero', 'zeros', 'Zeros', '0', 0]:
        oob = 0
    elif out_of_bounds in ['nan', 'nans', 'NaN', 'NaNs', np.nan]:
        oob = np.nan

    X = len(xo)
    Y = len(yo)
    zo = np.zeros([X, Y], dtype=dtype)

    Ntotal = X * Y
    header = 'interpolate 2D'
    silent = not info_on_screen
    with Chronometer(Ntotal, header=header, silent=silent) as chrono:
        if silent:
            chrono.exit()

        for i, j in itertools.product(range(X), range(Y)):
            if not silent:
                chrono.show().loop()

            z = f(xi, yi, zi, xo[i], yo[j], oob, assume_sorted=assume_sorted)
            zo[i, j] = dtype(z)
                
    return zo 
    
def interpolate_2d_to_1d(
        x_in,
        y_in,
        samples,
        x_out,
        y_out,
        x_tolerance=0.,
        y_tolerance=0.,
        out_of_bounds = 'nan',
        info_on_screen=True,
        prefix='interpolate_2D_linear: ',
        assume_sorted=False,
        dtype=float,
        ):
    """Linear 2D-interpolation for numericals or datetime.
    
        Parameters
        ----------
        x_in, y_in: the sample grid. They can both either be 1D-arrays of
            numercials or of instances of datetime
        samples: 2D-array of the function values on the (x_in,y_in)-grid.
            Must be given in the form samples[x_in, y_in]
        x_out, y_out: the grid on which samples will be interpolated. If the
            elements of x_in are instances of datetime, the those of x_out
            should be, too. The same holds for y_in and y_out.
        out_of_bounds ('nan' or 'zero'): if x_out or y_out are outside the
            input grid, then the output will be assigned this value
        dtype : datatype of the output array.
       
        Returns
        -------
        The function returns a 2D-array which is a linear interpolation
        between the sample values in the form samples_out[x_out,y_out]
        
        NOTE: implementation assume_sorted to be done!
       
        Written in 2014
        by Andreas Anhaeuser
        Insitute for Geophysics and Meteorology
        University of Cologne
        Germany
        <andreas.anhaeuser.data_analyst@posteo.net>
    """
    #### INPUT CHECK ###
    if not (len(x_in), len(y_in)) == np.shape(samples):
        raise AssertionError('Dimensions of x_in, y_in and samples do ' + \
            'not match.')
    xi = np.array(x_in)
    yi = np.array(y_in)
    xo = np.array(x_out)
    yo = np.array(y_out)
    zi = np.array(samples)
    
    # convert datetime_list into a numerical list:
    for c in [xi, yi, xo, yo]:
        if c[0].__class__ == dt.datetime:
            c[:] = np.array(du.datetime_to_seconds(c))

    # out of bounds value:
    if out_of_bounds in ['zero', 'Zero', 'zeros', 'Zeros', '0', 0]:
        oob = 0
    elif out_of_bounds in ['nan', 'nans', 'NaN', 'NaNs', np.nan]:
        oob = np.nan

    N = len(xo)
    assert len(yo) == N
    zo = oob * np.ones(N, dtype=dtype)

    f = interpolate_2d_one_point
    for n in range(N):
        z = f(xi, yi, zi, xo[n], yo[n], oob, assume_sorted=assume_sorted)
        zo[n] = dtype(z)

    return zo 
    
def interpolate_2d_linear_v1(
        x_in,
        y_in,
        samples,
        x_out,
        y_out,
        dtype=np.float,
        out_of_bounds = 'nan',
        info_on_screen=True,
        prefix='interpolate_2D_linear: ',
        ):
    """Linear 2D-interpolation for numericals or datetime.
    
    Parameters
    ----------
    * x_in, y_in: the sample grid. They can both either be 1D-arrays of
      numercials or of instances of datetime
    * samples: 2D-array of the function values on the (x_in,y_in)-grid.
      Must be given in the form samples[x_in, y_in]
    * x_out, y_out: the grid on which samples will be interpolated. If the
      elements of x_in are instances of datetime, the those of x_out
      should be, too. The same holds for y_in and y_out.
    * out_of_bounds ('nan' or 'zero'): if x_out or y_out are outside the
      input grid, then the output will be assigned this value
    * dtype : datatype of the output array.
   
    Returns
    -------
    The function returns a 2D-array which is a linear interpolation
    between the sample values in the form samples_out[x_out,y_out]
    
    NOTE: implementation assume_sorted to be done!
   
    Written in 2014
    by Andreas Anhaeuser
    Insitute for Geophysics and Meteorology
    University of Cologne
    Germany
    <andreasfrederik.anhaeuser@smail.uni-koeln.de> """
    #### INPUT CHECK ###
    if not (len(x_in), len(y_in)) == np.shape(samples):
        raise AssertionError('Dimensions of x_in, y_in and samples do ' + \
            'not match.')
    xi = np.array(x_in)
    yi = np.array(y_in)
    xo = np.array(x_out)
    yo = np.array(y_out)
    zi = np.array(samples)
    
    # convert datetime_list into a numerical list:
    for c in [xi, yi, xo, yo]:
        if c[0].__class__ == dt.datetime:
            c[:] = np.array(du.datetime_to_seconds(c))

    # 2D linear interpolation function:
    def f(xin, yin, zin, x, y, oob):
        if x < min(xin) or x > max(xin) or y < min(yin) or y > max(yin):
            return oob        
        
        x1, idx1 = nearest_neighbour(xin, x, direction='down')
        x2, idx2 = nearest_neighbour(xin, x, direction='up'  )
        y1, idy1 = nearest_neighbour(yin, y, direction='down')
        y2, idy2 = nearest_neighbour(yin, y, direction='up'  )
        z11 = zin[idx1, idy1]
        z12 = zin[idx1, idy2]
        z21 = zin[idx2, idy1]
        z22 = zin[idx2, idy2]
        
        # avoid problems if D == 0:
        if x1 == x2 and y1 == y2:
            result = z11
        elif x1 == x2:
            # interpolate only in y-direction
            result = ((y-y1)*z12 + (y2-y)*z11)/(y2-y1)
        elif y1 == y2:
            # interpolate only in x-direction
            result = ((x - x1) * z21 + (x2 - x) * z11) / (x2  -x1)
        else:
            # 2D interpolation:
            D = (x2 - x1) * (y2 - y1)
            result =  (1/D * (z11 * (x2 - x)  *(y2 - y) + \
                              z21 * (x - x1) * (y2 - y) + \
                              z12 * (x2 - x) * (y - y1) + \
                              z22 * (x - x1) * (y - y1) \
                              )  \
                       )
        return result
               
    # out of bounds value:
    if out_of_bounds in ['zero', 'Zero', 'zeros', 'Zeros', '0', 0]:
        oob = 0
    elif out_of_bounds in ['nan', 'nans', 'NaN', 'NaNs', np.nan]:
        oob = np.nan

    X = len(xo)
    Y = len(yo)
    zo = np.zeros([X, Y], dtype=dtype)

    performance_info = PI((X*Y), prefix=prefix)
    for i in range(X):
        for j in range(Y):
            if info_on_screen:
                performance_info.loop_and_show()

            zo[i, j] = dtype(f(xi, yi, zi, xo[i], yo[j], oob))
                
    return zo 

def intp_1d_get_indices_and_weights():
    raise NotImplementedError

def interpolate_2d_one_point(xin, yin, zin, x, y, oob, assume_sorted=False):
    if x < min(xin) or x > max(xin) or y < min(yin) or y > max(yin):
        return oob        
    
    kwargs = {'assume_sorted' : assume_sorted}
    x1, idx1 = nearest_neighbour(xin, x, direction='down', **kwargs)
    x2, idx2 = nearest_neighbour(xin, x, direction='up', **kwargs)
    y1, idy1 = nearest_neighbour(yin, y, direction='down', **kwargs)
    y2, idy2 = nearest_neighbour(yin, y, direction='up', **kwargs)
    z11 = zin[idx1, idy1]
    z12 = zin[idx1, idy2]
    z21 = zin[idx2, idy1]
    z22 = zin[idx2, idy2]

    # handle None
    # =================================================
    if x1 is None and x2 is None:
        return oob

    if y1 is None and y2 is None:
        return oob

    if x1 is None:
        x1 = x2
        idx1 = idx2

    if x2 is None:
        x2 = x1
        idx2 = idx1

    if y1 is None:
        y1 = y2
        idy1 = idy2

    if y2 is None:
        y2 = y1
        idy2 = idy1
    # =================================================
    
    # avoid problems if D == 0:
    if x1 == x2 and y1 == y2:
        return z11

    if x1 == x2:
        # interpolate only in y-direction
        return ((y-y1)*z12 + (y2-y)*z11)/(y2-y1)
    if y1 == y2:
        # interpolate only in x-direction
        return ((x - x1) * z21 + (x2 - x) * z11) / (x2  -x1)

    # 2D interpolation
    # --------------------------------------------
    # weights
    w11 = (x2 - x) * (y2 - y)
    w21 = (x - x1) * (y2 - y)
    w12 = (x2 - x) * (y - y1)
    w22 = (x - x1) * (y - y1)
    total_weights = (x2 - x1) * (y2 - y1)

    # weighted contributions
    # ll, lr, ul, ur : lower/upper left/right
    ll = z11 * w11
    lr = z21 * w21
    ul = z12 * w12
    ur = z22 * w22 

    # weighted mean
    result = (ll + lr + ul + ur) / total_weights
    
    return result
               
    
######################################
# INTERPOLATION OF data DICTIONARIES #
######################################
def temp_interpolate(
        data,
        time,
        time_tol,
        interpolate=True,
        time_key='time', 
        ignore_keys=[],
        messages=True,
        stop_on_warning=True,
        ):
    """Perform temporal linear interpolation of all time-dependent variables.

    Modifies the input dictionary. Interpolates all entries that have an
    axis with length equal to the time_key entry in data, except for those
    entries specified in ignore_keys.
    Exception: If time_tol is None, no interpolation is performed.

    Parameters
    ----------
    data : dict
        must contain at least the entry datetime
    time : datetime.datetime object
        the time on which the data will be interpolated
    time_tol : datetime.timedelta object or None
        time_tolerance, neighbours are only considered, if they are within
        this range around time.
    interpolate : bool, optional
        if False, the value of the nearest neighbour is taken, otherwise
        values are interpolated between the two neighbouring time steps.
        Default: True
    time_key : str, optional
        the name of the time variable in data. Default: 'datetime'
    ignore_keys : list of str
        these entries of data will not be interpolated even if they have
        axis that has a length equal to the time dimension. Default: []
    messages : bool
        show messages
    stop_on_warning : bool
        wait for confirmation when warning occurs.

    Returns
    -------
    int
        8 bit flags. 0 meens ok. 1 is an error or a warnung
        (0) If time_tol is None, the dictionary remains unchanged and
            0 is returned.
        (0) If time_tol is not None and the interpolation is successful,
            the dictionary variables are interpolated and 0 is returned.
        (1) If an ambiguous axis case occurs, bit 0 is set to 1.
        (128) If unsuccessful, then variables with a time axis are set to nan
              and bit 7 is set to 1.

    What happens to time?
    ---------------------
    The dimension and value(s) of data['time'] after to operation depend
    on the values of the parameters time_tol and interpolate:

        time_tol -----(None?)-----> RETURN ALL MEASUREMENT --(NO DATA?)---+
        |                             TIMES OF THE DAY                    |
        |                           (data['time'] is a list               |
        |                            of datetime.datetime)                |
        (not None?)                                                       |
        |                                                                 |
        interpolate --(True?)--> INTERPOLATE BETWEEN TWO NEIGHBOURS       |
        |                           IF BOTH ARE WITHIN time_tol           |
        (False?)                  (data['time'] is a dt.datetime)         |
        |                                          |                      |
        V                                          |                      | 
        USE VALUE OF NEAREST NEIGHBOUR             |                      |
           IF IT IS WITHIN time_tol                |                      |
        (data['time'] is a dt.datetime)            |                      |
                               |                   |                      |
                          (NO NEIGHBOUR WITHIN time_tol?)                 |
                                           |                              |
                                           +--> RETURN False <------------+
                                               
                                        +------------------+
                                        | ELSE RETURN True |
                                        | (REGULAR CASE)   |
                                        +------------------+

    Bugs
    ----
    Moderately tested. Not unlikely to contain some bugs. (AA, 2015-10-10)

    Author
    ------
    Written in 2015
    by Andreas Anhaeuser
    Insitute for Geophysics and Meteorology
    University of Cologne
    Germany
    <anhaeus@meteo.uni-koeln.de>
    """
    error = 0
    if time_tol is None:
        if data is not None:
            return True
        else:
            return False
    
    dtime = data[time_key]
    T = len(dtime)

    ###################################
    # NEIGHBOURS                      #
    ###################################
    (timelo, tlo) = nearest_neighbour(dtime, time, 'd')
    (timehi, thi) = nearest_neighbour(dtime, time, 'u')

    ###################################
    # WEIGHTING COEFFICIENTS          #
    ###################################
    a, b = time_weighting_coefficients(time, time_tol, timelo, timehi)
    if a is None and b is None:
        if messages:
            str_utils.printc('INFO: temporally_interpolate(): ' + 
                'No values found within time tolerance.', 'i')
        error |= 128
        return error

    ###################################
    # CASE: NEAREST NEIGHBOUR         #
    ###################################
    if not interpolate:
        if a < b:
            a, b = None, 1.
            data[time_key] = timehi
        else:
            a, b = 1., None
            data[time_key] = timelo

    ###################################
    # INTERPOLATE                     #
    ###################################
    if interpolate:
        data[time_key] = time
    for key in data.keys():
        if key in ignore_keys:
            continue
        shape = np.shape(data[key])
        if T not in shape:
            continue
        axis = shape.index(T)

        # check whether axis is ambiguous:
        if T in shape[axis+1:]:
            error |= 1  # change last bit to 1
            str_utils.printc('WARNING: Axis in %s ambiguous!' % key, 'w')
            if stop_on_warning:
                inp = raw_input('Press [ENTER] to continue, ' + \
                        ' type E + [ENTER] to raise an Exception...')
                if inp.lower() == 'e':
                    raise Exception('Axis ambiguous.')

        # low value
        if a:
            datalo = a * np.take(data[key], tlo, axis)
        else:
            datalo = 0

        # high value
        if b:
            datahi = b * np.take(data[key], thi, axis)
        else:
            datahi = 0

        data[key] = datalo + datahi

    return error

def time_weighting_coefficients(time, time_tol, timelo, timehi):
    """Return a pair of floats.

    Helping function to temp_interpolate. timelo and/or timehi may be None. For
    more details on the behavior of the function see the Return section.

    Parameters
    ----------
    time : dt.datetime
        must be between timelo and timehi
    time_tol : dt.timedelta
        must be non-negative
    timelo : dt.datetime or None
        lower neighbour
    timehi : dt.datetime or None
        upper neighbour

    Returns
    -------
    float, float
    float, None
    None, float
    None, None
        The regular case are two floats between 0 and 1, which are the
        weighting coefficients for the lower and higher neighbour,
        respectively. If one or both of them are non-existent or out of bounds
        (with respect to time_tol), the respective coefficient is None.
    """
    # check input
    if timelo is not None and timelo > time:
        raise ValueError('timelo must not be larger than time.')
    if timehi is not None and timehi < time:
        raise ValueError('timehi must not be smaller than time.')

    # which neighbours to use:
    lo, hi = True, True

    # check whether to use low neighbour
    if timelo is None:
        lo = False
    else:
        dello = time - timelo
        if dello < dt.timedelta():
            raise ValueError('time must not be smaller than timelo.')
        if dello > time_tol:
            lo = False

    # check whether to use high neighbour
    if timehi is None:
        hi = False
    else:
        delhi = timehi - time
        if delhi < dt.timedelta():
            raise ValueError('time must not be larger than timehi.')
        if delhi > time_tol:
            hi = False

    # coefficients
    if not lo and not hi:
        return None, None
    if not lo:
        return None, 1.
    if not hi:
        return 1., None
    a = delhi.total_seconds() / (dello + delhi).total_seconds()
    b = dello.total_seconds() / (dello + delhi).total_seconds()
    return a, b

#####################################
# MISCELLANEOUS                     #
#####################################
def convert_to_nan(data, key=None):
    """Convert data to NaN if key applies.

    Parameters
    ----------
    data : any object
    key  : a function that returns a boolean for input objects of which the
           data type is that of data

    Returns
    -------
    If data is not iterable, it returns data if key(data) yields False.
    Otherwise it returns a nan.

    If data is iterable, the output is a list.
    """
    if key is None:
        key = lambda x: np.isnan(x) or x < 0
    if isinstance(data, collections.Iterable):
        return [convert_to_nan(d, key) for d in data]
    return np.nan if key(data) else data

 
if __name__ == "__main__":
    N = 10000
    T = N/5.
    noise_level = 0.2
    hw = N/100
    per = T*2

    xin = np.arange(N)
    yin = np.sin(xin * 2*np.pi/T) + noise_level * np.random.random(N)


    yout, xout = running_mean(xin, yin, halfwidth=hw, periodicity=per,
            assume_equidistant=True, assume_sorted=True,
            contains_nans=True,)

    plt.plot(xin, yin, 'b-')
    plt.plot(xout, yout, 'r-')
    plt.show()
