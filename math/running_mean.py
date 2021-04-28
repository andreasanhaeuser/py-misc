# standard modules
from copy import deepcopy as copy
import warnings

# PyPI modules
import numpy as np

# misc
from misc.math.nearest_neighbour import nearest_neighbour
# from misc.chronometer import Chronometer

################################################################
# main                                                         #
################################################################
def running_mean(x, y, halfwidth, axis=0, assume_equidistant=False, **kwargs):
    """Return running mean an the corresponding x-positions.
    
        Parameters
        ----------
        x : flat array
            positions of y.
        y : array
            the values at the x-positions to interpolate. The length of `x`
            must be the same as the dimension of `y` along `axis`.
        halfwidth : float > 0 
            half width of the averaging window
        axis : int, optional
            The axis along which to take the mean. If not equidistant, the
            function only works on flat arrays. Default: 0
        include_left_bound : bool, optional
            use values at position x[i] - halfwidth for averaging or not.
            Default: True
        include_right_bound : bool, optional
            as include_left_bound. Default: False
        assume_equidistant : bool, optional
            set this to True if x is equidistant. Default: False
        assume_sorted : bool, optional
            set this to True if x is sorted. Default: False
        periodicity : {float, np.inf}, optional
            non-positive values (including 0) are interpreted as np.inf.
            Default: np.inf
        contains_nans : bool, optional
            Set this False, if you are sure that y does not contains nan's.
            This will speed up the function. Default: True
                
        Returns
        -------
        ymean : array
            the running mean of y
        x_out : array, same shape
            the x values of the ymean
        
        Note
        ----
        The values in ymean correspond to x_out. x_out may be different
        to x. This is the case if one or more of the following is true:
            a) periodicity is positive and finite
            b) assume_sorted is False
        
        Performance
        -----------
        The following will usually significantly speed up the function:
            a) assume_equidistant is True
            b) assume_sorted is True
            c) contains_nans is False (especially in combination with a)

        Bugs
        ----
        In the non-equadistant case, the include_bounds parameters do not work
        and intended.

        The multi-axes versions does not work.
        
        Author
        ------
        Written in 2015-2016
        by Andreas Anhaeuser
        Insitute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """ 
    # multi-axes version does not work
    if len(np.shape(y)) > 1:
        raise NotImplementedError('Sorry, the multiple axes version does ' + \
                'not work. Feel free to fix it.')

    ###################################################
    # ROLL AXIS                                       #
    ###################################################
    if axis == 0:
        y0 = y
    else:
        assert assume_equidistant
        y0 = np.rollaxis(y, axis, 0)

    assert len(x) == len(y0)
    y0m, xm = running_mean_axis_0(x, y0, halfwidth,
            assume_equidistant=assume_equidistant, **kwargs)

    if axis == 0:
        ym = y0m
    else:
        ym = np.rollaxis(y0m, 0, axis + 1)

    return ym, xm

def running_mean_2d(
        x, y, z, hw_x, hw_y,
        assume_equidistant_x=False,
        assume_equidistant_y=False,
        periodicity_x=np.inf,
        periodicity_y=np.inf,
        assume_sorted_x=False,
        assume_sorted_y=False,
        contains_nans=True,
        info_on_screen=True,
        prefix='running mean 2d: ',
        ):
    """Works as running_mean() but in 2 dimensions.
    
        x and y can be given as vectors.
        z must be a 2d-matrix and its shape must be (len(x), len(y))
            
        Author
        ------
        Written in 2015
        by Andreas Anhaeuser
        Insitute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """ 
    #### INPUT CHECK ###
    I, J  = np.shape(z)
    if not len(x) == I:
        raise AssertionError('')
    if not len(y) == J:
        raise AssertionError('')

    performance_info = Chronometer(J, header='(ydir) ' + prefix)
    for j in range(J):
        if info_on_screen:
            performance_info.loop_and_show()

        rm = running_mean(
            x, z[:,j], hw_x, 
            assume_equidistant=assume_equidistant_x,
            assume_sorted=True,
            periodicity=periodicity_x,
            contains_nans=contains_nans,
            )
        if j == 0:
            x_new = rm[1]
            I = len(x_new)
            meanx = np.zeros([I, J])
            meanx[:,j] = rm[0]
        else:
            meanx[:,j] = rm[0]
            
    performance_info = Chronometer(I, header='(xdir) ' + prefix)
    for i in range(I):
        if info_on_screen:
            performance_info.loop_and_show()

        rm = running_mean(
            y, meanx[i,:], hw_y, 
            assume_equidistant=assume_equidistant_y,
            assume_sorted=True,
            periodicity=periodicity_y,
            contains_nans=contains_nans,
            )
        if i == 0:
            y_new = rm[1]
            J = len(y_new)
            mean = np.zeros([I, J])
            mean[i,:] = rm[0]
        else:
            mean[i,:] = rm[0]
            
    return (mean, x_new, y_new)

################################################################
# helper functions                                             #
################################################################
def running_mean_axis_0(
        x, y, halfwidth,
        include_left_bound=True,
        include_right_bound=False,
        assume_equidistant=False,
        assume_sorted=False,
        periodicity=np.inf,
        contains_nans=True,
        ):
    """Return running mean an the corresponding x-positions.
    
        Parameters
        ----------
        x : flat array
            positions of y.
        y : array
            the values at the x-positions to interpolate. The length of `x`
            must be the same as the dimension of `y` along the first axis.
        halfwidth : float > 0 
            half width of the averaging window
        include_left_bound : bool, optional
            use values at position x[i] - halfwidth for averaging or not.
            Default: True
        include_right_bound : bool, optional
            as include_left_bound. Default: False
        assume_equidistant : bool, optional
            set this to True if x is equidistant. Default: False
        assume_sorted : bool, optional
            set this to True if x is sorted. Default: False
        periodicity : {float, np.inf}, optional
            non-positive values (including 0) are interpreted as np.inf.
            Default: np.inf
        contains_nans : bool, optional
            Set this False, if you are sure that y does not contains nan's.
            This will speed up the function. Default: True
                
        Returns
        -------
        ymean : array
            the running mean of y
        x_out : array, same shape
            the x values of the ymean
        
        Note
        ----
        The values in ymean correspond to x_out. x_out may be different
        to x. This is the case if one or more of the following is true:
            a) periodicity is positive and finite
            b) assume_sorted is False
        
        Performance
        -----------
        The following will usually significantly speed up the function:
            a) assume_equidistant is True
            b) assume_sorted is True
            c) contains_nans is False (especially in combination with a)

        Bugs
        ----
        In the non-equadistant case, the include_bounds parameters do not work
        and intended.

        The multi-axes versions does not work.
        
        Author
        ------
        Written in 2015-2016
        by Andreas Anhaeuser
        Insitute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """ 

    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    assert 0 <= periodicity <= np.inf
    assert len(x) == len(y)

    ###################################################
    # CASE: UNSORTED                                  #
    ###################################################
    # (recursively call function with sorted arrays):
    if not assume_sorted:
        idx = np.argsort(x)
        inv_idx = np.argsort(idx)
        return running_mean_axis_0(
            x[idx], y[idx], halfwidth,
            include_left_bound=include_left_bound,
            include_right_bound=include_right_bound,
            assume_equidistant=assume_equidistant,
            assume_sorted=True,
            periodicity=periodicity,
            contains_nans=contains_nans,
            )

    ###################################################
    # ABBREVIATIONS                                   #
    ###################################################
    N = len(x)
    hw = abs(halfwidth)
    per = periodicity

    ###################################################
    # CASE: NON-EQUIDISTANT                           #
    ###################################################
    if not assume_equidistant:
        return running_mean_non_equidistant(
            x,
            y,
            halfwidth=hw,
            include_left_bound=include_left_bound,
            include_right_bound=include_right_bound,
            periodicity=per,
            contains_nans=contains_nans,
            )

    ###################################################
    # CASE: PERIODIC                                  #
    ###################################################
    periodic = (0 < per < np.inf)

    if periodic:
        # convert periodicity to number of indices
        periodic = True
        Nperiod = nearest_neighbour(x, x[0] + per)[1]
        
        #=============== FOLD ARRAYS ======================
        # case: contains nans 
        if contains_nans and Nperiod != N:
            y_new = np.zeros(Nperiod)
            count = np.zeros(Nperiod)
            for n in range(N):
                val = y[n]
                if np.isnan(val):
                    continue
                idx = n % Nperiod
                y_new[idx] += y[n]
                count[idx] += 1
            y = y_new / count
            x = x[:Nperiod]
            contains_nans = False

        # case: nan-free
        elif Nperiod != N:
            y_new = np.zeros(Nperiod)
            count = np.emtpy(Nperiod)
            beg = N % Nperiod
            count[:beg] = N//Nperiod + 1
            count[beg:] = N//Nperiod
            for n in range(N):
                y_new[n % Nperiod] += y[n]
            y = y_new / count
            x = x[:Nperiod]
        #==================================================

    ###################################################
    # CONVERT HALFWIDTH TO NUMBER OF ELEMENTS         #
    ###################################################
    xhi = x[0] + hw
    xn, n = nearest_neighbour(x, xhi, direction='d', assume_sorted=True)[:2]
    nlo = n
    nhi = n
    # include bounds or not:
    if xn == xhi:
        if not include_left_bound and nlo > 0:
            nlo -= 1
        if not include_right_bound and nhi > 0:
            nhi -= 1
        
    ###################################################
    # CALL SUB-FUNCTION                              #
    ###################################################
    mean = running_mean_equidistant(
                    y,
                    nlo=nlo,
                    nhi=nhi,
                    periodic=periodic,
                    contains_nans=contains_nans,
                    )

    return mean, x

def running_mean_non_equidistant(
        x,
        y,
        halfwidth,
        include_left_bound=True,
        include_right_bound=False,
        periodicity=np.inf,
        contains_nans=False,
        ):
    """*Return the running mean of a sorted but not equidistant array.
    
        -- HELPER FUNCTION -- 
        
        This is written to be a helper function for running_mean().
        
        Parameters
        ----------
        see running_mean()
        
        Returns
        -------
        as running_mean()
        
        Note
        ----
        If x is equidistant, running_mean_equidistant is much much faster!

        Bugs
        ----
        include_left_bound/include_right_bound don't work properly.  If
        periodicity hits exactly the x-spacing, the function yields wrong
        results.

        Author
        ------
        Written in 2015
        by Andreas Anhaeuser
        Insitute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """ 
    #### INPUT CHECK ###

    #### ABBREVIATIONS ###
    hw = halfwidth
    N = len(x)
    per = periodicity
    periodic = 0 < per < np.inf
    mean = np.zeros(N)
    imessage = 2**18
    x = np.array(x)
    y = np.array(y)
    
    #### DEFINE HELPER FUNCTIONS ###
    f = np.mean if not contains_nans else np.nanmean
    if include_left_bound:
        def increase_left(i, ilo):
            return (x[i] - x[ilo]) % per >= hw
    else:
        def increase_left(i, ilo):
            return (x[i] - x[ilo]) % per > hw
    if include_right_bound:
        def increase_right(i, ihi):
            return (x[(ihi+1) % N] - x[i]) % per <= hw
    else:
        def increase_right(i, ihi):
            return (x[(ihi+1) % N] - x[i]) % per < hw

    #### IF PERIODIC, THEN "FOLD" ARRAYS ###
    if periodic:
        x = x % per
        idx = np.argsort(x)
        x = x[idx]
        y = y[idx]

    #### COMPUTE ###
    if not periodic:
        ilo = 0
    else:
        ilo = nearest_neighbour(x, (x[0]-hw) % per, 'u')[1]
    ihi = 0
    if ilo is None:
        ilo = 0

    for i in range(N):
        if not i % imessage:
            pc = ('%1.1f' % (100.*(i+1.)/N)) + '%'
#            print(pc)        
        while increase_left(i, ilo):
            ilo = (ilo + 1) % N
        while increase_right(i, ihi):
            ihi = (ihi + 1) % N
        mean[i] = f(y[ilo:ihi+1])
    return (mean, x)

def running_mean_equidistant(
        x,
        nlo,
        nhi=None,
        periodic=False,
        contains_nans=True,
        ):
    """*Return the running mean of a sorted and equidistant array.  

        -- HELPER FUNCTION --
        This is written to be a helper function for running_mean().
        
        Parameters
        ----------
        x : array
            the values. operation is performed only on the first axis.
        nlo : int >= 0
            number of points to the left to consider. 0 means only the center
            point.
        nhi : int >= 0, optional.
            as nlo, but to the right. Default: same as nlo
        contains_nans : bool, optional
            turn it on, if you are not sure, otherwise you receive wrong
            results.  Default: True
        periodic : bool, optional
            Consider the edges of x the be connected or not. Default: False 

        Returns
        -------
        running_mean : array
            array of the same shape a x
        
        Complexity
        ----------
        Linear in all cases.

        Performance
        -----------
        contains_nans slows the function down by a factor of 2.
        periodic has no effect on the performance.

        Author
        ------
        Written in 2016
        by Andreas Anhaeuser
        Insitute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """ 
    ###################################################
    # DEFAULT                                         #
    ###################################################
    nhi = nlo if nhi is None else int(nhi)

    ###################################################
    # CALL SUB-FUNCTIONS                              #
    ###################################################
    if not contains_nans:
        return running_mean_equidistant_no_weights(x, nlo, nhi, periodic)

    else:
        nans = np.isnan(x)
        weights = ~nans
        xmod = copy(x)
        xmod[nans] = 0.
        return running_mean_equidistant_with_weights(
                xmod, nlo, nhi, weights, periodic)

def running_mean_equidistant_no_weights(
        x,
        nlo,
        nhi,
        periodic=False,
        axis=0,
        ):
    """*Helper function.
    
    Parameters
    ----------
    x : array
        the operation is only performed along the first axis.
    nlo : int >= 0
        number of points to the left to consider. 0 means only the center
        point.
    nhi : int >= 0
        as nlo, but to the right.
    periodic : bool, optional
        Consider the edges of x the be connected or not. Default: False 

    Returns
    -------
    xmean : array, same shape as x
        if not periodic, a lower number of points is taken into account on the
        edges (nlo is smaller on the lower edge, nhi get smaller on the top
        edge).
    """
    ############################################################
    # axis != 0                                                #
    ############################################################
    if axis != 0:
        self = running_mean_equidistant_no_weights
        x_swapped = np.moveaxis(x, axis, 0)
        x_mean_swapped = self(x_swapped, nlo, nhi, periodic, axis=0)
        x_mean = np.moveaxis(x_mean_swapped, 0, axis)
        return x_mean

    shape = np.shape(x)
    N = shape[axis]
    n = nlo + nhi + 1
    if not n <= N:
        text = 'Averaging window larger than input values: %i > %i' % (n, N)
        raise ValueError(text)

    shape_zeros = (1,) + shape[1:]
    zeros = np.zeros(shape_zeros)

    # broadcast range
    range_lo = np.arange(n-nlo, n)
    range_hi = np.arange(n-1, n-nhi-1, -1)
    axes = range(len(shape))
    for axis in axes[1:]:
        range_lo = np.expand_dims(range_lo, axis)
        range_lo = np.repeat(range_lo, shape[axis], axis)
        range_hi = np.expand_dims(range_hi, axis)
        range_hi = np.repeat(range_hi, shape[axis], axis)


    if not periodic:
        xm = np.nan * np.ones_like(x)
        xc = np.cumsum(np.insert(x, 0, zeros, 0), dtype=float, axis=0)

        # lower edge:
        xm[:nlo] = xc[n-nlo : n] / range_lo
        # center:
        xm[nlo : -nhi] = (xc[n:] - xc[:-n]) / n
        # upper edge:
        xm[-nhi:] = (xc[N] - xc[-n : -nlo-1]) / range_hi
        return xm

    else:
        xp = 1. * np.concatenate(([0.], x[-nlo:], x, x[:nhi]), 0)
        xc = np.cumsum(xp, 0)
        xm = (xc[n:] - xc[:-n]) / n
        return xm

def running_mean_equidistant_with_weights(
        x,
        nlo,
        nhi,
        w,
        periodic=False,
        contains_nans=True,
        ):
    """*Helper function.
    
    Parameters
    ----------
    x : array
        the operation is only performed along the first axis.
    nlo : int >= 0
        number of points to the left to consider. 0 means only the center
        point.
    nhi : int >= 0
        as nlo, but to the right.
    w : array, same shape as x
        the weight of each value in x to the running mean
    periodic : bool, optional
        Consider the edges of x the be connected or not. Default: False 

    Returns
    -------
    xmean : array, same shape as x
        if not periodic, a lower number of points is taken into account on the
        edges (nlo is smaller on the lower edge, nhi get smaller on the top
        edge).
    """
    N = len(x)
    n = nlo + nhi + 1
    assert n <= N

    if contains_nans:
        idx_nan = np.isnan(x)
        x[idx_nan] = 0
        w[idx_nan] = 0


    ############################################################
    # periodic                                                 #
    ############################################################
    if periodic:
        # weights:
        wp = 1. * np.concatenate(([0.], w[-nlo:], w, w[:nhi]))
        wc = np.cumsum(wp)

        # weighted x:
        xw = x * w 
        xwp = 1. * np.concatenate(([0.], xw[-nlo:], xw, xw[:nhi]))
        xwc = np.cumsum(xwp)

        xm = (xwc[n:] - xwc[:-n]) / (wc[n:] - wc[:-n])
        return xm

    ############################################################
    # non-periodic                                             #
    ############################################################
    xm = np.empty(N)

    # weights:
    wc = np.cumsum(np.insert(w, 0, 0), dtype=float)

    # weighted x:
    xw = x * w 
    xwc = np.cumsum(np.insert(xw, 0, 0), dtype=float)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)

        # lower edge:
        lo = n - nlo
        hi = n
        xm[:nlo] = xwc[lo:hi] / wc[lo:hi] 

        # center:
        xm[nlo : -nhi] = (xwc[n:] - xwc[:-n]) / (wc[n:] - wc[:-n])

        # upper edge:
        lo = - (n + 1)
        hi = - (nlo + 2)
        xm[-nhi:] = (xwc[N] - xwc[lo:hi]) / (wc[N] - wc[lo:hi])

    return xm

def running_mean_equidistant_with_weights_multiple_axes(
        x,
        nlo,
        nhi,
        w,
        periodic=False,
        ):
    """*Helper function.
    
    NOTE: DOES NOT WORK

    Parameters
    ----------
    x : array
        the operation is only performed along the first axis.
    nlo : int >= 0
        number of points to the left to consider. 0 means only the center
        point.
    nhi : int >= 0
        as nlo, but to the right.
    w : array, same shape as x
        the weight of each value in x to the running mean
    periodic : bool, optional
        Consider the edges of x the be connected or not. Default: False 

    Returns
    -------
    xmean : array, same shape as x
        if not periodic, a lower number of points is taken into account on the
        edges (nlo is smaller on the lower edge, nhi get smaller on the top
        edge).
    """
    raise NotImplementedError('function is buggy.')

    N = len(x)
    S = np.shape(x)
    n = nlo + nhi + 1
    assert n <= N

    Szeros = (1,) + S[1:]
    zeros = np.zeros(Szeros)

    if not periodic:
        xm = np.empty(S)

        # weights:
        wc = np.cumsum(np.insert(w, 0, 0), 0, dtype=float)

        # weighted x:
        xw = x * w 

        xwc = np.cumsum(np.insert(xw, 0, zeros, 0), 0, dtype=float)

        # lower edge:
        xm[:nlo] = xwc[n-nlo : n] / wc[n-nlo : n] 
        # center:
        xm[nlo : -nhi] = (xwc[n:] - xwc[:-n]) / (wc[n:] - wc[:-n])
        # upper edge:
        xm[-nhi:] = (xwc[N] - xwc[-n-1 : -nlo-2]) / (wc[N] - wc[-n-1 : -nlo-2])
        return xm

    else:
        # weights:
        wp = 1. * np.concatenate((zeros, w[-nlo:], w, w[:nhi]))
        wc = np.cumsum(wp, 0)

        # weighted x:
        xw = x * w 
        xwp = 1. * np.concatenate((zeros, xw[-nlo:], xw, xw[:nhi]))
        xwc = np.cumsum(xwp, 0)

        xm = (xwc[n:] - xwc[:-n]) / (wc[n:] - wc[:-n])
        return xm
    
if __name__ == '__main__':
    shape = (3, 100)
    axis = 1
    nlo = 5
    nhi = 5
    periodic = False
    ideal = 10 * np.sin(np.arange(shape[1])/10)
    noise = np.random.random(shape)
    x = (ideal + noise)
    f = running_mean_equidistant_no_weights
    x_mean = f(x, nlo, nhi, periodic, axis)
    print(x)
    print('')
    print(x_mean)

    import matplotlib.pyplot as plt
    plt.close()
    plt.plot(x[1], 'k.')
    plt.plot(x_mean[1], 'r-')
    plt.show()
