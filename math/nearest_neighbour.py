# standard modules
import datetime as dt

# PyPI modules
import numpy as np

# IGMK modules
import misc.date_time.utils as dt_utils

def nearest_neighbour(
        x_array, x, direction='both', values=None, skip_nans=False,
        assume_sorted=False
        ):
    """Return the nearest neighbour in the given array.
    
        Parameters
        ----------
        x_array : array
        x : float
        direction : str
            'both', 'up' or 'down'. specifies in which direction(s) to look for
            the nearest neighbour
        values : array
            a 1d array or list of the same length as x_array if values is
            given, then the function also returns the value at the nearest
            neighbour position
        skip_nans : bool
            if True, then positions in x_array are ignored if the corresponding
            element in values is nan.
        assume_sorted : bool
            if True, complexity is O(log(n)), else O(n)
      
        Returns
        -------
        (float, int, float) or
        (float, int)
            (x_nearest, index, value)
            x_nearest :
                is the element in x_array closest to x, respecting the
                restrictions given by direction an skip_nans.
            index : 
                the position of x_nearest in a_array
            value : 
                the element of values corresponding to x_nearest.
        
        Author
        ------
        Written in 2014
        by Andreas Anhaeuser
        Insitute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """    
    ###################################################
    # ABBREVIATIONS                                   #
    ###################################################
    direction = direction[:1].lower()
    xa = x_array

    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    if values is not None:
        assert len(values) == len(xa)
    
    assert direction in ['b', 'd', 'u']
    
    #######################################################
    # EMPTY INPUT ARRAY                                   #
    #######################################################
    if not any(xa):
        if values is None:
            return None, None
        else:
            return None, None, None
    
    #######################################################
    # DATETIME OBJECT HANDLING                            #
    #######################################################
    if isinstance(xa[0], dt.datetime):
        xa_new = dt_utils.datetime_to_seconds(xa)
        out = nearest_neighbour(
                xa_new,
                x,
                direction=direction,
                values=values,
                skip_nans=skip_nans,
                assume_sorted=assume_sorted,
                )
        nn = out[0]
        if nn is not None:
            nndt = dt_utils.seconds_to_datetime(nn)
        else:
            nndt = None
        return (nndt,) + out[1:]

    if isinstance(x, dt.datetime):
        x_new = dt_utils.datetime_to_seconds(x)
        return nearest_neighbour(
                xa,
                x_new,
                direction=direction,
                values=values,
                skip_nans=skip_nans,
                assume_sorted=assume_sorted,
                )

    ###################################################
    # CAST TO ARRAY                                   #
    ###################################################
    xa = np.array(xa)
    if values is not None:
        values = np.array(values)


    ###################################################
    # SORT                                            #
    ###################################################
    # if array unsorted then recursive function call:
    if not assume_sorted:
        idx = np.argsort(xa)
        xa = xa[idx]
        if values is not None:
            values = values[idx]
        
        # recursively call function:
        result = nearest_neighbour(
                    x_array=xa,
                    x=x,
                    direction=direction,
                    values=values,
                    skip_nans=skip_nans,
                    assume_sorted=True
                    )
        
        xn = result[0]
        i  = result[1]
        
        # output:
        if i is None and values is None:
            return None, None
        
        elif i is None:
            return None, None, None
        
        elif values is None:
            return (xn, idx[i])
        
        else:
            return (xn, idx[i], result[2])
            
    ###################################################
    # CALL FAST FUNCTION                              #
    ###################################################
    result = nearest_neighbour_sorted(x=xa, x0=x, direction=direction)
    
    ###################################################
    # SELECT ilo OR ihi                               #
    ###################################################
    if direction in 'du':
        i = result

    if direction == 'b':
        ilo, ihi = result
        if ilo is None:
            i = ihi
        elif ihi is None:
            i = ilo
        else:
            xlo = xa[ilo]
            xhi = xa[ihi]
            dlo = x - xlo
            dhi = xhi - x
            if dlo <= dhi:
                i = ilo
            else:
                i = ihi

    ###################################################
    # RETURN                                          #
    ###################################################
    # regular case:
    if i is not None:
        if values is not None:
            return xa[i], i, values[i]
        else:
            return xa[i], i

    # no nearest neighbour:
    else:
        if values is not None:
            return None, None, None
        else:
            return None, None

def nearest_neighbour_sorted(x, x0, direction='both'):
    """Return the nearest neighbour in the given array.

        Helper function to nearest_neighbour
    
        Parameters
        ----------
        x : array
        x0 : float
            the value to find
        direction : str
            'both', 'up' or 'down'. specifies in which direction(s) to look for
            the nearest neighbour
      
        Returns
        -------
        index : int or None or pair of such
            if direction is 'up' or 'down', int or None is returned
            if direction is 'both', a pair of such is returned
        
        Author
        ------
        Written 2014-2016
        by Andreas Anhaeuser
        Insitute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """    
    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    direction = direction[:1].lower()
    assert direction in ['b', 'd', 'u']
        
    #######################################################
    # EMPTY INPUT ARRAY                                   #
    #######################################################
    if not any(x):
        return None

    # initialize
    N = len(x)
    n = N // 2
    xn = x[n]
    dn = n

    ###################################################
    # FIND NEAREST INDEX                              #
    ###################################################
    while dn > 1:
        dn = (dn + 1) // 2
        if x0 < xn:
            n = n - dn
        else:
            n = n + dn
        if n >= N:
            n = N - 1
        xn = x[n]

    if x0 < xn:
        n = n - 1
        xn = x[n]

    ###################################################
    # nlo, nhi                                        #
    ###################################################
    # Note: n == nlo

    # x0 < x[0]
    if n == 0 and  x0 < xn:
        n = None

    # nhi
    if direction in 'ub':
        # x0 < x[0]
        if n is None:
            nhi = 0

        # x[-1] < x0
        elif n == N - 1:
            nhi = None

        # xlo == xhi
        elif xn == x0:
            nhi = n

        # regular case
        else:
            nhi = n + 1

    ###################################################
    # RETURN                                          #
    ###################################################
    if direction == 'u':
        return nhi
    elif direction == 'd':
        return n 
    else:
        return n, nhi
