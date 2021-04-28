###################################################
# QUANTILES                                       #
###################################################
def quantile_from_histogram(c, f, q, axis=0, mode='linear'):
    """Return a float or an array.

    Parameters
    ----------
    c : array of floats
        classes
    f : array of floats, same shape as c
        frequencies
    q : float
        quantile (between 0 and 1)
    axis : int, optional
        Default: 0
    mode : {'linear', 'lower', 'higher'}, optional
        'linear' interpolates

    Returns
    -------
    class : float or array of such
    """

    cs = np.cumsum(f, axis)
    cs = np.rollaxis(cs, axis, 0)

    warnings.simplefilter('ignore', RuntimeWarning)
    s = 1. * cs / np.sum(f, axis)
    warnings.resetwarnings()

    s = np.rollaxis(s, 0, axis + 1)


    X = threshold_position(pos=c, val=s, thr=q, axis=axis)
    clo, cin, chi, ilo, ihi = X
    if mode == 'linear':
        return cin
    elif mode == 'lower':
        return clo
    elif mode == 'higher':
        return chi

def threshold_position(pos, val, thr, axis=0, allow_borders=True):
    """Return position where threshold is first crossed.

    Parameters
    ----------
    pos : array
        must be sorted in ascending order
    val : array, same shape as pos
    thr : float
    axis : int, optional
        Default: 0
    allow_borders : bool, optional
        see `Special cases` for description. Default: True

    Returns
    -------
    plo : array or float
        last pos below where thr is crossed
    pin : array or float
        interpolated pos where thr is crossed
    phi : array or float
        first pos above thr is crossed
    ilo : int or array of such
        indices of plo in pos
    ihi : int or array of such
        indices of phi in pos

    Special cases
    -------------
    1) if thr is never crossed:
       ========================
        if allow_borders:
            plo == pos[-1]  (at the respective axis)
            pin == pos[-1]  (at the respective axis)
            phi == pos[-1]  (at the respective axis)
            ilo == highest  (at the respective axis)
            ihi == highest  (at the respective axis)
        if not allow_borders:
            plo == pos[-1]  (at the respective axis)
            pin == np.nan
            phi == np.nan
            ilo == highest  (at the respective axis)
            ihi == highest  (at the respective axis)

    2) if the first pos is already above thr:
       ======================================
        if allow_borders:
            plo == pos[0]  (at the respective axis)
            pin == pos[0]  (at the respective axis)
            phi == pos[0]  (at the respective axis)
            ilo == 0       (at the respective axis)
            ihi == 0       (at the respective axis)
        if not allow_borders:
            plo == np.nan
            pin == np.nan
            phi == pos[0]  (at the respective axis)
            ilo == 0       (at the respective axis)
            ihi == 0       (at the respective axis)
    """
    # Note
    # ====
    # In the context of this function, 'flat' means an array that has the same
    # dimensions as pos reduced by the vertical dimension. This usage the term
    # is different than the numpy usage. 
    #
    # `above` is the 'flat' array which contains the index of the layer where
    # val is first above thr. `below` is the 'flat' array with indices just
    # below `above`
    ###################################################
    # SHAPES                                          #
    ###################################################
    S = np.shape(pos)
    Ndims = len(S)
    S_flat = tuple(list(S)[:axis] + list(S)[axis+1:])
    N = S[axis]

    ###################################################
    # FIND LEVEL WHERE THRESHOLD IS CROSSED           #
    ###################################################
    # full size array of booleans:
    smaller = val < thr
    warnings.simplefilter('ignore', RuntimeWarning)
    above = np.sum(np.cumprod(smaller, axis), axis)    # flat array of indices
    warnings.resetwarnings()
    below = above - 1                                  # flat array of indices

    # If val >= thr on the bottom-most level, then `below` == -1.. This is
    # corrected here. The following function yields 0 if below < 0, otherwise
    # `below`. In this formulation, in works both for scalars and arrays.
    below += below < 0

    # likewise is thr is never crossed:
    above -= above >= N
        
    # special cases (treated later):
    always_smaller = (above == N)

    ###################################################
    # CREATE 3-TUPLES OF INDICES FOR INTERPOLATION    #
    ###################################################
    above_idx = np.expand_dims(above, 0)
    below_idx = np.expand_dims(below, 0)
    other_idx = np.indices(above.shape)

    ihi = tuple(np.concatenate(
        (other_idx[:axis], above_idx, other_idx[axis:]), 0))
    ilo = tuple(np.concatenate(
        (other_idx[:axis], below_idx, other_idx[axis:]), 0))

    too_hi = ilo[axis] == N - 1
    too_lo = ihi[axis] == 0


    ###################################################
    # INTERPOLATE                                     #
    ###################################################
    phi = pos[ihi]
    plo = pos[ilo]
    vhi = val[ihi]

    vlo = val[ilo]
    vdiff = vhi - vlo

    # surpress nan-warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        a = (vhi - thr) / vdiff
        b = (thr - vlo) / vdiff

    # special cases
    if Ndims == 1:
        if too_hi:
            a = 1
            b = 0
        elif too_lo:
            a = 0
            b = 1
    elif Ndims > 1:
        a[too_hi] = 1
        b[too_hi] = 0
        a[too_lo] = 0
        b[too_lo] = 1

    pin = a * plo + b * phi        # p interpolated

    ###################################################
    # EXCLUDE INVALID VALUES                          #
    ###################################################
    if not allow_borders:
        if Ndims == 1:  # case: return value is scalar
            if too_hi:
                phi = np.nan
                pin = np.nan
            if too_lo:
                plo = np.nan
                pin = np.nan

        else:
            phi[too_hi] = np.nan
            plo[too_lo] = np.nan

            pin[too_hi] = np.nan
            pin[too_lo] = np.nan

    return plo, pin, phi, ilo, ihi


