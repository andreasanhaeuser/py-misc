#####################################
# SMOOTHING                         #
#####################################
def smooth_distribution(
        data,
        kernel_function,
        grid,
        window_halfwidth=np.inf,
        normalize=False,
        ):
    """Smooth a discrete occurence density function.
    
    INPUT:
    * data (list or array)
    * kernel_function (function handle):
        the function that the ODF should be convolved with
    * grid (list or array): output sample grid
    
    Parameter
    * normalize (boolean): if True, then the integral over the output 
        function is 1.
        otherwise the integral is roughly (but not exactly) the number of 
        elements in occurence_values.
    
    OUTPUT:
    an array that is the smoothed ODF/PDF over grid.
    """
    f  = kernel_function
    # initialize:
    y = np.array([0. for x in grid])
    # convolvle:
    for d in data:
        y += np.array([f(d-float(x)) for x in grid])
    # normalize:
    if normalize:
        y = y/sum(y)        
    return y


#####################################
# CONVOLUTIONS                      #
#####################################
def convolve_binned_distribution_1d_equidistant(
        data,
        kernel,
        periodic=False,
        contains_nans=True,
        ):
    """Return the convolved 1D array..
    
    Convolvse (smooth) the data given in z by a 1D-kernel function.
    
    INPUT
    * data: a 2D array
    * kernel: a 2D array. data will be colvolved with this kernel
    * periodic: (boolean) denoting whether the bounds are periodic or not.
    * contains_nans: (boolean) setting this to True will speed up the
      function by about a factor 2.
      NOTE: Setting this to False while data or kernel actually contain NaN's
      may cause false results!
      
    OUTPUT:
    A 1D array: the values of this convolution.
    
    Written in 2015
    by Andreas Anhaeuser
    Insitute for Geophysics and Meteorology
    University of Cologne
    Germany
    <andreasfrederik.anhaeuser@smail.uni-koeln.de> """ 
    Ik = len(kernel)
    I  = len(data)
    c = np.zeros(I)
        
    def index_range_function(N, Nk, periodic=False):
        if not periodic:
            def func(n):
                nlo = max(0, n - (Nk//2))
                nhi = min(N, n + (Nk - 1) // 2 + 1)
                return range(nlo, nhi)
        else:
            def func(n):
                nlo = (n - (Nk//2)) % N
                nhi = (n + (Nk - 1) // 2 + 1) % N
                if nlo >= nhi:
                    nlo -= N
                return range(nlo, nhi)
        return func
        
    def kernel_index_range(n, nn, N, Nk):
        shift = n - Nk//2
        first = (nn[0] - shift) % N
        L = len(nn)
        return range(first, first + L)

    get_ii = index_range_function(I, Ik, periodic)
    
    if periodic and not contains_nans:
        C = 1./np.sum(kernel)    
        for i in range(I):
            ii = get_ii(i)
            d = data[ii]
            c[i] = C * np.sum(kernel * d)
    elif not contains_nans:
        for i in range(I):
            ii = get_ii(i)
            ik = kernel_index_range(i, ii, I, Ik)
            d = data[ii]
            k = kernel[ik]
            C = 1./np.sum(k)
            c[i] = C * np.sum(k * d)  
    else:
        for i in range(I):
            ii = get_ii(i)
            ik = kernel_index_range(i, ii, I, Ik)
            d = data[ii]
            k = kernel[ik]
            try:
                C = 1./np.sum(k[-np.isnan(d)])
            except ZeroDivisionError:
                C = np.nan
            d[np.isnan(d)] = 0.
            c[i] = C * np.sum(k * d)  
    return c
    
def convolve_binned_distribution_1d(
        x,
        z,
        kernel_function,
        window_halfwidth=np.inf,
        assume_equidistant=False,
        periodicity=np.inf,
        assume_sorted=False,
        contains_nans=True,
        ):
    """Return a 2D array.
    
    Convolvse (smooth) the data given in z by a 2D-kernel function.
    
    INPUT
    * 
      
    OUTPUT:
    (x_sorted, y_sorted, z_conv)
    
    Written in 2015
    by Andreas Anhaeuser
    Insitute for Geophysics and Meteorology
    University of Cologne
    Germany
    <andreasfrederik.anhaeuser@smail.uni-koeln.de> """ 
    
    if not assume_sorted:
        idx = np.argsort(x)
        inv_idx = np.argsort(idx)
        return convolve_binned_distribution_1d(
            x[idx],
            z[idx,:],
            kernel_function,
            window_halfwidth=window_halfwidth,
            assume_equidistant=assume_equidistant,
            periodicity=periodicity,
            assume_sorted=True,
            contains_nans=contains_nans,
            )
        
    if assume_equidistant:
        per = 0 < periodicity < np.inf
        #### find kernel x grid:
        hw = window_halfwidth
        delta = x[1] - x[0]
        # half number of elements in x direction:
        hn = hw//delta
        xk = np.arange(-hn, hn+1) * hw/hn
        
        kernel = kernel_function(xk)
        zc = convolve_binned_distribution_1d_equidistant(
            z,
            kernel,
            periodic=per,
            contains_nans=contains_nans,
            )
        return (x, zc)
        
    else:
        raise NotImplementedError('')
    
def convolve_distribution_2d(
        data,
        xgrid,
        ygrid,
        kernel_function,
        x_periodicity=np.inf,
        y_periodicity=np.inf,
        x_window_halfwidth=np.inf,
        y_window_halfwidth=np.inf,
        performance_info_interval=30,
        ):
    """Return a 2D array corresponding to xgrid and ygrid.
    
    Convolve a discrete distribution of data points with a (smoothing)
    2D-kernel function and return the values of this convolution on the xgrid
    and ygrid points as a 2D array.
    
    INPUT
    * data: a list of pairs (lists or tuple of length 2) that specifies x and
      y coordinates of the data points.
    * xgrid, ygrid: list or 1D arrays of the grid coordinates on which the
      convolution will be calculated.
    * kernel_function: the convolution kernel f(x, y)
    * x_periodicity: (positive value) if this is finite, the 'edges' in x
      direction are 'glued' together with a 'cylinder circumference' of
      x_periodicity, all coordinates (x+N*x_periodicity) are equivalently
      treated.
    * y_periodicity: analogously to x_periodicity
    * x_window_halfwidth: a positive value not larger than x_periodicity/2.
      To save computation time, coordinate differences larger than this value
      are not considered (set this value reasonably larger than your kernel
      width)
    * y_window_halfwidth: as x_window_halfwidth
    * performance_info_interval (seconds): the frequency with which
      performance infromation is printed on stdout.
      
    OUTPUT:
    A 2D array: the values of this convolution on the xgrid and ygrid points.
    """
    
    start = dt.datetime.now()
    
    X = len(xgrid)
    Y = len(ygrid)
    f = kernel_function
    c = np.zeros([X, Y])
    
    xwin = x_window_halfwidth
    ywin = y_window_halfwidth
    
    N = len(data)
    n = 0
    last_message_time = start
    tolerance = dt.timedelta(seconds=performance_info_interval)
    for d in data:
        n += 1
        time = dt.datetime.now()
        if time - last_message_time > tolerance:
            last_message_time = time
            print('loop ' + str(n) + '/' + str(N) + 
                ', time elapsed = ' + str(time-start))
            
        xd, yd = d
        # initialize arrays of indices of the coordinates within the window:
        iarr = []
        jarr = []
        # ... and the corresponding coordiante differences:
        dx   = []
        dy   = []
        # search indices and coordinate differences:
        for i in range(X):
            delta_x = min(
                    np.mod(xgrid[i]-xd, x_periodicity),
                    np.mod(xd-xgrid[i], x_periodicity)
                    )
            if delta_x < x_window_halfwidth:
                iarr.append(i)
                dx.append(delta_x)
                
        for j in range(Y):
            delta_y = min(
                    np.mod(ygrid[j]-yd, y_periodicity),
                    np.mod(yd-ygrid[j], y_periodicity)
                    )
            if delta_y < y_window_halfwidth:
                jarr.append(j)
                dy.append(delta_y) 
        # calculate the convulution:
        for ii in range(len(iarr)):
            for jj in range(len(jarr)):
                i = iarr[ii]
                j = jarr[jj]
                c[i,j] += f(dx[ii], dy[jj])
    return c

def convolve_binned_distribution_2d_equidistant(
        data,
        kernel,
        periodic_x=False,
        periodic_y=False,
        contains_nans=True,
        ):
    """Return the convolved 2D array..
    
    Convolvse (smooth) the data given in z by a 2D-kernel function.
    
    INPUT
    * data: a 2D array
    * kernel: a 2D array. data will be colvolved with this kernel
    * periodic_x, periodic_y: (boolean) denoting whether the bounds in x- or
      y-direction are periodic or not.
    * contains_nans: (boolean) setting this to True will speed up the
      function by about a factor 2.
      NOTE: Setting this to False while data or kernel actually contain NaN's
      may cause false results!
      
    OUTPUT:
    A 2D array: the values of this convolution.
    
    Written in 2015
    by Andreas Anhaeuser
    Insitute for Geophysics and Meteorology
    University of Cologne
    Germany
    <andreasfrederik.anhaeuser@smail.uni-koeln.de> """ 
    Ik, Jk = np.shape(kernel)
    I, J = np.shape(data)
    c = np.zeros([I, J])
        
    def index_range_function(N, Nk, periodic=False):
        if not periodic:
            def func(n):
                nlo = max(0, n - (Nk//2))
                nhi = min(N, n + (Nk - 1) // 2 + 1)
                return range(nlo, nhi)
        else:
            def func(n):
                nlo = (n - (Nk//2)) % N
                nhi = (n + (Nk - 1) // 2 + 1) % N
                if nlo >= nhi:
                    nlo -= N
                return range(nlo, nhi)
        return func
        
    def kernel_index_range(n, nn, N, Nk):
        shift = n - Nk//2
        first = (nn[0] - shift) % N
        L = len(nn)
        return range(first, first + L)

    get_ii = index_range_function(I, Ik, periodic_x)
    get_jj = index_range_function(J, Jk, periodic_y)
    
    if periodic_x and periodic_y and not contains_nans:
        C = 1./np.sum(kernel)    
        for j in range(J):
            jj = get_jj(j)
            for i in range(I):
                ii = get_ii(i)
                d = data[ii,:][:,jj]
                c[i,j] = C * np.sum(kernel * d)
    elif not contains_nans:
        for j in range(J):
            jj = get_jj(j)
            jk = kernel_index_range(j, jj, J, Jk)
            for i in range(I):
                ii = get_ii(i)
                ik = kernel_index_range(i, ii, I, Ik)
                d = data[ii,:][:,jj]
                k = kernel[ik,:][:,jk]
                C = 1./np.sum(k)
                c[i,j] = C * np.sum(k * d)  
    else:
        for j in range(J):
            jj = get_jj(j)
            jk = kernel_index_range(j, jj, J, Jk)
            for i in range(I):
                ii = get_ii(i)
                ik = kernel_index_range(i, ii, I, Ik)
                d = data[ii,:][:,jj]
                k = kernel[ik,:][:,jk]
                try:
                    C = 1./np.sum(k[-np.isnan(d)])
                except ZeroDivisionError:
                    C = np.nan
                d[np.isnan(d)] = 0.
                c[i,j] = C * np.sum(k * d)  
    return c

def convolve_binned_distribution_2d(
        x,
        y,
        z,
        kernel_function,
        x_window_halfwidth=np.inf,
        y_window_halfwidth=np.inf,
        assume_equidistant_x=False,
        assume_equidistant_y=False,
        periodicity_x=np.inf,
        periodicity_y=np.inf,
        assume_sorted_x=False,
        assume_sorted_y=False,
        contains_nans=True,
        ):
    """Return a 2D array.
    
    Convolvse (smooth) the data given in z by a 2D-kernel function.
    
    INPUT
    * 
      
    OUTPUT:
    (x_sorted, y_sorted, z_conv)
    
    Written in 2015
    by Andreas Anhaeuser
    Insitute for Geophysics and Meteorology
    University of Cologne
    Germany
    <andreasfrederik.anhaeuser@smail.uni-koeln.de> """ 
    
    if not assume_sorted_x:
        idx = np.argsort(x)
        inv_idx = np.argsort(idx)
        return convolve_binned_distribution_2d(
            x[idx],
            y,
            z[idx,:],
            kernel_function,
            x_window_halfwidth=x_window_halfwidth,
            y_window_halfwidth=y_window_halfwidth,
            assume_equidistant_x=assume_equidistant_x,
            assume_equidistant_y=assume_equidistant_y,
            periodicity_x=periodicity_x,
            periodicity_y=periodicity_y,
            assume_sorted_x=True,
            assume_sorted_y=assume_sorted_y,
            contains_nans=contains_nans,
            )
    if not assume_sorted_y:
        idx = np.argsort(y)
        inv_idx = np.argsort(idx)
        return convolve_binned_distribution_2d(
            x,
            y[idx],
            z[:,idx],
            kernel_function,
            x_window_halfwidth=x_window_halfwidth,
            y_window_halfwidth=y_window_halfwidth,
            assume_equidistant_x=assume_equidistant_x,
            assume_equidistant_y=assume_equidistant_y,
            periodicity_x=periodicity_x,
            periodicity_y=periodicity_y,
            assume_sorted_x=True,
            assume_sorted_y=True,
            contains_nans=contains_nans,
            )
        
    if assume_equidistant_x and assume_equidistant_y:
        per_x = 0 < periodicity_x < np.inf
        per_y = 0 < periodicity_y < np.inf
        #### find kernel x grid:
        xhw = x_window_halfwidth
        xdelta = x[1] - x[0]
        # half number of elements in x direction:
        hn = xhw//xdelta
        xk = np.arange(-hn, hn+1) * xhw/hn
        #### find kernel y grid:
        yhw = y_window_halfwidth
        ydelta = y[1] - y[0]
        # half number of elements in x direction:
        hn = yhw//ydelta
        yk = np.arange(-hn, hn+1) * yhw/hn
        
        mky, mkx = np.meshgrid(yk, xk)
        kernel = kernel_function(mkx, mky)
        zc = convolve_binned_distribution_2d_equidistant(
            z,
            kernel,
            periodic_x=per_x,
            periodic_y=per_y,
            contains_nans=contains_nans,
            )
        return (x, y, zc)
        
    else:
        raise NotImplementedError('')
