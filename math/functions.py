# PyPI
import numpy as np

def gauss(fwhm=None, center=0., amplitude=1., sigma=None):
    """Return a Gauss function.
    
        Parameters
        ----------
        fwhm : float or None
            full width at half maximum
        sigma : float or None
            mutually exclusive with fwhm

        center : float
            (default: 0.)
        amplitude : float
            integral over the function (default: 1.)

        Returns
        -------
        callable
    """
    if sigma is None and fwhm is None:
        raise ValueError('Either fwhm or sigma must be given, but not both.')

    if sigma is not None and fwhm is not None:
        raise ValueError('Either fwhm or sigma must be given, but not both.')

    if sigma is None:
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))

    A  = amplitude / (sigma * np.sqrt(2*np.pi)) # normalization factor
    C  = -1 / (2 * sigma**2)                    # exponent factor
    x0 = center
    f = lambda x: A * np.exp(C * (x-x0)**2)
    return f
    
def gauss2d(xwidth, ywidth, xcenter=0., ycenter=0., amplitude=1.):
    """Return a 2 dimensional Gauss function.
    
    Parameters:
    * xcenter (default: 0.)
    * ycenter (default: 0.)
    * amplitude: the integral over the function (default: 1.)"""
    A  = float(amplitude)/(xwidth * ywidth * 2*np.pi) # amplitude factor
    Cx  = -1/(2*float(xwidth)**2) # x-exponent factor
    Cy  = -1/(2*float(ywidth)**2) # y-exponent factor
    x0 = float(xcenter)
    y0 = float(ycenter)
    def f(x,y):
        return A * np.exp(Cx*(x-x0)**2 + Cy*(y-y0)**2)
    return f    
    
def triangle(width=1., center=0., amplitude=1.):
    """Return a 1D triangle function.
    
    Parameters:
    * center (default: 0.)
    * amplitude: the integral over the function (default: 1.)
    
    Returns the function
    
    f(x) = amplitude/width * (1 - |x-center|/width) if |x-center| < width
         = 0.                                          else
    """
    A  = float(amplitude)
    w  = float(width)
    x0 = float(center)
    """
    For speeding up the function, it is evaluated in a different, but
    mathematically equivalent way:
    f(x) = A/w - A*x0/w**2 + A/w**2 * x   if x < x0
    f(x) = A/w + A*x0/w**2 - A/w**2 * x   if x > x0 
            |       |           |     |
            V       V           V     V
            c  -/+  d     +/-   m   * x
    """
    c  = A/w
    d  = A*x0/w**2
    m  = A/w**2
    elarge = c + d
    esmall = c - d
    def f(x):
        # trivial cases (x out of range):
        if x-x0 >= width:
            return 0.
        elif x0-x >= width:
            return 0.
        # non-trivial cases:
        elif x < x0:
            return esmall + m * x
        else:
            return elarge - m * x
    return f

def cone(xwidth=1., ywidth=1., xcenter=0., ycenter=0., volume=None,
        height=None):
    """Return a 2 dimensional elliptic cone function.
    
    The returned function z(x,y) gives the height of the cone at the point
    (x,y) if they are within the cone and 0 otherwise.
    
    Parameters
    ----------
    * xwidth (default: 1.) width of the cone base in x-direction
    * ywidth (default: 1.) width of the cone base in y-direction
    * xcenter (default: 0.)
    * ycenter (default: 0.)
    * you can either specify volume OR height. If don't specify any of them,
    the volume will be set to 1.
    
    xwidth and ywidth are allowed to be negative. One of them being negative
    and one positve results in negative function values.
    
    Negative height or volume will alter the sign of the function values, too.
    
    Written in 2015
    by Andreas Anhaeuser
    Insitute for Geophysics and Meteorology
    University of Cologne
    Germany
    <andreasfrederik.anhaeuser@smail.uni-koeln.de> """
    
    x0 = float(xcenter)
    y0 = float(ycenter)
    wx = float(xwidth)
    wy = float(ywidth)

    if wx == 0 or wy == 0:
        raise ValueError('xwidth and ywidth must not be 0.')
    # Volume/height:
    if not height or np.isnan(height):
        if not volume or np.isnan(volume):
            V = 1.
        else:
            V = float(volume)
        h  = 3 * V / (np.pi * wx * wy)
    elif not volume or np.isnan(volume):
        h = float(height) * np.sign(wx) * np.sign(wy)
    else:
        raise Exception('Not allowed to specify both volume AND height.')

    def f(x, y):
        dist = np.sqrt(((x - x0) / wx)**2 + ((y - y0) / wy)**2)
        if dist > 1:
            return 0.
        else:
            return h * (1. - dist)
    return f    
