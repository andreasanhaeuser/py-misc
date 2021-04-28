#! /usr/bin/env python
"""A collection of functions related to meteorology.

    Units
    -----
    Preferred units are SI base units:

    use      | avoid
    ---------+---------------
    K        | deg C
    Pa       | mbar, bar, hPa
    kg       | g
    kg/kg    | g/kg
    kg/m^3   | g/m^3
    kg/mol   | g/mol
    fraction | %


    Variable names
    --------------
    t : time (s)
    z : elevation (m)
    p : total pressure (Pa)
    e : partial pressure (Pa)
     
    u : zonal wind speed (m/s)
    v : meridional wind speed (m/s)
    w : vertical wind speed (m/s)
    U : horizontal wind speed (m/s)  [== sqrt(u**2 + v**2)]
    V : total wind speed (m/s)  [== sqrt(u**2 + v**2 + w**2)]

    T : temperature (K)
    Td : dew point temperature (K)
    th : potential temperature (K)
    thv : virtual potential temperature (K)
    L : lapse rate (K/m)

    q : volumetric humidity (kg/m3)
    r : relative humidity  [1: saturated]
    s : specific humidity
    m : mixing ratio


    Arrays
    ------
    Most functions which take scalars also work with array albeit not
    explicitly stated. Just try!


    Author
    ------
    Written in 2015-2016
    by Andreas Anhaeuser
    Insitute for Geophysics and Meteorology
    University of Cologne
    Germany
    <anhaeus@meteo.uni-koeln.de>
"""

# standard modules
from bisect import bisect
import collections
from copy import deepcopy as copy
import warnings

# PyPI modules
import numpy as np

#######################################
# CONSTANTS                           #
#######################################
# physical constants:
_R = 8.3144598  # (J K-1 mol-1) universal gas constant

# geophysical constants:
_g = 9.80655  # (m s-2) gravitational acceleration

# molecular masses (kg mol-1)
_Md = 28.9645e-3   # dry air
_Mw = 18.01528e-3  # water
_epsilon = _Mw / _Md

# specific gas constants (J K-1 kg-1):
_Rd = _R / _Md   # dry air
_Rw = _R / _Mw   # water 

# isobaric heat capacity (J K-1 kg-1):
_cp_d = 1.005e3  # dry air
_cp_w = 1.860e3  # water vapor

# isochoric heat capacity (J K-1 kg-1):
_cv_d = 0.718e3  # dry air

# latent heat, specific (J kg-1):
_L = 2264.76e3   # water, liquid-vapor
_L_f = 334e3     # water, solid-liquid

#######################################
# SATURATION PRESSURE                 #
#######################################
def water_vapor_saturation_pressure(T, phase='liquid'):
    """Return water vapor saturation pressure.

        Uses Arden Buck equation. See
        - Buck (1996), Buck Research CR-1A User's Manual, Appendix 1
        - https://en.wikipedia.org/wiki/Arden_Buck_equation

        Parameters
        ----------
        T : float or array
            temperature (K)
        phase : str, optional
            'liquid' or 'solid'. Default is 'liquid'

        Returns
        -------
        float or array
            water vapor saturation pressure (Pa)

        Accuracy 
        --------
        unknown
    """
    ph = phase[0].lower()
    # liquid
    if ph in ['l', 'w']:
        a = 6.1121
        b = 18.678
        c = 257.14
        d = 234.5
    # solid
    elif ph in ['s', 'i']:
        a = 6.1115
        b = 23.036
        c = 279.82
        d = 333.7
    else:
        raise ValueError("phase must be 'liquid' or 'solid'.")
    
    TC = T - 273.15  # temperature in deg C

    # Arden Buck Formula:
    es_hPa = a * np.exp((b - TC/d) * (TC / (c + TC)))

    # convert from hPa to Pa:
    return es_hPa * 1e2

def r_from_Td(T, Td):
    """Return relative humidity from dew point temperature.

        Parameters
        ----------
        T : float or array of such
            (K) temperature
        Td : float or array of such
            (K) dew point temperature

        Returns
        -------
        r : float or array of such
            relative humidity (1: saturated)

        Author
        ------
        Andreas Anhaeuser (AA) <anhaeus@meteo.uni-koeln.de>
        Institute for Geophysics and Meteorology
        University of Cologne, Germany

        History
        -------
        2016-09-30 (AA): Created
    """
    es_Td = water_vapor_saturation_pressure(Td)
    es_T  = water_vapor_saturation_pressure(T)
    return es_Td / es_T

#######################################
# PRESSURE                            #
#######################################
def standard_pressure(z, z0=0., T0=300., p0=1013e2, L=6.5e-3,
        z_tropopause=11e3, z_stratosphere1=20e3, T_tropopause=217.):
    """Return a standard atmosphere pressure array.

        Parameters
        ----------
        z : array
            (m) geopotential height
        z0 : float
            (m) reference geopotential height
        T0 : float
            (K) temperature at z0
        p0 : float
             (Pa) pressure at z0
        L : array, optional
            (K/m) temperature lapse rate -dT/dz
            default is 6.5e-3

        Returns
        -------
        array
            pressure (Pa) corresponding to the height levels in z, using the
            barometric formula.
    """
    if np.nanmax(z) > z_tropopause:
        raise NotImplementedError(
                'Only troposphere implemented yet.')
    p = p0 * (1 - L * (z - z0) / T0)**(_g * _Md / (_R * L))
    return p

def equal_mass_range(zmin, zmax, N, z0=0., T0=300., p0=1013e2, L=6.5e-3):
    """Return levels which devide the atmosphere into layers of equal mass.

        Parameters
        ----------
        zmin : float
            (m) lowest level of the return array
        zmax : float
            (m) highest level of the return array (inclusive)
        N : int
            length of the return array
        z0 : float, optional
            (m) reference level (for T0 and p0), default is 0.
        T0 : float, optional
            (K) temperature at z0, default is 300.
        p0 : float, optional
            (Pa) pressure at z0, default is 1013e2
        L : float, optional
            (K/m) lapse rate, default is 6.5e-3

        Returns
        -------
        array
            (m) height levels
    """
    ###################################
    # FIND START AND END PRESSURE     #
    ###################################
    
    p_beg = standard_pressure(zmin, z0, T0, p0, L) 
    p_end = standard_pressure(zmax, z0, T0, p0, L) 

    # create pressure array:
    p_inc = (p_end - p_beg) / (N - 1.)
    p = np.arange(p_beg, p_end + p_inc / 2, p_inc)
    """
    The division by 2 in p_end + p_inc/2 has a numercial reason: In some cases,
    rounding errors would lead to an array of length (N + 1) without this
    adjustment.
    """

    # invert standard atmosphere pressure formula:
    alpha =  (_R * L) / (_g * _Md)
    return z0 + T0 / L * (1. - (p/p0)**alpha)

def equal_transmission_range(zmin, zmax, N, alpha, z0=0., T0=300, L=6.5e-3,
        res=None):
    """Return levels that devide the atmosphere into equal transmission layers.

        It is assumed that each layer emits the same radiance and that the
        absorption only depends on the mass of air that is traversed.

        Furthermore, the atmosphere is assumed to be in hydrostatic and
        thermodynamic equilibrium with a constant lapse rate within the
        considered range.

        The devision is made such that a ground-based radiometer is equally
        sensitive to radiation from each layer. Since lower layers absorb
        radiation from higher layers, the higher layers must be thicker in
        order to achieve this goal.  The fact that the air density decreases
        with altitude requires an additional increase of the spacing. Both
        effects, the lower air density and the absorption from lower layers are
        taken into account in the algorithm.

        Parameters
        ----------
        zmin : float
            (m) lowest level of the return array
        zmax : float
            (m) highest level of the return array (inclusive)
        N : int
            length of the return array
        alpha0 : float
            (m-1) absorption coefficient at z0.
        z0 : float, optional
            (m) reference level (for T0 and alpha0), default is 0.
        T0 : float, optional
            (K) temperature at z0, default is 300.
        L : float, optional
            (K/m) lapse rate, default is 6.5e-3
        res : float, optional
            (m) internal height resolution. Default: (zmax - zmin) / (1000 * N)

        Returns
        -------
        array
            (m) height levels

        Note
        ----
        The absorption coefficient alpha is defined as:
        I(x) = I(0) * exp(- alpha * x)
        I : radiance,
        x : path length
    """
    ###################################
    # INCREMENTAL RADIANCE            #
    ###################################
    N = int(N)
    deltaz = zmax - zmin
    if res is None:
        res = deltaz / (N * 1000.)
    dz = res 
    z = np.arange(zmin, zmax + dz / 2, dz)
    Nz = len(z)
    # dI(z) is the fraction of the radiance emitted in height z that reaches
    # the ground.  I is (proportinal to) the integral.
    a = alpha * _Rd / _g * T0  # pre-factor
    b = L / T0
    gamma = _g * _Md / (_R * L)
    dI = np.exp(a * ((1 - b * (z - z0))**gamma - 1),
            dtype=np.longdouble)
    I = np.cumsum(dI) * dz

    A = I[-1]
    dA = A / (N - 1)

    zout = np.nan * np.ones(N)
    for n in range(1, N - 1):
        Iexact = n * dA
        i = bisect(I, Iexact)
        if i == Nz:
            i = Nz - 1
        # interpolate:
        DI = I[i] - I[i-1]
        Dz = z[i] - z[i-1]
        zout[n] = (Iexact * Dz + z[i-1] * I[i] - z[i] * I[i-1]) / DI

    zout[0]  = zmin
    zout[-1] = zmax
    return zout

def barometric_altitude(p, z0=0., T0=300., p0=1013e2, L=6.5e-3,
        p_tropopause=200e2,):
    """Return a standard atmosphere altitude array.

        Parameters
        ----------
        p : array
            (Pa) pressure
        z0 : float
            (m) reference geopotential height
        T0 : float
            (K) temperature at z0
        p0 : float
             (Pa) pressure at z0
        L : array, optional
            (K/m) temperature lapse rate -dT/dz
            default is 6.5e-3

        Returns
        -------
        array
            altitude (m) corresponding to the pressure `p`, using the
            barometric formula.
    """
#    assert np.sum(p < p_tropopause) == 0
    z = z0 + T0/L * (1 - (p/p0)**(_R * L / (_g * _Md)))
    return z

#######################################
# TEMPERATURE                         #
#######################################
def potential_temperature(p, T, p0=1e5, mV=0., mL=0., mS=0., mode='normal'):
    """Works for cloudy and cloud-free conditions.

        Parameters
        ----------
        p : float
            (Pa) pressure
        T : float
            (K) temperature
        mV : float
            mixing ratio of water vapor
        mL : float, optional
            mixing ratio of liquid water, default: 0.
        mS : float, optional
            mixing ratio of solid water, default: 0.
        mode : str, optional
            'normal' or 'cosmo'. If 'cosmo' then it is computed to exactly
            match the value in COSMO code (albeit less precise). Default:
            'normal'

        Returns
        -------
        float
            (K) potential temperature

        Note
        ----
        Works also with arrays.
    """
    if mode == 'normal':
        m_tot = mV + mL + mS            # total water mixing ratio (all phases)
        R  = _Rd   * (1 - mV) + _Rw   * mV
        cp = _cp_d * (1 - m_tot)   + _cp_w * mV
        kappa = R / cp
    elif mode == 'cosmo':
        # CONSTANTS AS COSMO USES THEM:
        p0   = 1e5     # (Pa)
        Rd   = 287.05 
        cp_d = 1005.0  
        kappa = Rd / cp_d
    else:
        raise ValueError('Unrecognized mode.')

    th = T * (p0 / p) ** kappa
    return th

def virtual_temperature(p, T, mV, mL=0., mS=0.):
    """Works for cloudy and cloud-free conditions.

        Parameters
        ----------
        p : float
            (Pa) pressure
        T : float
            (K) temperature
        mV : float
            mixing ratio of water vapor
        mL : float, optional
            mixing ratio of liquid water, default: 0.
        mS : float, optional
            mixing ratio of solid water, default: 0.

        Returns
        -------
        float
            (K) virtual temperature

        Note
        ----
        Works also with arrays.
    """
    Tv = T * (1 + mV / _epsilon) / (1 + mV + mL + mS)
    return Tv

def virtual_potential_temperature(p, T, mV, mL=0., mS=0., p0=1e5, mode='normal'):
    """Works for cloudy and cloud-free conditions.

        Parameters
        ----------
        p : float
            (Pa) pressure
        T : float
            (K) temperature
        mV : float
            (kg/kg) mixing ratio of water vapor
        mL : float, optional
            (kg/kg) mixing ratio of liquid water, default: 0
        mS : float, optional
            (kg/kg) mixing ratio of solid water, default: 0
        p0 : float, optional
            (Pa) reference pressure, default: 1e5
        mode : str, optional
            'normal' or 'cosmo'. If 'cosmo' then it is computed to exactly
            match the value in COSMO code (albeit less precise). Default:
            'normal'

        Returns
        -------
        float
            (K) virtual potential temperature
     
        Note
        ----
        Works also with arrays.
    """
    th = potential_temperature(p, T, mV=mV, mL=mL, mS=mS, mode=mode)
    if mode == 'normal':
        thv = virtual_temperature(p, th, mV=mV, mL=mL, mS=mS)
    elif mode == 'cosmo':
        rvdmo = 0.60776868141438767  # cosmo value for (rw/rd - 1)
        thv = th * (1 + rvdmo * mV)

    return thv

# ALIASES:
pt  = potential_temperature
vt  = virtual_temperature
vpt = virtual_potential_temperature

#######################################
# RICHARDSON NUMBER                   #
#######################################
def bulk_richardson_number(
        thv, U, z, ref_thv, ref_U=0., ref_z=0., g=_g,
        mode='caporaso', assume_sorted=False,
        ):
    """Return a float or an array.

        mode == 'caporaso' uses Eq. 3.2 in Caporaso 2012.

        Parameters
        ----------
        thv : float or array
            (K) virtual potential temperature
        U : float or array
            (m/s) horizontal wind speed
        z : float or array
            (m) elevation
        ref_thv : float
            (K) virtual potential temperature at reference level
        ref_U : float, optional
            (m/s) horizontal wind speed at reference level. Default: 0.
        ref_z : float, optional
            (m) elevation of reference level. Default: 0.
        g : float, optional
            (m/s2) gravitational acceleration. Default is the value of _g of
            this module
        mode : str, optional
            'caporaso', 'cosmo' or 'shao'. Default: 'caporaso'
        assume_sorted : bool, optional
            set this to True, if z is sorted in ascending order. Default: False

        Returns
        -------
        float or array
            The bulk Richardson number

        Note
        ----
        If thv, U, z are given as an array, they must all have the same
        length. mode 'cosmo' only works with arrays.
    """
    if isinstance(thv, collections.Iterable):
        assert len(thv) == len(U) == len(z)

    # ALIASES 
    if mode == 'igmk':
        mode = 'caporaso'

    # sorting and averaging:
    if mode == 'cosmo' and not assume_sorted:
        # sorting:
        idx = np.argsort(z, 0)  # sort index

        # handle case if z is 2D:
        if len(np.shape(z)) == 2:
            idx = idx[:, 0]

        ridx = np.argsort(idx)   # reverse index
        z    = z[idx]
        U    = U[idx]
        thv  = thv[idx]

        # mean virt. pot. temp:
        use = -np.isnan(thv)        # only non-nan position are used
        thv_mod = copy(thv)
        thv_mod[-use] = 0.          # replace nan's by 0
        thv_mean = np.cumsum(thv_mod, 0) / np.cumsum(use, 0)   # mean

        
    # finite differences:
    D_thv = thv - ref_thv  # Delta thv
    D_U   = U - ref_U      # Delta U
    D_z   = z - ref_z      # Delta z
    

    # bulk Richardson number:
    if   mode == 'caporaso':  Rib = g * D_thv / ref_thv  * D_z /   U**2
    elif mode == 'cosmo':     Rib = g * D_thv / thv_mean * D_z /   U**2
    elif mode == 'shao':      Rib = g * D_thv / thv      * D_z / D_U**2
    else: raise NotImplementedError('unknown mode: ' + mode)

    # reverse sorting:
    if mode == 'cosmo' and not assume_sorted:
        Rib = Rib[ridx]

    return Rib

#######################################
# TESTING                             #
#######################################
if __name__ == '__main__':
    L = 6.5e-3
    zmin = 0
    zmax = 10000
    N    = 40
    alpha = 1e-3
    zt = equal_transmission_range(zmin=zmin, zmax=zmax, N=N, alpha=alpha)
    zm = equal_mass_range(zmin=zmin, zmax=zmax, N=N)
    plt.plot(zt, 'bo')
    plt.plot(zm, 'ro')
    plt.show()

#######################################
# OUTDATED FUNCTIONS                  #
#######################################

def water_vapor_saturation_pressure_Buck1981(T, p=1e5, phase='liquid'):
    """*Return water vapor saturation pressure.

        DEACTIVATED ON 2016-01-08.
        
        Takes and returns SI base units.
        
        Uses e_w4 and f_w5 for water and e_i3 and f_i4 for ice as described in
        Arden L. Buck "New Equations for Computing Vapor Pressure and
        Enhancement Factor", Am.  Met. Soc. 1981.
        
        Parameters
        ----------
        T : float
            temperature (K)
        p : float, default: 1e5
            pressure (Pa)
        phase : str
                'liquid' or 'solid'. Default is 'liquid'

        Returns
        -------
        float : water vapor pressure (Pa).

        Accuracy 
        --------
        With respect to liquid water:
        Less than 0.26% within
        T in [-40, 50] deg C and 
        p in [250, 1000] hPa.

        With respect to ice:
        Less than 0.08% within
        T in [-80, 0] deg C and
        p in [250, 1000] hPa.

        Chemistry
        ---------
        Assumes Earth chemical composition of the Earth's troposphere.
    """
    raise NotImplementedError('Function does not yield the correct values.')
    ##############################
    # CONSTANTS                  #
    ##############################
    try:
        ph = phase[0].lower()
    except Exception:
        raise TypeError("phase must be 'liquid' or 'solid'.")
    # liquid
    if ph in ['l', 'w']:
        a = 6.1121
        b = 18.729
        c = 257.87
        d = 227.3
        A = 4.1e-4
        B = 3.48e-6
        C = 7.4e-10
        D = 30.6
        E = -3.8e-2
    # solid
    elif ph in ['s', 'i']:
        a = 6.1115
        b = 23.036
        c = 279.82
        d = 333.7
        A = 2.2e-4
        B = 3.83e-6
        C = 6.4e-10
        D = 0.
        E = 0.
    else:
        raise TypeError("phase must be 'liquid' or 'solid'.")

    ##############################
    # COMPUTATIONS               #
    ##############################
    """
        Units are converted to SI base units (Pa, K). Buck uses hPa and deg C.
        This is compensated for by multiplying pressure by 100 or 1/100 and
        adding +/- 273.15 to temperature.
    """
    TC = T - 273.15  # temperature in deg C
    P  = p * 1e-2    # pressure in hPa
    # enhancement factor:
    f = 1 + A + p * (B + C * (TC + D + E * P)**2)
    # pure water vapor pressure (Pa):
    e = a * np.exp((b - TC/d) * TC / (TC + c)) * 1e2
    # water vapor pressure in air:
    return  e * f
