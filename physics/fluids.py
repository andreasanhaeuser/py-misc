#!/usr/bin/python
"""A collection of functions related to fluid dynamics.

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
     
    u : zonal wind speed (m/s)
    v : meridional wind speed (m/s)
    w : vertical wind speed (m/s)
    U : horizontal wind speed (m/s)  [== sqrt(u**2 + v**2)]
    V : total wind speed (m/s)  [==sqrt(u**2 + v**2 + w**2)]

    T : temperature (K)

    Arrays
    ------
    Most functions which take scalars also work with array albeit not explicitly
    states. Just try!


    Author
    ------
    Written in 2016
    by Andreas Anhaeuser
    Insitute for Geophysics and Meteorology
    University of Cologne
    Germany
    <anhaeus@meteo.uni-koeln.de>
"""

# PyPI modules
import numpy as np

# local modules
from . import meteo as met

###################################################
# CONSTANTS                                       #
###################################################
# geophysical constants
_g = met._g        # (m s-2) Earth's graviational acceleration
_RE = 6371e3     # (m) Earth's radius

# thermodynamic constants
_Rd = met._Rd      # (J K-1 kg-1) specific gas constant for dry air
_cp_d = met._cp_d  # (J K-1 kg-1) isochoric specific heat capacity of dry air
_L = met._L        # (J kg-1) specific latent heat for vaporization of water

###################################################
# FUNCTIONS                                       #
###################################################
def mass_stream_function_yz(v, p, lat=None, mode='local'):
    """Return an array.
        
        Computes the mass stream function in meridional-vertical direction
        (assuming symmetry in zonal direction).

        Shapes and axes
        ---------------
        - Vertical axis must be last axis in v and p.
        - Shape of p must be compatible for multiplication with v.
        - Shape of lat must be compatible for multiplication with p.

        Parameters
        ----------
        v : array, at least 2d
            (m/s) meridional wind
        p : array, at least 1d
            (Pa), pressure
        lat : array, at least 2d, optional
            (deg) latitude. Must be given in 'global' mode.
        mode : {'local', 'global'}, optional
            if 'global', then the result is integrated in zonal direction,
            otherwise not. Default: 'local'.

        Returns
        -------
        psi : array, same shape as `v`
            mass stream function
    
        Units
        -----
        In 'local' mode, units of psi are (kg m-1 s-1), in 'global' mode, they
        are (kg s-1).

        Author
        ------
        Written in 2016
        by Andreas Anhaeuser
        Institute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """
    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    Sv = np.shape(v)
    Sp = np.shape(p)
    assert len(Sp) == 1
    assert len(np.shape(v)) > 1
    assert Sv[-1] == Sp[-1]
    assert mode in ['local', 'global']

    ###################################################
    # INTEGRATE                                       #
    ###################################################
    # pressure increment
    dp = np.nan * np.empty(Sp)
    dp[1:] = (p[1:] - p[:-1])
    dp[0] = 0
    
    assert np.max(dp) <= 0

    # mean wind between levels
    v_mean = level_mean(v)

    integral = 1./_g *  np.cumsum(v_mean * dp, -1)

    ###################################################
    # RETURN                                          #
    ###################################################
    if mode == 'local':
        return integral
    elif mode == 'global':
        return 2 * np.pi * _RE * np.cos(np.pi / 180 * lat) * integral
    
def energy_stream_function_yz(v, T, s, p, z=0, U=0, lat=None, mode='local'):
    """Return an array.
        
        Computes the energy stream function in meridional-vertical direction
        (assuming symmetry in zonal direction). Only sensible and latent heat
        are taken into account.

        Shapes and axes
        ---------------
        - Vertical axis must be last axis in v and p.
        - p, T and q must have same shapes.
        - Shape of p must be compatible for multiplication with v.
        - Shape lat must be compatible for multiplication with p.

        Parameters
        ----------
        v : array, at least 2d
            (m s-1) meridional wind
        T : array, same shape as v
            (K) temperature
        s : array, same shape as v
            (kg kg-1) specific humidity
        z : array, same shape as v, optional
            (m) altitude of the pressure level. If not given, potential energy
            is neglected.
        U : array, same shape as v, optional
            (m s-1) total wind speed. If not given, kinetic energy is
            neglected, same shape as v, optional
                (m s-1) total wind speed. If not given, kinetic energy is
                neglected.
        p : array, at least 1d
            (Pa), pressure
        lat : array, at least 2d, optional
            (deg) latitude. Must be given in 'global' mode.
        mode : {'local', 'global'}, optional
            if 'global', then the result is integrated in zonal direction,
            otherwise not. Default: 'local'.

        Returns
        -------
        psi : array, same shape as `v`
            energy stream function
    
        Units
        -----
        In 'local' mode, units of psi are (J m-1 s-1), in 'global' mode, they
        are (J s-1).

        Author
        ------
        Written in 2016
        by Andreas Anhaeuser
        Institute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """
    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    Sv = np.shape(v)
    ST = np.shape(T)
    Ss = np.shape(s)
    SU = np.shape(U)
    Sz = np.shape(z)
    Sp = np.shape(p)

    if SU == ():
        U = 0 * v
    if Sz == ():
        z = 0 * v

    assert Sv == ST == Ss == SU == Sz
    assert len(np.shape(v)) > 1
    assert Sv[-1] == Sp[-1]
    assert mode in ['local', 'global']

    ###################################################
    # INTEGRATE                                       #
    ###################################################
    # pressure increment
    dp = np.nan * np.empty(Sp)
    dp[1:] = (p[1:] - p[:-1])
    dp[0] = 0

    assert np.max(dp) <= 0

    # mean values between levels:
    v_mean = level_mean(v)
    T_mean = level_mean(T)
    s_mean = level_mean(s)
    z_mean = level_mean(z)
    U_mean = level_mean(U)

    # energy components (devided by rho)
    E_sens = _cp_d * T_mean      # sensible heat
    E_lat = _L * s_mean          # latent heat
    E_pot = _g * z_mean          # potential energy
    E_kin = 0.5 * U_mean**2      # kinetic energy

    # total energy
    E_tot = E_sens + E_lat + E_pot + E_kin

    # energy flux
    integrand = E_tot * v_mean

    # integral
    integral = 1./_g *  np.cumsum(integrand * dp, -1)

    ###################################################
    # RETURN                                          #
    ###################################################
    if mode == 'local':
        return integral
    elif mode == 'global':
        return 2 * np.pi * _RE * np.cos(np.pi / 180 * lat) * integral

###################################################
# HELPER FUNCTIONS                                #
###################################################
def level_mean(x):
    """*Helper function."""
    x_mean = np.swapaxes(x, -1, 0)
    x_mean[1:] = 0.5 * (x_mean[1:] + x_mean[:-1]) 
    x_mean[0] = 0.
    x_mean = np.swapaxes(x_mean, 0, -1)
    return x_mean
