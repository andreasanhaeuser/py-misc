#!/usr/bin/python
"""Humidity measure conversions.

    q : float
        (kg m-3) volumetric (absolute) humidity
    r : float
        relative humidty (1: saturated)
    s : float 
        (kg kg-1) specific humidity (mass of water per total air mass)
    m : float
        (kg kg-1) mixing ratio (mass of water devided by remaining mass)

        from --->
      t     | q  r  s  m
      o     +------------
      |   q |    1  5  7  
      |   r | 2     9 11
      V   s | 6 10     3
          m | 8 12  4

    Accuracy
    --------
    Precision is only limited by water vapor saturation pressure.
    All other computations are precise (within numerical precision).

    Arrays
    ------
    All functions also work with arrays.

    Tested
    ------
    All humidity conversion functions have been tested and yield correct
    results.
"""

# PyPI modules
import numpy as np

# local modules
from . import meteo as met

_epsilon = met._epsilon
sat_press = met.water_vapor_saturation_pressure

def q_from_r(T, r):                    # 1
    """Convert relative humidity to volumetric humidity."""
    es = sat_press(T)
    return es * r / (_Rw * T)

def r_from_q(T, q):                    # 2
    """Convert volumetric humidity to relative humidity."""
    es = sat_press(T)
    return _Rw * T * q / es

def s_from_m(m):                       # 3
    """Convert mixing ratio to specific humidity."""
    return m / (1 + m)

def m_from_s(s):                       # 4 
    """Convert to specific humidity to mixing ratio."""
    return s / (1 - s)

def q_from_s(p, T, s):                 # 5
    """Convert to specific humidity to volumetric humidity."""
    r = r_from_s(p, T, s)
    return q_from_r(T, r)

def s_from_q(p, T, q):                 # 6
    """Convert volumetric humidity to specific humidity."""
    r = r_from_q(T, q)
    return s_from_r(p, T, r)
        
def q_from_m(p, T, m):                 # 7
    """Convert mixing ratio to volumetric humidity."""
    r = r_from_m(p, T, m)
    return q_from_r(T, r)

def m_from_q(p, T, q):                 # 8
    """Convert volumetric humidity to mixing ratio."""
    s = s_from_q(p, T, q)
    return m_from_s(s)

def r_from_s(p, T, s):                 # 9
    """Convert specific humidity to relative humidity."""
    es = sat_press(T)
    return p / es  / (_epsilon / s - _epsilon + 1)

def s_from_r(p, T, r):                # 10
    """Convert relative humidity to specific humidity."""
    es = sat_press(T)
    return _epsilon / (p / (r * es) - 1 + _epsilon)
   
def r_from_m(p, T, m):                # 11
    """Convert mixing ratio to relative humidity."""
    s = s_from_m(m)
    return r_from_s(p, T, s)

def m_from_r(p, T, r):                # 12
    """Convert relative humidity to mixing ratio."""
    es = sat_press(T)
    s = s_from_r(p, T, r)
    return m_from_s(s)
