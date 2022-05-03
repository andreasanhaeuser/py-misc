# PyPI modules
import numpy as np

# constants
_avogadro_constant = 6.02214076e23  # mol-1
_boltzmann_constant = 1.380649e-23  # (J K-1)

# defaults
_pressure = 1013.25e2    # N m-2
_temperature = 288.15    # K
_dobson_unit = 2.687e20  # m-2

def ppb_to_ug(ppb, substance, p=None, T=None):
    """Convert ppb to ug/m3."""
    M = molecular_mass(substance)
    v = ppb * 1e-9
    mass_si = volume_to_mass_density(v, M, p, T)
    mass_ug = mass_si * 1e9
    return mass_ug

def ug_to_ppb(mass_ug, substance, p=None, T=None):
    """Convert ug/m3 to ppb."""
    factor = ppb_factor_ug(substance, p, T)
    ppb = mass_ug * factor
    return ppb

def volume_to_mass_density(v, M, p=None, T=None):
    """

        Parameters
        ----------
        v : float
            (m3 m-3) volume density
        M : float
            (kg) molecular mass
        p : float
            (N m-2) pressure
        T : float
            (K) temperature

        Returns
        -------
        rho : float
            (kg m-3) mass density
    """
    # defaults ============================================
    if p is None:
        p = _pressure
    if T is None:
        T = _temperature
    # =====================================================

    rho = p * M / (T * _boltzmann_constant) * v

    return rho

def si_to_dobson_unit(value, substance):
    """Convert from kg m-2 to Dobson unit."""
    m = molecular_mass(substance)
    N_molec = value / m             # (m-2) molecules per surface area
    return N_molec / _dobson_unit

def dobson_unit_to_si(value, substance):
    """Convert from Dobson unit to kg m-2."""
    m = molecular_mass(substance)
    return value * m * _dobson_unit

def dobson_unit_factor(substance):
    """Return ratio of DU/SI numeric values."""
    m = molecular_mass(substance)
    return 1 / (m * _dobson_unit)

def ppb_factor_ug(substance, p=None, T=None):
    factor_si = ppb_factor(substance, p, T)
    factor_ug = factor_si * 1e-9
    return factor_ug

def ppb_factor(substance, p=None, T=None):
    """Factor for ppb = mass_si * ppb_factor."""
    M = molecular_mass(substance)
    mass_to_volume = 1. / volume_to_mass_density(1, M, p, T)
    mass_to_ppb = 1e9 * mass_to_volume
    return mass_to_ppb

def molecular_mass(substance):
    """Return molecular mass in kg."""
    name_short = get_chemical_formula(substance).lower()
    masses_mol = {
            'co'  : 28.010e-3,
            'no' : 30.006e-3,
            'no2' : 46.006e-3,
            'o3'  : 47.997e-3,
            'no3' : 62.004e-3,
            'so2' : 64.066e-3,
            'so4' : 96.06e-3,
            }

    if name_short not in masses_mol:
        raise NotImplementedError(
                'Substance not implemented: %s (%s)'
                % (substance, name_short)
                )

    return masses_mol[name_short] / _avogadro_constant

def get_chemical_formula(substance):
    aliases = {
            'carbonmonoxide' : 'CO',
            'nitricoxide' : 'NO',
            'nitric oxide' : 'NO',
            'nitrogendioxide' : 'NO2',
            'ozone' : 'O3',
            'sulfurdioxide' : 'SO2',
            }
    substance_lower = substance.lower()
    if substance_lower in aliases:
        return aliases[substance_lower]
    return substance
