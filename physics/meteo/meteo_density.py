# PyPI modules
import numpy as np

# constants
_avogadro_constant = 6.02214076e23  # mol-1
_boltzmann_constant = 1.380649e-23  # (J K-1)
_dobson_unit = 2.2687e20            # m-2

# defaults
_pressure = 1013.25e2   # N m-2
_temperature = 288.15   # K

def ppb_to_ug(ppb, substance, p=None, T=None):
    M = molecular_mass(substance)
    v = ppb * 1e-9
    mass_si = volume_to_mass_density(v, M, p, T)
    mass_ug = mass_si * 1e9
    return mass_ug

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

def dobson_unit_factor(substance):
    """Return ratio of DU/SI numeric values."""
    m = molecular_mass(substance)
    return 1 / (m * _dobson_unit)

def ppb_factor(substance, p=None, T=None):
    M = molecular_mass(substance)
    mass_to_volume = 1. / volume_to_mass_density(1, M, p, T)
    mass_to_ppb = 1e9 * mass_to_volume
    return mass_to_ppb

def molecular_mass(substance):
    """Return molecular mass in kg."""
    name_short = get_chemical_formula(substance).lower()
    masses_mol = {
            'co'  : 28.010e-3,
            'no2' : 46.006e-3,
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
            'nitrogendioxide' : 'NO2',
            'sulfurdioxide' : 'SO2',
            }
    substance_lower = substance.lower()
    if substance_lower in aliases:
        return aliases[substance_lower]
    return substance
