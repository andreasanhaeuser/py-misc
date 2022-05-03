def convert_to_si(value, units):
    factor = get_scale_factor(units)
    return value * factor

def convert_from_si(value, units):
    factor = get_scale_factor(units)
    return value / factor

def get_scale_factor(units):
    """Return a float."""
    known_units = {
            # distance
            'nm' : 1e-9,
            'um' : 1e-6,
            'mm' : 1e-3,
            'm' : 1,
            'km' : 1e3,

            # mass
            'ug' : 1e-9,
            'mg' : 1e-6,
            'g' : 1e-3,
            'kg' : 1,
            't' : 1e3,
            'kt' : 1e6,
            'Mt' : 1e9,
            'Gt' : 1e12,

            # time
            'fs' : 1e-15,
            'ps' : 1e-12,
            'ns' : 1e-9,
            'us' : 1e-6,
            'ms' : 1e-3,
            's' : 1,
            'min' : 60,
            'h' : 3600,
            'd' : 86400,

            # concentration
            '[\\mug/m^3]' : 1e-9,
            'ug m-3' : 1e-9,
            'kg m-3' : 1.,
            }

    if units in known_units:
        return known_units[units]

    raise ValueError('Unknown units: %s' % units)

def si_units(units):
    """Return a str."""
    known_units = {
            # distance
            'nm' : 'm',
            'um' : 'm',
            'mm' : 'm',
            'm' : 'm',
            'km' : 'm',

            # mass
            'ug' : 'kg',
            'mg' : 'kg',
            'g' : 'kg',
            'kg' : 'kg',
            't' : 'kg',
            'kt' : 'kg',
            'Mt' : 'kg',
            'Gt' : 'kg',

            # time
            'fs' : 's',
            'ps' : 's',
            'ns' : 's',
            'us' : 's',
            'ms' : 's',
            's' : 's',
            'min' : 's',
            'h' : 's',
            'd' : 's',

            # concentration
            '[\\mug/m^3]' : 'kg m-3',
            'ug m-3' : 'kg m-3',
            'kg m-3' : 'kg m-3',
            }

    if units in known_units:
        return known_units[units]

    raise ValueError('Unknown units: %s' % units)
