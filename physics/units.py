def convert_to_si(value, units):
    factor = get_scale_factor(units)
    return value * factor

def convert_from_si(value, units):
    factor = get_scale_factor(units)
    return value / factor

def get_scale_factor(units):
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
            'ug m-3' : 1e-9,
            'kg m-3' : 1.,
            }

    if units in known_units:
        return known_units[units]

    raise ValueError('Unknown units: %s' % units)
