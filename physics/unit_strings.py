def get_nice_unit_str(unit):
    convert = {
            'ug m-3' : '\mug m^{-3}',
            }
    if unit in convert:
        return convert[unit]
    return unit
