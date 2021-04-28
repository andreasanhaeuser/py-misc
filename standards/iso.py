#!/usr/bin/python2
"""Iso-formats."""

# local modules
import string

# A0 satisfies two conditions:
# 1)  long * short == 1        [area is 1 square meter]
# 2)  long / short == sqrt(2)  [when cutting in half, edge ratio remains equal]
_a0_long_edge = 2**0.25             # (m) long edge of A0 paper

_exponent = {
        'a' : 0,
        'e' : 1. / 16,
        'c' : 1. / 8, 
        'g' : 3. / 16,
        'b' : 1. / 4,
        'f' : 5. / 16,
        'd' : 3. / 8,
        'h' : 7. / 16,
        }

def paper_size(series='A', size=4, orientation='portrait', units='m'):
    """Return ISO 216 A-series format.

        Parameters
        ----------
        series : str of length 1, optional
            'A'...'H'. Default: 'A'
        size : int or float, optional
            position in the series. Default: 4
        orientation : str, optional
            'p', 'portrait', 'l', 'landscape'. Default: 'portrait'
        units : str, optional
            {'m', 'mm', 'inch'} units of the returned values. Default: 'm'

        Returns
        -------
        width : float
            width in m
        height : float
            height in m

        Notes
        -----
        Caution: By default, returned values are in meters!
        Function works also with negative and float input for `size`.

        Author
        ------
        Andreas Anhaeuser (AA) <anhaeus@meteo.uni-koeln.de>
        Institute for Geophysics and Meteorology
        University of Cologne, Germany

        History
        -------
        2017-11-15 (AA): Created
    """
    # input check
    if not isinstance(series, str):
        raise TypeError('`series` must be str.')
    series = series.lower()
    allowed = list(string.ascii_lowercase)[:8]
    if not series in allowed:
        raise ValueError('`series` must be in %s', str(allowed))
    if not isinstance(orientation, str):
        raise TypeError('`orientation` must be str.')
    orientation = orientation[:1].lower()
    if orientation not in ('p', 'l'):
        raise ValueError()

    # unit conversion
    if units == 'm':
        uc = 1.
    elif units == 'mm':
        uc = 1000.
    elif units == 'inch':
        uc = 1000 / 25.4
    else:
        raise ValueError()

    # edge lengths / series
    long_edge = _a0_long_edge / 2**(0.5 * size) * 2**_exponent[series] * uc
    short_edge = long_edge / 2**0.5


    # orientation
    if orientation == 'p':
        return (short_edge, long_edge)
    else:
        return (long_edge, short_edge)


###################################################
# TESTING                                         #
###################################################
if __name__ == '__main__':
    for series in string.ascii_lowercase[:8]:
        print ('%s4' % series), paper_size(series, 4)
    print(paper_size('A', 4.5))
    print(paper_size('A', -4.5))
