# standard modules
import collections

# PyPI modules
import numpy as np

def round_digits(x, digits=1, method='round'):
    """

        Parameters
        ----------
        x : float
        digits : int, optional
            (default: 1) number of significant digits
        method : {'round', 'floor', 'ceil'}, optional
            (default: 'round')

        Returns
        -------
        float
    """
    if digits <= 0:
        raise ValueError('`digits` must be a positive integer.')

    known_methods = ('round', 'floor', 'ceil')
    if method not in known_methods:
        raise ValueError('Unknown method: %s' % method)

    if isinstance(x, collections.Iterable):
        N = len(x)
        return np.array([round_digits(x[n], digits) for n in range(N)])

    if x == 0:
        return x

    magnitude_fractional = np.log10(np.abs(x))
    magnitude = int(np.floor(magnitude_fractional))
    decimals = digits - 1 - magnitude
    rounded = np.round(x, decimals)

    if method == 'round':
        return rounded

    if method == 'floor' and rounded <= x:
        return rounded

    if method == 'ceil' and rounded >= x:
        return rounded

    correction = 10**(-decimals)

    if method == 'floor':
        return rounded - correction

    if method == 'ceil':
        return rounded + correction
