# standard
from math import sqrt

def rectangular(x):
    """Return two int of which the product is x.
    
        Returns to ints `a` and `b`, with a >= b and a*b == x.
        They are those which are closest to a square, but can be as extreme as
        `x` and `1` if `x` is prime.

        Parameters
        ----------
        x: int > 0

        Returns
        -------
        a: int > 0
        b: int
            b > 0 and b <= a
    """
    # input check
    # =========================================
    if not isinstance(x, int):
        try:
            x = int(x)
        except TypeError:
            raise TypeError('x must be int, got %s' % type(x))
        except ValueError:
            raise ValueError('Cannot be converted to int: %s' % str(x))

    if not x > 0:
        raise ValueError('x must be positive, got %i' % x)
    # =========================================

    b = int(sqrt(x))
    while True:
        a = x // b
        if a*b == x:
            break
        b -= 1
    return a, b
