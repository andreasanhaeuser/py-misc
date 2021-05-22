from collections import Iterable

def is_increasing(x):
    """Return a bool if strictly increasing."""
    if not isinstance(x, Iterable):
        raise TypeError('x must be iterable, got %s' % type(x))

    N = len(x)
    for n in range(N-1):
        if x[n] < x[n+1]:
            continue
        return False
    return True
