def is_prime(x):
    assert_positive_int(x)

    if x == 1:
        return False

    factor = 2
    while factor**2 <= x:
        if x % factor == 0:
            return False
        factor += 1
    return True

def prime_factors(x):
    """Return prime factors of x."""
    assert_positive_int(x)

    factor = 2
    remainder = x
    factors = []
    while factor*factor <= remainder:
        if remainder % factor == 0:
            factors.append(factor)
            remainder /= factor
        else:
            factor += 1
    if remainder > 1:
        factors.append(remainder)
    return sorted(factors)

################################################################
# helpers                                                      #
################################################################
def assert_positive_int(x):
    if not isinstance(x, int):
        raise TypeError('x must be int, got %s.' % x)

    if x < 1:
        raise ValueError('x must be positive, got %i.' % x)
