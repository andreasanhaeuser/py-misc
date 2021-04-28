def is_prime(x):
    if not isinstance(x, int):
        raise TypeError('x must be int, got %s.' % x)

    if x < 1:
        raise ValueError('x must be positive, got %i.' % x)

    if x == 1:
        return False

    factor = 2
    while factor**2 <= x:
        if x % factor == 0:
            return False
        factor += 1
    return True
