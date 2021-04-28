"""List manipulation tools."""

def unique(x, sort=False):
    """Return list of unique elements.

        Parameters
        ----------
        x : iterable
        sort : bool, optional
            Default: False. Return a sorted list

        Returns
        -------
        list
            contains the unique elements of `x`
    """
    x_unique = []
    for element in x:
        if element in x_unique:
            continue
        x_unique.append(element)

    if sort:
        x_unique = sorted(x_unique)

    return x_unique
