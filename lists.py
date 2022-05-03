"""List manipulation tools."""

from collections import Iterable

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

def occurences(choice, element):
    """Return number of times that element appears in choice."""
    count = 0
    elements = [x for x in choice if x == element]
    return len(elements)

def contiguous_interval_indices(iterable, criterion=True):
    """Return indices of contigous subintervals with iterable[n]==criterion."""
    assert isinstance(iterable, Iterable)

    N = len(iterable)
    if criterion == True:
        good = iterable
    else:
        good = [element == criterion for element in iterable]

    intervals = []
    nlo = None
    for n in range(N):
        if not good[n]:
            continue

        if nlo is None:
            nlo = n

        complete = False
        if n == N-1:
            complete = True
        elif not good[n+1]:
            complete = True

        if not complete:
            continue

        interval = (nlo, n+1)
        nlo = None
        intervals.append(interval)

    return intervals

def contiguous_intervals(iterable, criterion=True):
    """Return a list of contigous subintervals with iterable[n]==criterion."""
    interval_indices = contiguous_interval_indices(iterable, criterion)
    intervals = []
    for nlo, nhi in interval_indices:
        interval = iterable[nlo:nhi]
        intervals.append(interval)
    return intervals
