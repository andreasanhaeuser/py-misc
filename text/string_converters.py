"""String converters.

    Author
    ------
    Andreas Anhaeuser (AA) <andreas.anhaeuser@posteo.net>
"""

# standard modules
import collections

# PyPI modules
import numpy as np

# constants
_missing_str = ('nan', '')

def string_to_number(s, dtype=float, missing_str=None, missing_value=np.nan):
    """Convert a str to number.

        Parameters
        ----------
        s : str or iterable of such
        dtype : class, optional
            Default: float. number type
        nan_str : None or str or iterable of such, optional
            Default: 'nan'. str to be interpreted as nan

        Returns
        -------
        number or array of such

        Raises
        ------
        ValueError
            At least one element is not convertible
    """
    if missing_str is None:
        missing_str = _missing_str

    ############################################################
    # scalar                                                   #
    ############################################################
    # recursively call function
    if isinstance(s, str):
        numbers = string_to_number([s], dtype)
        return numbers[0]

    ############################################################
    # input check                                              #
    ############################################################
    if not isinstance(s, collections.Iterable):
        raise TypeError('Need str or iterable of such, got %s.' % s)

    ############################################################
    # standardize                                              #
    ############################################################
    if isinstance(missing_str, str):
        missing_words = [missing_str]
    else:
        missing_words = missing_str

    ############################################################
    # main                                                     #
    ############################################################
    N = len(s)
    out = np.empty(N, dtype=dtype)

    for n in range(N):
        word = s[n]

        # nan
        if word in missing_words:
            out[n] = missing_value
            continue

        out[n] = dtype(word)

    return out
