#!/usr/bin/python3
"""Translation utils."""

# standard modules
import os

# local modules
from misc.in_out import tables as tab

def translate(
        term, filename=None, dictionary=None, term_end=':', sep='|',
        ignore_case=False, raise_error_if_no_result=False, **kwargs
        ):
    """Return translation as str or list of such.

        Parameters
        ----------
        term : object
        filename : str, optional
            a text file where each line is a dictionary.
        dictionary : dict, optional
            keys are source language, items are the translations
        term_end : str, optional
            separates the source language term and it(s) translation in the
            dictionary file
        sep : str, optional
            separator between alternative translations in the dictionary file
        ignore_case : bool, optional
        raise_error_if_no_result: bool, optional
            Default: False
        **kwargs : keyword arguments
            they are passed directly to `utils_igkm.tables.get_namelist()`


        Returns
        -------
        (str) or (list of str) or (None)
            translation
            None if there isn't any


        Raises
        ------
        KeyError
            If `raise_error_if_no_result` is True and nothing is found


        Example
        -------
        An example dictionary file with `sep='|'` and `term_end=':'` :

            term one : translation
            term two : another translation
            term three : a list | of | three translations
    """
    # load -------------------------------------------
    if dictionary is None:
        if filename is None:
            raise ValueError('Both dictionary and filename are None.')
        args = filename, term_end, sep
        dictionary = load_dictionary(*args, **kwargs)

    # find -------------------------------------------
    for key in dictionary:
        # non case-sensitive
        if ignore_case and key.lower() == term.lower():
            return dictionary[key]

        # case-sensitive
        elif key == term:
            return dictionary[key]
    # -----------------------------------------------

    # if this line is reached, nothing has been found
    # -----------------------------------------------
    if raise_error_if_no_result:
        raise KeyError(
            "Cannot find term '%s' in dictionary file %s"
            % (term, filename)
            )

    return None

def load_dictionary(filename, term_end, sep, **kwargs):
    """Return a dict."""
    if not os.path.isfile(filename):
        raise OSError('File does not exist: %s' % filename)
    return tab.read_namelist(filename, name_end=term_end, sep=sep, **kwargs)
