"""Tools for handling command-line arguments."""
# standard
import os
import sys
import argparse
from copy import deepcopy as copy

# misc
import misc.in_out.tables as tab

def init_parser(setup_file=None):
    """Return an ArgumentParser with 'setup_file' as first argument."""
    if setup_file is None:
        calling_script = os.path.realpath(sys.argv[0])
        dirname = os.path.dirname(calling_script)
        setup_file = dirname + '/setup.txt'

    # initialize
    fmt_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(
            description=__doc__, formatter_class=fmt_class,
            )


    # setup file
    if setup_file is not None:
        setup_help = (
                'Setup file with default parameters. Command-line parameters'
                + ' override setup file parameters.'
                )
        parser.add_argument(
                dest='setup_file', nargs='?', default=setup_file,
                help=setup_help,
                )

    return parser

def get_setup(parser, *args, convert_to_number=True, **kwargs):
    """Parse command-line args and extend by setup file if present.

        If parser contains 'setup_file', this is loaded. Command-line options
        override setup file options.

        Parameters
        ----------
        parser : argparse.ArgumentParser or argparse.Namespace or dict
            command-line argument parser
            or parsed commmand-line arguments
            or conversion thereof into a dict

        Returns
        -------
        setup : dict
    """
    # ArgumentParser -> Namespace
    # ------------------------------------------------
    if isinstance(parser, argparse.ArgumentParser):
        cl_args = parser.parse_args()
    else:
        cl_args = parser
    # ------------------------------------------------

    # Namespace -> dict
    # ------------------------------------------------
    if isinstance(cl_args, argparse.Namespace):
        setup_cl = cl_args.__dict__
    else:
        setup_cl = cl_args
    # ------------------------------------------------

    # Check type (dict)
    # ------------------------------------------------
    if not isinstance(setup_cl, dict):
        raise TypeError('Cannont handle parse type: %s' % type(parser))
    # ------------------------------------------------

    # get base setup
    # ------------------------------------------------
    setup_base = {}
    setup_file = None
    if 'setup_file' in setup_cl:
        setup_file = setup_cl['setup_file']
    if setup_file is not None:
        setup_base = tab.read_namelist(
                setup_file, *args, convert_to_number=convert_to_number,
                **kwargs
                )
    # ------------------------------------------------

    # supersede with command-line arguments
    setup = supersede_setup(setup_base, setup_cl)

    return setup

################################################################
# helpers                                                      #
################################################################
def supersede_setup(
        setup_base, setup_cl, skip_none=True, skip_nan=True, inplace=True,
        ):
    """Override values from setup file with those in command_line_setup.

        Parameters
        ----------
        setup_base : dict
        setup_cl : dict
            items in here override those of `setup_base`
        skip_none : bool, optional
            default: True. Ignore None-valued items in `setup_cl`
        skip_nan : bool, optional
            default: True. Ignore NaN-valued items in `setup_cl`
        inplace : bool, optinal
            default: True. If True, operate on setup_base, otherwise return a
            new dict

        Returns
        -------
        setup : dict
            Combination of `setup_base` and `setup_cl` with precedence to the
            latter.
    """
    # ArgumentParser -> Namespace
    # ------------------------------------------------
    if isinstance(setup_cl, argparse.ArgumentParser):
        setup_cl = parser.parse_args()
    # ------------------------------------------------

    # Namespace -> dict
    # ------------------------------------------------
    if isinstance(setup_cl, argparse.Namespace):
        setup_cl = setup_cl.__dict__
    # ------------------------------------------------

    if inplace:
        setup = setup_base
    else:
        setup = copy(setup_base)

    for key in setup_cl:
        value = setup_cl[key]

        # None
        if skip_none and value is None:
            continue

        # NaN
        if skip_nan and (value != value):
            continue

        setup[key] = value

    return setup
