#!/usr/bin/python2
"""Identify mesoscale convective systems (MCS) in model output.

    Uses an intensity threshold for outgoing longwave radiation (OLR).

    Areas are identified as an MCS if:
        - brightness temperature (Tb) computed from OLR is lower than a certain
          threshold AND
        - a minimum number of such pixels are connected AND
        - the feature persist for a minimum number of time steps


    Dependencies
    ============
    - python 2 (possibly also works with python 3) with these modules:
        * numpy
        * scipy
        * scipy.ndimage
        * netCDF4
    - netCDF
    - cdo (Climate Data Operators) [tested with version 1.6.2]


    Usage
    =====
    - If you want to use it from command line, type <scriptname> -h
    - If you want to use it as a module inside python, see docstring of main()

    
    Reliability
    ===========
    The code has not been tested heavily and most probably contains bugs.

    Works fine with WRF output files.


    Contribute
    ==========
    Some constraints are in principle not necessary and only exist because
    they are not implemented in the script yet.
    
    If you want to help improving the program, please feel free to do so! The
    author would be glad
        - to be informed
        - to receive the new (commented) code

        Ideas for development
        ---------------------
        - MCS detection across file borders
        - allow time axis to be at arbitray position


    Author
    ======
    Written in 2016
    by Andreas Anhaeuser
    Institute for Geophysics and Meteorology
    University of Cologne
    Germany
    <anhaeus@meteo.uni-koeln.de>
"""

# standard modules
import os
import sys
import argparse

# PyPI modules
import numpy as np
from scipy import ndimage
from scipy import constants as cst
from netCDF4 import Dataset

#========= DEFAULT ============================#
_sigma = cst.sigma   # (W/(m2 K4)) Stefan-Boltzmann constant
_Tmax = 213          # (K)
_min_size = 1        # minimum size in pixels
_min_time_steps = 0  # minimum number of time steps
#==============================================#

_doc = {}
_doc['min_size'] = """(in number of pixels) Use this if you want the area
covered by the MCS to have a minimum size. Default: %1.0f""" % _min_size

_doc['min_time'] = """(in number of time steps) Use this if you want the life
time of the MCS to have a minimum length. 0 means that even features that flash
up at only one time step are identified as MCS (i. e. no time information).  1
means that features that are detected at one time step have to PERSIST for at
least one more time step. Default: %1.0f""" % _min_time_steps

_doc['max_temp'] = """ (Kelvin) Maximum brightness temperature of the MCS.
Default: %1.0f""" % _Tmax

_doc['allow_diagonal_connection'] = """By default, pixels are considered to be
connected, only if they share a common line. Using this option will also
consider them to be connected if their share a common corner."""

_doc['file_in'] = """A netcdf file. It must contain: (1) a time axis ('time' /
'times' / 'Time' / 'Times'), (2) two horizontal axes, and (3) a 3d-variable
called 'OLR' with time as first axis. If you use the `min_time` option, it is
advisable to merge as many time steps as possible into one file. MCS which
extend over several files (temporally or spatially) are only detected in those
files where they fulfill the detection conditions independently of the other
files."""

_doc['file_out'] = """A netcdf file. May be identical to input file (use with
caution!). File only contains one variable 'MCS'. Output directory must
exist."""

_doc['description'] = __doc__.split('\n')[0]


###################################################
# MAIN                                            #
###################################################
def main(
        fni,
        fno,
        Tmax=_Tmax,
        min_size=_min_size,
        min_time_steps=_min_time_steps,
        allow_diagonal=False,
        ):
    """Identify MCS areas and save to file.

        Parameters
        ----------
        fni : str
            name of input file
        fno : str
            name of output file
        Tmax : float, optional
            (K) maximum brightness temperature of outgoing longwave radiation.
            Default: `_Tmax`
        min_size : int, optional
            minimum number of connected pixels. Default: `_min_size`
        min_time_steps : int, optional
            minimum number of time intervals that the feature has to persist
            after its first appearance. Default: `_min_time_steps`
        allow_diagonal : bool, optional
            if True, diagonnaly neighbouring pixels also count as 'connected'.
            Default: False

        Returns
        -------
        None
    """
    ###################################################
    # INITIALIZE                                      #
    ###################################################
    # Here, it is defined what is a neighbour.
    # 1's are cells that count as 'neighbours'.
    if allow_diagonal:
        pattern = np.array([[1, 1, 1,],
                            [1, 1, 1,],
                            [1, 1, 1,]])
    else:
        pattern = np.array([[0, 1, 0,],
                            [1, 1, 1,],
                            [0, 1, 0,]])

    # brightness temperature --> outgoing longwave radiation (OLR):
    olr_max = _sigma * Tmax**4

    ###################################################
    # PATHS                                           #
    ###################################################
    # input file
    assert os.path.isfile(fni), 'Input file does not exist.'

    # output directory
    idx = fno.rfind('/')
    if idx > 0:
        path_out = fno[:idx]
        assert os.path.isdir(path_out), 'Output directory does not exist.'

    # temp files
    fn_tmp = fno + '.tmp'

    ###################################################
    # CREATE TEMP FILE                                #
    ###################################################
    # extract OLR -- temp file
    command = 'cdo -s selname,OLR %s %s' % (fni, fn_tmp)
    status = os.system(command)
    assert status == 0

    ###################################################
    # READ OLR                                        #
    ###################################################
    with Dataset(fn_tmp, 'a') as fid:
        vars = fid.variables
        varid = vars['OLR']

        # read OLR data
        olr = varid[:]

        # shape of OLR
        message = \
        '''OLR must currently have exactly three dimensions. If you assume this
        is an unnecessary constraint in the program and you know how to code in
        python, please feel free to extend the program and contact the
        author!'''
        S = np.shape(olr)
        assert len(S) == 3, message

        ###################################################
        # FIND TIME AXIS                                  #
        ###################################################
        # This section is intended for a generalisation of the program.
        # `taxis`, the position of the time axis (which is determined here) may
        # be used to allow the time axis of OLR at arbitrary position in future
        # versions of this program. Contributions welcome!

        timekey = None
        for key in vars.keys():
            if 'time' in key.lower():
                timekey = key
                times = vars[timekey][:]
                break
        if timekey is None:
            raise NotImplementedError('')

        # time axis in OLR
        assert timekey in varid.dimensions
        taxis = varid.dimensions.index(timekey)

        message = \
        '''Currently, time needs to be the first dimension. If you find this is
        an unnecessary constraint in the program and you know how to code in
        python, please feel free to extend the program and contact the
        author!'''
        assert taxis == 0, message

        # number of time steps
        N = S[taxis]

        ###################################################
        # FIND MCS                                        #
        ###################################################
        # ========== INITIALIZE ====================== #
        # (potential) MCS areas are individually labeled by numbers
        # label 0 means 'not an MCS'
        labels = np.zeros(S, dtype=int)

        # number of MCS candidates
        number = 0

        # List of age of each MCS candidate (will be extended in the for-loop).
        # The position in this list refers to the label number of the MCS. MCS
        # label numbers start at 1.  For easier indexing lateron, the 'no
        # mcs'-label (position 0) also gets assigned an age.  This remains
        # unused.
        age = [0]
        # ============================================ #

        # loop over time
        for n in range(N):
            # mask candidate areas
            mask = olr[n] <= olr_max
            labels_t, number_t = ndimage.label(mask, pattern)

            # loop over candidate areas
            for label_t in range(1, number_t+1):
                ###################################################
                # SIZE                                            #
                ###################################################
                idx = labels_t == label_t
                if np.sum(idx) < min_size:
                    labels_t[idx] = 0
                    continue

                ###################################################
                # FIND OVERLAP                                    #
                ###################################################
                # overlap with previous time step:
                if n == 0:
                    overlap = 0
                else:
                    overlap = labels[n-1] * idx

                label = np.max(overlap)  # either the label or zero (if no
                                         # overlap)
                # create new label:
                if label == 0:
                    number += 1
                    label = number
                    age.append(0)
                else:
                    age[label] += 1

                labels[n][idx] = label
                
        ###################################################
        # CHECK AGE                                       #
        ###################################################
        mcs = np.zeros(S)
        label_min = 1
        label_max = 1
        for n in range(N):
            slab = labels[n]
            mcs_slab = np.zeros(S[1:])
            labels_slab = np.unique(slab)
            for label in labels_slab:
                if label == 0:
                    continue
                if age[label] < min_time_steps:
                    continue
                idx = slab == label
                mcs_slab[idx] = 1

            mcs[n] = mcs_slab

        vars['OLR'][:] = mcs

        ###################################################
        # MODIFY VARIABLE ATTRIBUTES                      #
        ###################################################
        varid.setncattr('units', '')
        varid.setncattr('long_name', 'mesoscale convective system')
        varid.setncattr('standard_name', '')
        varid.setncattr('description', 'MCS detected from OLR')
        varid.setncattr('comment', '1: yes, 0: no')
        varid.setncattr('comment_a', 'minimum size: %1.0f px' % min_size)
        varid.setncattr('comment_b', 'minimum life time: %1.0f step(s)' %
                min_time_steps)
        varid.setncattr('comment_c', 'maximum brightness temperature: %1.0f px' % 
                Tmax)

    ###################################################
    # RENAME VARIABLE                                 #
    ###################################################
    # change variable name
    command = 'cdo -s chname,OLR,MCS %s %s' % (fn_tmp, fno)

    status = os.system(command)
    assert status == 0


if __name__ == '__main__':
    ###################################################
    # PARSE COMMAND LINE ARGUMENTS                    #
    ###################################################
    parser = argparse.ArgumentParser(
            description=_doc['description'],
            epilog=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            )
    parser.add_argument('file_in', nargs=1, help=_doc['file_in'])
    parser.add_argument('file_out', nargs=1, help=_doc['file_out'])
    parser.add_argument('-D', '--allow_diagonal_connection',
            action='store_true', dest='allow_diagonal',
            help=_doc['allow_diagonal_connection'])
    parser.add_argument('-S', '--min_size', nargs='?', type=int,
            default=_min_size, help=_doc['min_size'])
    parser.add_argument('-t', '--min_time', nargs='?', type=int,
            default=_min_time_steps,
            help=_doc['min_time'])
    parser.add_argument('-T', '--max_temp', nargs='?', type=float,
            default=_Tmax, help=_doc['max_temp'])
    args = parser.parse_args()

    ###################################################
    # CALL MAIN                                       #
    ###################################################
    fni = args.file_in[0]
    fno = args.file_out[0]

    if fni == fno:
        proceed = raw_input(
            'Caution: input and output files are identical. ' +
            'If you continue, you will lose\n' + ' '*9 +
            'all other data in the input file. ' +
                'Do you want to continue? [y/n]: ')
        if proceed.lower() != 'y':
            print('Operation aborted by user.')
            sys.exit(1)

        print('Ok, assuming you know what you are doing...')

    main(fni=fni, fno=fno, Tmax=args.max_temp,
            min_size=args.min_size, min_time_steps=args.min_time,
            allow_diagonal=args.allow_diagonal)

