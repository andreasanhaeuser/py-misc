#!/usr/bin/python3

# standard modules
import os
import sys
import argparse
import datetime as dt
import glob
import shutil

# misc
from misc import datetime_utils as dt_utils
from misc import tables as tab
from misc.string_utils import highlighted_string_list as hsl
from misc.chronometer import Chronometer
from misc.record_file import Record
from misc.filelock import FileLock as FL

################################################################
# main                                                         #
################################################################
def main(setup):
    calling_script = os.path.realpath(sys.argv[0])
    workdir = os.path.dirname(calling_script)

    table_name = workdir + '/' + setup['filename_jobs']
    record_file_copy = workdir + '/' + setup['filename_copy_record']
    timeout = setup['timeout_filelock']

    jobs = tab.read_vartable(table_name)
    Njobs = len(jobs['pattern_in'])

    header = 'File transfer services'
    Ndeleted = 0
    Ncopied = 0
    Ntotal = 0
    for njob in range(Njobs):
        stats = run_one_job(jobs, njob, setup)
        Ndeleted += stats['Ndeleted']
        Ncopied += stats['Ncopied']
        Ntotal += stats['Ntotal']



    print('\nSUMMARY')
    print(('%i' % Ntotal).rjust(4) + ' files')
    print(('%i' % Ncopied).rjust(4) + ' copied')
    print(('%i' % Ndeleted).rjust(4) + ' deleted')

def run_one_job(jobs, njob, setup):
    times_copy = get_times_copy(jobs)
    times_delete = get_times_delete(jobs)
    record_copy = Record(record_file_copy)
    patterns_in = jobs['pattern_in']

    actions = ('copy', 'delete')

    substance = 'job_%02i' % njob
    log_file = cron_utils.get_name_logfile(args, name, substance)

    hours_copy = float(jobs['hours_copy'][njob])
    hours_delete = float(jobs['hours_delete'][njob])

    pattern_in = patterns_in[njob]
    filenames_in = sorted(glob.glob(os.path.expanduser(pattern_in)))

    Nfiles = len(filenames_in)

    print('*' * 40)
    print('------------------')
    print('START JOB #%i.' % njob)

    Ndeleted = 0
    Ncopied = 0
    Ntotal = 0

    with Chronometer(Nfiles, file=log_file) as chrono:
        if log_file is None:
            chrono.use_object_print()
        sep = '\n                    '
        ages = ['%sh' % s for s in (hours_copy, hours_delete)]

        info = (''
                + 'Record file (copy): ' + record_file_copy
                + '\nInput pattern:      ' + hsl(patterns_in, njob, sep=sep)
                + '\nAction:             ' + hsl(actions, None)
                + '\nMinimum file age:   ' + hsl(ages, None)
                )
        chrono.set_header(header)
        chrono.set_info(info)
        chrono.set_item_name('file')
        chrono.show()

        # copy ================================================
        for filename_in in filenames_in:
            Ntotal += 1
            filename_in = os.path.expanduser(filename_in)
            filename_tmp = filename_in + '.ft.tmp'
            filename_out = get_filename_out(filename_in)
            print('---------------')
            print('Source:      %s' % filename_in)
            print('Temporary:   %s' % filename_tmp)
            print('Destination: %s' % filename_out)

            copy_one_file(filename_in, filename_out, hours_copy, record_copy)
            delete_one_file(filename_in, filename_out, hours_copy, record_copy)

            
            copy_file = check_whether_to_copy(
                    filename_in, filename_out, hours_copy, record_copy,
                    )
            delete_file = check_whether_to_delete(filename_in, hours_delete)

            if not copy_file and not delete_file:
                chrono.loop()
                continue

            with FL(filename_in, timeout) as fli, FL(filename_out, timeout) as flo:
                if not fli.attempt_to_lock():
                    print('Source locked -> continue')
                    chrono.loop()
                    continue

                if not flo.attempt_to_lock():
                    print('Destination locked -> continue')
                    chrono.loop()
                    continue

                #########################################################
                # copy                                                  #
                #########################################################
                if copy_file:
                    input_pattern_str = hsl(patterns_in, njob, sep=sep)
                    info = (''
                            + 'Record file (copy): ' + record_file_copy
                            + '\nInput pattern:      ' + input_pattern_str
                            + '\nAction:             ' + hsl(actions, 0)
                            + '\nMinimum file age:   ' + hsl(ages, 0)
                            )
                    chrono.set_info(info).show()

                    # create output directory
                    dirname_out = os.path.dirname(filename_out)
                    if not os.path.isdir(dirname_out):
                        print('Create directory: %s' % dirname_out)
                        os.makedirs(dirname_out)

                    # copy file
                    print('Copy source -> temp')
                    shutil.copyfile(filename_in, filename_tmp)
                    fli.release()

                    print('Copy temp -> dest')
                    shutil.copyfile(filename_tmp, filename_out)
                    record_copy.append(filename_out)
                    flo.release()

                    print('Remove temp')
                    os.remove(filename_tmp)

                    Ncopied += 1

                #############################################################
                # delete                                                    #
                #############################################################
                if not delete_file:
                    chrono.loop()
                    continue

                if not fli.attempt_to_lock():
                    print('Source locked -> continue')
                    chrono.loop()
                    continue

                input_pattern_str = hsl(patterns_in, njob, sep=sep)
                info = (''
                        + 'Record file (copy): ' + record_file_copy
                        + '\nInput pattern:      ' + input_pattern_str
                        + '\nAction:             ' + hsl(actions, 1)
                        + '\nMinimum file age:   ' + hsl(ages, 1)
                        )
                chrono.set_info(info).show()

                # delete file
                print('Delete source')
                os.remove(filename_in)

                # recursively delete directories if empty
                dirname_in = os.path.dirname(filename_in)
                delete_empty_dirs(dirname_in)

                chrono.loop()

                Ndeleted += 1

def copy_one_file(*args):
    pass

def delete_one_file(*args):
    pass

################################################################
# helpers                                                      #
################################################################
def get_size(filename):
    return os.path.getsize(filename)

def check_whether_to_copy(
        filename_in, filename_out, min_age_hours, record_copy,
        ):
    """Return a bool."""
    # insanely long time
    min_age_years = min_age_hours / 24 / 365.24
    if not min_age_years < 100:
        return False

    # input file non-existant
    if not os.path.isfile(filename_in):
        print("Source does not exists anymore -> don't copy")
        return False

    # check age (input)
    time_infile = get_time_file(filename_in)
    time_act = dt.datetime.now() - dt.timedelta(hours=min_age_hours)
    if time_infile is None:
        return False
    if time_infile > time_act:
        print("Source too young for copying")
        return False

    if not os.path.isfile(filename_out):
        return True

    if filename_out not in record_copy:
        return True
        
    # output file younger
    time_outfile = get_time_file(filename_out)
    if time_outfile == time_infile:
        print("Destination existing and same age as source")
        return False

    if time_outfile > time_infile:
        print("Destination existing and younger than source")
        return False

def check_whether_to_delete(filename_in, min_age_hours):
    """Return a bool."""
    # insanely long time
    min_age_years = min_age_hours / 24 / 365.24
    if not min_age_years < 100:
        return False

    # input file non-existant
    if not os.path.isfile(filename_in):
        return False

    # check age (input)
    time_infile = get_time_file(filename_in)
    time_act = dt.datetime.now() - dt.timedelta(hours=min_age_hours)
    if time_infile is None:
        print("Source does not exists anymore -> don't delete")
        return False
    if time_infile > time_act:
        print("Source too young for deletion")
        return False

    return True

if __name__ == '__main__':
    setup = command_line.get_setup(None)
    main(setup)
