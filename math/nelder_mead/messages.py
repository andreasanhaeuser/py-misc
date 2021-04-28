#!/usr/bin/python
"""Interactive message display sub-module to nelder_mead."""

# standard modules
import datetime as dt

# local modules
from misc.text import string_utils as su

_col_widths = [25, 21, 21, 12]

def hum_fmt(x, digits=2):
    """Convert float to human-friendly string."""
    return su.human_format(x, digits=digits, mode='power')

def display_intermediate_message(record, options):
    """Display progress overview.

        Display something like this:

                        current (unc.)  | initial (unc.)  | tolerance
        function    :   4.9 (0.4)       | 12.1 (2.0)      | 0.01          |||         |
        parameter 0 :   4.54            | 8.64            | 0.1           ||||||||    |
        paramater 1 :
        paramater 2 :
        paramater 3 :
        parameter 4 :    -43e-4 (43e-6) | -0.1 (0.02)     | 0.01          |||||||||||||

        iterations  :     45 /  200     |||||||      |    ( 34s per iteration )
        func. eval. :     80 /  400     |||||||||    |    ( 18s per evaluation )
        time        :   2:43 / 5:00     |||||||      |    ( 2:17 remaining )
    """

    lines = []
    lines += '\n'
    lines += get_parameter_section(record, options)

    lines += '\n'
    lines += get_exit_limits_section(record, options)

    lines += '\n'
    lines += get_action_history_section(record, options)

    text = ''.join(lines)
    remove_old_message(record, options)
    record_message(record, options, text)
    print(text)
    return text

def remove_old_message(record, options):
    """Clear display from previously issued message."""
    if 'last_message' not in record.keys():
        record['last_message'] = ''
    text = record['last_message']
    Nlines = text.count('\n') + 1

    for nline in range(Nlines):
        print(su.move_cursor('up') + su._CLEARLINE + su.move_cursor('up'))

def record_message(record, options, text):
    """Remember latest issued message."""
    record['last_message'] = text

def get_exit_limits_section(record, options):
    """Return a list of str."""
    lines = []
    parameters = ('maxiter', 'maxfev', 'maxtime')
    for parameter in parameters:
        line = get_exit_limit_line(record, options, parameter)
        lines.append(line)
    return lines

def get_exit_limit_line(record, options, parameter):
    """Return a str."""
    secs_passed = record['time_passed'].total_seconds()

    # maxiter
    if parameter == 'maxiter':
        name = 'iterations'
        val = record['Niter']
        val_max = options['maxiter']
        val_str = '%i' % val
        val_max_str = str(val_max)
        rate = 1. * secs_passed / val
        rate_str = hum_fmt(rate)
        annotation = rate_str + 's / iter.'

    # maxfev
    elif parameter == 'maxfev':
        name = 'func. eval.'
        val = record['Nfev']
        val_max = options['maxfev']
        val_str = '%i' % val
        val_max_str = str(val_max)
        rate = 1. * secs_passed / val
        rate_str = hum_fmt(rate)
        annotation = rate_str + 's / eval.'

    # maxtime
    elif parameter == 'maxtime':
        name = 'time'
        passed = record['time_passed']
        maxtime = options['maxtime']
        val = secs_passed
        val_max = maxtime.total_seconds()
        val_str = nice_timedelta_string(passed)
        val_max_str = nice_timedelta_string(maxtime)
        if passed > maxtime:
            rem = dt.timedelta(0)
        else:
            rem = maxtime - passed
        rem_str = max(0, nice_timedelta_string(rem))
        annotation = rem_str + ' remaining'

    # column width
    cws = _col_widths

    # progress bars
    fraction = min(1, 1. * val / val_max)
    if fraction == 1:
        color = su.color('red')
    else:
        color = ''
    progress_bar = su.progress_bar(fraction, cws[2]+1, color)

    name = name[:cws[0]-1]

            # '{:>{hcw}s} / {:>{hcw}s}'.format(val_str, val_max_str, hcw=hcw) + \
    # build line
    line = (su.color('bold') + color +
            ('%s:' % name).ljust(cws[0]) + su.color(None) +
            color + '{:>{cw}s}'.format(
                ('%s / %s ' % (val_str, val_max_str)), cw=cws[1]-2) +
            ' %s' % progress_bar +
            ' ( %s )' % annotation +
            '\n'
            )
    return line

def get_parameter_section(record, options):
    """Return a list of str."""
    header = get_parameter_header_line()
    lines = [header]

    N = options['Nparam']

    for n in range(-1, N):
        line = get_parameter_line(record, options, n)
        lines.append(line)

    return lines

def get_parameter_header_line():
    """Return a str."""
    cws = _col_widths
    words = [su.color('bold')]
    words.append(''.ljust(cws[0]))
    words.append('initial (unc.) |'.rjust(cws[1]))
    words.append('current (unc.) |'.rjust(cws[2]))
    words.append('  tol.')
    words.append(su.color(None))
    line = ''.join(words) + '\n'
    return line

def get_parameter_line(record, options, n=-1):
    """Return a str.
    
        If n == -1 return line for forward function.
    """
    if n == -1:
        name = 'function'
        val = record['f']
        unc = record['f_unc']
        val0 = record['f0']
        unc0 = record['f0_unc']
        tol = options['fatol']
    else:
        name = record['parameter_names'][n]
        val = record['x'][n]
        unc = record['x_unc'][n]
        val0 = record['x0'][n]
        unc0 = record['x0_unc'][n]
        tol = options['xatol'][n]

    if unc <= tol:
        fraction = 1.
    else:
        fraction = tol / unc

    # column width
    cws = _col_widths

    name = name[:cws[0]-1]

    # progress bars
    if fraction == 1:
        color = su.color('g')
    else:
        color = ''
    progress_bar = su.progress_bar(fraction, cws[3], color)

    # build line
    line = (# Name:
            su.color('bold') + color + (name + ':').ljust(cws[0]) +
            # Color:
            su.color(None) + color +
            # Initial:
            u'{:>{cw}s} '.format(val_unc_str(val0, unc0) + ' |', cw=cws[1]) +
            # Current:
            u'{:>{cw}s} '.format(val_unc_str(val, unc) + ' |', cw=cws[2]-1) +
            # Progress bar:
            u'{:>{w}s} '.format(hum_fmt(tol), w=6) + ' %s' % progress_bar +
            # Newline:
            '\n'
            )
    return line

def val_unc_str(value, unc):
    """Return str."""
    val_str = hum_fmt(value, 4)
    unc_str = hum_fmt(unc, 2)
    return '%s (%s)' % (val_str, unc_str)

def nice_timedelta_string(d):
    secs = int(d.total_seconds())
    if secs < 600:
        i = 3
    elif secs < 3600:
        i = 2
    else:
        i = 0
    return str(dt.timedelta(seconds=secs))[i:]

def display_termination_message(record, options=None):
    message = 'Succes                  : %s\n' % record['success'] + \
              'Terminatation reason    : %s\n' % record['message'] + \
              'Result                  : %s\n' % record['x'] + \
              '- uncertainty           : %s\n' % record['x_unc'] + \
              'Function value (result) : %f\n' % record['f'] + \
              '- uncertainty           : %f\n' % record['f_unc'] + \
              'Ititial guess           : %s\n' % record['x0'] + \
              '- uncertainty           : %s\n' % record['x0_unc'] + \
              'Ititial function value  : %f\n' % record['f0'] + \
              'Iterations              : %i\n' % record['Niter'] + \
              'Function evaluations    : %i\n' % record['Nfev'] + \
              ''
    print(message)

def get_action_history_section(record, options):
    """Return a str."""
    # constants
    colors = [''] * 5
    endc = su.color(None)

    chars = 'ISCRE'
    symbols = (
            unichr(0x2584),
            unichr(0x2581),
            unichr(0x2582),
            unichr(0x2584),
            unichr(0x2588),
            )

    length = sum(_col_widths) + 2 * len(_col_widths)
    actions = record['action'][-length:]

    line_header = su.color('bold') + 'Simplex modification:\n' + endc
    line_footer = (
            '(%sI%s)nitialize' % (colors[0], endc) +
            ', (%sS%s)hrink' % (colors[1], endc) +
            ', (%sC%s)ontract' % (colors[2], endc) +
            ', (%sR%s)eflect' % (colors[3], endc) +
            ', (%sE%s)xpand' % (colors[4], endc) +
            su.color(None) + '\n'
            )

    line_char = ''
    line_sym = ''
    for action in actions:
        char = chars[action]
        sym = symbols[action]
        color = colors[action]
        line_char = line_char + color + char + endc
        line_sym = line_sym + color + sym + endc

    finish = lambda s: s + su.color(None) + '\n'

    line_char = finish(line_char)
    line_sym = finish(line_sym)

    return line_header + line_sym + line_char + line_footer

def get_termination_message(record, options):
    """Return a str."""
    # convergence
    if record['convergence_reached']:
        return 'Convergence reached.'

    # iteration limit
    itlim, message = record['it_lim_reached']
    if itlim:
        return message

    # other
    return 'Termination reason unknown.'
