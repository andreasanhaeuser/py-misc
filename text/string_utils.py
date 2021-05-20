#!/usr/bin/python
"""String utils.

    Author
    ------
    Andreas Anhaeuser (AA) <andreas.anhaeuser@posteo.net>
"""

# standard modules
import sys
from collections import Iterable

# PyPI modules
import numpy as np

# misc
# from misc.math import round as round_utils

###################################################
# COLORS                                          #
###################################################
# foreground colors
_ENDC = '\033[0m'
_BOLD = '\033[1m'
_BLUE = '\033[94m'
_CYAN = '\033[96m'
_GREEN = '\033[92m'
_MAGENTA = '\033[95m'
_ORANGE = '\033[93m'
_RED = '\033[91m'
_YELLOW = '\033[93m'

# background colors
_BG_BLACK = '\033[40m'
_BG_BLUE = '\033[44m'
_BG_CYAN = '\033[46m'
_BG_RED = '\033[41m'
_BG_GREEN = '\033[42m'
_BG_ORANGE = '\033[43m'
_BG_MAGENTA = '\033[45m'
_BG_GREY = '\033[47m'

# motions
_UPWARD = '\033[1A'
_DOWNWARD = '\033[1B'
_FORWARD = '\033[1C'
_BACKWARD = '\033[1D'

_CLRSCR = '\033[2J'     # clear screen, move to (0, 0)
_CLEARLINE  = '\033[K'  # erase to end of line
_SAVECRS = '\033[s'     # save cursor position
_RESTORECRS = '\033[u'  # restore cursor position

###################################################
# HUMAN READABLE FORMATS                          #
###################################################
def human_format(num, digits=2, sep='', mode='prefix', type='float'):
    """Convert large and small numbers to human with SI prefixes.

        Parameters
        ----------
        num : float
            the number to convert
        digits : int, optional
            total (minimum) number of decimal digits to include. Default: 2
        sep : str, optional
            separator between the number and the prefix. Default: ''
        mode : {'prefix', 'power', 'tex'}, optional
            'prefix' : SI prefix
            'power' : power of 10, in regular format [ '1.4e5' ]
            'tex' : power of 10, in tex format [ '$\\times 10^{power}$' ]
            Default: 'prefix'
        type : {'float', 'int'}, optional
            if 'int', small numbers don't have decimal digits (i. e. '4' rather
            than '4.0')

        Returns
        -------
        str

        Example
        -------
        >>> human_format(234567, 1, ' ')
        '234.6 k'
    """
    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    valid_types = ['float', 'int']

    assert isinstance(digits, int)
    assert isinstance(sep, str)
    assert type in valid_types

    ###################################################
    # SPECIAL CASES                                   #
    ###################################################
    if np.isnan(num):
        return 'NaN'

    if num == 0:
        return '0'

    if np.isposinf(num):
        if mode == 'tex':
            return '$\\infty$'
        else:
            return 'inf'

    if np.isneginf(num):
        if mode == 'tex':
            return '$-\\infty$'
        else:
            return '-inf'

    ###################################################
    # SETUP                                           #
    ###################################################
    if mode == 'prefix':
        P = ['', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']
        p = ['', 'm', u'\u03bc', 'n', 'p', 'f', 'a', 'z', 'y']
    elif mode == 'power':
        powers = range(0, 30, 3)
        P = [''] + ['e%1.0f' % power for power in powers[1:]]
        p = [''] + ['e-%1.0f' % power for power in powers[1:]]
    elif mode == 'tex':
        powers = range(0, 30, 3)
        P = ['\\times\, 10^{%1.0f}' % power for power in powers]
        p = ['\\times\, 10^{-%1.0f}' % power for power in powers]
    elif mode in ('short scale', 'short_scale'):
        P = [
                '', 'thousand', 'million', 'billion', 'trillion',
                'quadrillion', 'quintillion', 'sextillion', 'septillion'
                ]
        p = [
                '', 'thousandth', 'millionth', 'billionth', 'trillionth',
                'quadrillionth', 'quintillionth', 'sextillionth',
                'septillionth'
                ]
    else:
        raise ValueError('Unknown mode: %s' % mode)


    ###################################################
    # MAGNITUDE                                       #
    ###################################################
    mag = 0

    # large numbers
    while abs(num) >= 1000:
        if mag >= len(P) - 1:
            break
        mag += 1
        num /= 1000.

    # small numbers
    while abs(num) < 1:
        if abs(mag) >= len(p) - 1:
            break
        mag -= 1
        num *= 1000

    ###################################################
    # NUMBER OF DIGITS                                #
    ###################################################
    # before '.'
    if abs(num) >= 100:
        n_before = 3
    elif abs(num) >= 10:
        n_before = 2
    else:
        n_before = 1

    # after '.'
    n_after = max(0, digits - n_before)

    ###################################################
    # FORMATTER                                       #
    ###################################################
    if mag == 0 and type == 'int':
        # special case: int close to 1:
        fmt = '%1.0f' + sep + '%s'
    else:
        # regular case:
        fmt = '%.' + ('%1.0f' % n_after) + 'f' + sep + '%s'

    ###################################################
    # MATH MODE                                       #
    ###################################################
    if mode == 'tex':
        fmt = '$' + fmt + '$'

    ###################################################
    # BUILD STRING                                    #
    ###################################################
    # select prefix
    if mag >= 0:
        prefix = P[mag]
    else:
        prefix = p[-mag]

    return fmt % (num, prefix)

def human_format_time(seconds):
    """Return N.Mu or NNu or NNNu."""
    value = seconds
    units = 's'
    if value > 300:
        value /= 60
        units = 'm'
        if value > 300:
            value /= 60
            units = 'h'
            if value > 100:
                value /= 24
                units = 'd'

    if value >= 10:
        fmt = '%1.0f'
    else:
        fmt = '%1.1f'

    return fmt % value + units

def percentage_string(fraction, sep=' ', sign=False):
    """Return a str with a sensible number of digits.

        The closer to 0 or 100%, the more the number of digits is increased.

        Parameters
        ----------
        fraction : float
            usually, but not necessarily, between 0 and 1
        sep : str
            separator between number and '%' sign in the output
        sign : bool, optional
            Preceed positive numbers with '+'. Default: False

        Returns
        -------
        str
            something like '0.024 %', '1.2 %', '63 %' or '99.9936 %'
    """
    if fraction < 0:
        return '-' + percentage_string(-fraction, sep=sep, sign=False)

    if fraction == 0:
        return '0%%'

    if fraction == np.inf:
        return 'infinite'

    if not np.isfinite(fraction):
        return 'NaN'

    pc = fraction * 100

    # > 100% : zero digits after .
    if pc >= 100:
        Ndig = 0

    # regular cases
    elif pc < 10:
        Ndig = 1 - np.floor(np.log10(pc))
    elif pc > 90:
        Ndig = 1 - np.floor(np.log10(100 - pc))
    else:
        Ndig = 0

    Ndig = min(6, int(Ndig))

    fmt = '%1.' + str(Ndig) + 'f' + sep + '%%'
    text = (fmt % pc)

    # sign
    if sign:
        text = '+' + text

    return text


def ordinal_str(number):
    """Return '2nd'/'14th'/'0th', etc."""
    assert isinstance(number, int)
    number = int(number)
    s = str(number)
    last_digit = s[-1]
    if s[-2:-1] == '1':
        suffix = 'th'
    elif last_digit == '1':
        suffix = 'st'
    elif last_digit == '2':
        suffix = 'nd'
    elif last_digit == '3':
        suffix = 'rd'
    else:
        suffix = 'th'
    return s + suffix

def group_thousands(number, decimals=0, round_digits=False, sep=','):
    """Return a str, where exponentials of 1000 are grouped.

        Parameters
        ----------
        number : float or int or Iterable of such
            the number to be grouped
        decimals : int >= 0, optional
            Default: 0. Number of decimal places to display.
        round_digits : bool or int, optional
            Default: False. Number of significant (leading) digits to round.
        sep : str, optional
            Default: ','. Thousands separator

        Returns
        -------
        str or list of such
            e. g. 12,345,678.9
            list of str if `number` is Iterable
    """
    # str -> number
    # =====================================================
    if isinstance(number, str):
        number = float(number)
    # =====================================================

    # `numbers` is iterable
    # =====================================================
    # recursively call function
    if isinstance(number, Iterable):
        results = []
        for num in number:
            result = group_thousands(num, decimals, round_digits, sep)
            results.append(result)
        return results
    # =====================================================

    # Special cases: non-numeric values
    # =====================================================
    if np.isnan(number):
        return 'NaN'
    if number == np.inf:
        return 'infinite'
    if number == -np.inf:
        return 'neg. infinite'
    # =====================================================

    assert decimals >= 0

    # Negative
    # =====================================================
    if number < 0:
        text = group_thousands(-number, decimals, sep)
        return '-' + text
    # =====================================================

    # Fractional part
    # =====================================================
    if decimals == 0:
        fractional_str = ''
    else:
        fractional = number % 1
        fmt = '%1.' + str(int(decimals)) + 'f'
        fractional_str = (fmt % fractional).lstrip('0')
    # =====================================================

    # Round
    # =====================================================
    if round_digits:
        number = round_utils.round_digits(number, round_digits)
    # =====================================================

    # Integer part
    # =====================================================
    if int(number) == 0:
        str_int = '0'
    else:
        factor = 1000
        position = 1
        str_ints = []
        while number >= position:
            value = number // position % factor
            str_ints.append('%03i' % value)
            position *= factor
        str_int = sep.join(str_ints[::-1])

    if number >= 1:
        str_int = str_int.lstrip('0')
    # =====================================================

    # combine int and fractional parts
    text = str_int + fractional_str

    return text

def leading_digits(number, ndigits=1, base=10, fmt='int'):
    """Return as int."""
    if base != 10:
        message = 'Only base 10 implemented so far, got %s' % str(base)
        raise NotImplementedError(message)

    if number < 0:
        return get_leading_digits(-number, ndigits, base)

    exp = int(np.floor(np.log10(number)))
    leaders = number // 10**(exp-(ndigits-1))
    return int(leaders)

###################################################
# COLOR FORMATTING                                #
###################################################
def highlighted_string_list(
        words, n=None, sep='|', start_hl=_BOLD, end_hl=_ENDC,
        index_error=False,
        ):
    """Return a string with one element highlighted.

        Parameters
        ----------
        words : list of str
        n : int or None
            element to be highlighted. If None, no element will be highlighted.
        sep : str, optional
            separator between list elements
        start_hl : str, optional
            highlight string
        end_hl : str, optional
            un-highlight string
        index_error : bool, optional
            (default: False) raise error if n out of bounds

        Returns
        -------
        str
    """
    # cast to str
    words = [str(word) for word in words]

    # highlight one elemnt
    if n is not None:
        if not isinstance(n, int):
            raise TypeError('n must be int.')
        try:
            words[n] = start_hl + words[n] + end_hl
        except IndexError as e:
            if index_error:
                raise e

    # join list to str
    line = sep.join(words)

    return line

###################################################
# PROGRESS BAR                                    #
###################################################
def progress_bar(
        fraction, length=79, fillcolor='', bgcolor='', delim_color=_BOLD,
        use_color=True,
        ):
    """Return a nice progress bar as str.

        Parameters
        ----------
        fraction : float
            fraction of bar to be filled, function optimized for
            0 <= fraction <= 1, but works with other values, too.
            Must be finite.
        length : int > 2
            full length of the progress bar including bounds
        fillcolor : str, optional
        bgcolor : str, optional
        delim_color : str, optional
        use_color : bool, optional

        Returns
        -------
        str
    """
    # input check
    # ===========================================================
    if not np.isfinite(fraction):
        raise ValueError('Need finite value, got %s' % fraction)

    if not isinstance(length, int):
        raise TypeError('length must be int, got %s' % type(length))
    if length < 3:
        raise ValueError('length must be larger than 2, got %s' % length)
    # ===========================================================

    if use_color:
        endcolor = _ENDC
    else:
        endcolor = ''
        fillcolor = ''
        bgcolor = ''
        delim_color=''

    beg = delim_color + u'\u2595' + endcolor
    end = delim_color + u'\u258f' + endcolor

    Nbar = length - 2

    # fraction > 1
    if fraction > 1:
        left_inner = progress_bar_inner(
                fraction=1, length=Nbar, fillcolor=fillcolor, bgcolor=bgcolor,
                use_color=use_color,
                )
        right_inner = progress_bar_inner(
                fraction=fraction-1, length=Nbar, fillcolor=_YELLOW,
                bgcolor=bgcolor, use_color=use_color,
                )
        line = beg + left_inner + right_inner

    # fraction < 0
    elif fraction < 0:
        return progress_bar(
                -fraction, length=length, fillcolor=_YELLOW,
                use_color=use_color,
                )

    # 0 <= fraction <= 1 (regular case)
    else:
        filling = progress_bar_inner(
                fraction, length=Nbar, fillcolor=fillcolor, bgcolor=bgcolor,
                use_color=use_color,
                )
        line = beg + filling + end

    return line

def progress_bar_inner(
        fraction, length=77, fillcolor='', bgcolor='', use_color=True,
        ):
    """*Helper function."""
    # input check
    # =========================================================
    if not isinstance(length, int):
        raise TypeError('length must be int, got %s' % type(length))
    if length < 1:
        raise ValueError('length must be positive, got %s' % length)
    # =========================================================


    # N : integers
    # F : float
    # full : number of (entirely) filled columns
    # blank : number of (entirely) clear columns
    Nbar = length

    if not use_color:
        fillcolor = ''
        bgcolor = ''
        endcolor = ''
    else:
        endcolor = _ENDC

    Ffull = fraction * Nbar
    Nfull = int(Ffull)
    Nblank = Nbar - Nfull - 1
    fractional_part = Ffull - Nfull
    Neighths = int(round(fractional_part * 8))

    block_full = get_frac_block(8)
    block_frac = get_frac_block(Neighths)

    # special case
    if Nfull == Nbar and fraction <= 1:
        block_frac = ''

    filling = bgcolor + fillcolor + \
            Nfull * block_full + \
            block_frac + \
            Nblank * ' ' + \
            endcolor
    return filling

def get_frac_block(Neighths):
    """Return a unicode block element as str."""
    # ========== INPUT CHECK  ============================ #
    if not isinstance(Neighths, int):
        raise TypeError('Neighths must be int.')
    if not 0 <= Neighths <= 8:
        raise ValueError('Neighths must be between 0 and 8.')

    # ========== SPECIAL CASE: EMPTY BLOCK ELEMENT  ====== #
    if Neighths == 0:
        return ' '

    # ========== COMPUTE UNICODE POSITION  =============== #
    # starting point: completely filled block
    pos0 = 0x2588

    # increase position for each missing eighth
    N_eighths_missing = 8 - Neighths
    pos = pos0 + N_eighths_missing

    # ========== CONVERT TO UNICODE CHARACTER =============================== #
    #
    # The function that does this has different names in python2 and python3.
    # --> Select the appropriate function and call it `convert_to_unicode`.

    # major version
    python_version = sys.version_info[0]    # (int)
    if python_version <= 2:
        convert_to_unicode = unichr
    else:
        convert_to_unicode = chr

    character = convert_to_unicode(pos)
    # ======================================================================= #

    return character

def get_frac_block_vertical(Neighths):
    """Return a unicode block element as str."""
    # input check
    if not isinstance(Neighths, int):
        raise TypeError('Neighths must be int.')
    if not 0 <= Neighths <= 8:
        raise ValueError('Neighths must be between 0 and 8.')

    # special case: empty block element
    if Neighths == 0:
        return ' '

    # regular case
    pos0 = 0x2580 + int(Neighths)
    return unichr(pos)

###################################################
# COLORS                                          #
###################################################
def color(col=None):
    if col is None:
        col = 'none'
    if not isinstance(col, str):
        raise TypeError('col must be str or None.')
    col = col.lower()

    if col == 'none':
        return _ENDC
    elif col in ('b', 'blue'):
        return _BLUE
    elif col in ('c', 'cyan'):
        return _CYAN
    elif col in ('g', 'green'):
        return _GREEN
    elif col in ('o', 'orange'):
        return _ORANGE
    elif col in ('m', 'magenta'):
        return _MAGENTA
    elif col in ('r', 'red'):
        return _RED
    elif col in ('w', 'white', 'bold'):
        return _BOLD
    elif col in ('y', 'yellow'):
        return _YELLOW
    else:
        raise ValueError('Unknown color: %s' % col)

def bg_color(col):
    if col is None:
        col = 'none'
    if not isinstance(col, str):
        raise TypeError('col must be str or None.')
    col = col.lower()

    if col == 'none':
        return _BG_BLACK
    elif col in ('b', 'blue'):
        return _BG_BLUE
    elif col in ('c', 'cyan'):
        return _BG_CYAN
    elif col in ('g', 'green'):
        return _BG_GREEN
    elif col in ('m', 'magenta'):
        return _BG_MAGENTA
    elif col in ('o', 'orange'):
        return _BG_ORANGE
    elif col in ('r', 'red'):
        return _BG_RED
    elif col in ('w', 'white', 'bold'):
        return _BG_BOLD
    elif col in ('y', 'yellow'):
        return _BG_YELLOW
    else:
        raise ValueError('Unknown color: %s' % col)

###################################################
# CURSOR                                          #
###################################################
def move_cursor(m):
    """u(pward), d(ownward), f(orward) or b(ackward)."""
    if not isinstance(m, str):
        raise TypeError('m must be str or None.')
    m = m.lower()[0]

    if m == 'u':
        return _UPWARD
    elif m == 'd':
        return _DOWNWARD
    elif m == 'f':
        return _FORWARD
    elif m == 'b':
        return _BACKWARD
    else:
        raise ValueError('Unknown motion: %s' % col)

def save_cursor_pos():
    """Return a str."""
    return _SAVECRS

def restore_cursor_pos():
    """Return a str."""
    return _RESTORECRS

###################################################
# CLEAR                                           #
###################################################
def clear_screen():
    """Return a str."""
    return _CLRSCR

def clear_line():
    """Return a str."""
    return _CLEARLINE

################################################################
# crop                                                         #
################################################################
def crop(s, word):
    """Remove leading and trailing appearences of `word`."""
    s = lcrop(s, word)
    s = rcrop(s, word)
    return s

def lcrop(s, word):
    """Remove leading appearences of `word`."""
    if s == '':
        return s

    L = len(word)
    while s.startswith(word):
        s = s[L:]
    return s

def rcrop(s, word):
    """Remove trailing appearences of `word`."""
    if s == '':
        return s

    L = len(word)
    while s.endswith(word):
        s = s[:-L]
    return s
