#! /usr/bin/env python
"""Create an overview of python scripts and functions in a directory."""

import os
import sys

def scripts_in_dir(path):
    """Return a list of pairs (script name & docstring header)."""
    introwords = [' ', '\n', '#', 'import', 'from']

    listing = sorted(os.listdir(path))
    scripts = []
    for name in listing:
        ###################################################
        # RETRIEVE SCRIPT NAME                            #
        ###################################################
        longname = path + '/' + name
        if len(name) < 3:
            continue
        if name[0] == '.':
            continue
        if name[:2] == '__':
            continue
        if name[-3:] != '.py':
            continue
        if not os.path.isfile(longname):
            continue

        ###################################################
        # RETRIEVE DOCSTRING HEADER                       #
        ###################################################
        f = open(longname, 'r')
        lines = f.readlines()
        f.close()

        header = ''
        for line in lines:
            # skip intro lines
            skip = False
            for word in introwords:
                L = len(word)
                if line[:L] == word:
                    skip = True
                    break
            if skip:
                continue

            # check whether line is start of a docstring
            if line[:3] not in ['"""', "'''"]:
                break

            # retrieve header of this string
            header = line[3:]

            # delete trailing closing quotes
            if header[-3:] in ['"""', "'''"]:
                header = header[:-3]

        scripts.append([name, header])

    return scripts

def funcs_in_script(filename):
    """Return a list of pairs (function name & docstring header)."""
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    N = len(lines)
    funcs   = []
    for n in range(N):
        line = lines[n]

        ###################################################
        # RETRIEVE FUNCTION NAME                          #
        ###################################################
        if not line[:4] == 'def ':
            continue
        if not '(' in line:
            continue
        end = line.index('(')
        name = line[4:end]

        ###################################################
        # RETRIEVE DOCSTRING HEADER                       #
        ###################################################
        header = ''
        for m in range(n, N - 1):
            line = lines[m]

            # this should not happen (when coded in python syntax, a closing
            # parenthesis must appear first)
            if m > n and line[:4] == 'def ':
                break

            # this marks the end of the function definition
            if '):' in line:
                hline = lines[m + 1]    # potential docstring header line
                                        # if it exists, then here


                # remove leading white spaces:
                while hline[0] == ' ':
                    hline = hline[1:]

                # check whether it is in fact (the start of) a docstring
                if hline[:3] not in ['"""', "'''"]:
                    break

                # take the first line of this docstring
                header = hline[3:-1]

                # remove docstring closing:
                if header[-3:] in ['"""', "'''"]:
                    header = header[:-3]

                # ignore outdated functions if labelled as such:
                if header.lower()[:10] == '[outdated]':
                    name = None
                if header.lower()[:1] == '*':
                    name = None
                break

        if name is None:
            continue

        funcs.append([name, header])

    return funcs

def dirs_in_dir(path):
    """Return list of subdirectories not starting with '.' ."""
    listing = sorted(os.listdir(path))

    dirs = []
    for name in listing:
        longname = path + '/' + name
        if name[0] == '.':
            continue
        if not os.path.isdir(longname):
            continue
        dirs.append(name)

    return dirs

def text_for_funcs_in_script(filename, prefix):
    """Return a str (one line for each function)."""
    funcs = funcs_in_script(filename)

    ###################################################
    # FIND LENGTH OF LONGEST FUNCTION NAME            #
    ###################################################
    maxlen = 0
    for func in funcs:
        name, header = func
        length = len(name)
        if length > maxlen:
            maxlen = length

    ###################################################
    # CREATE ONE LINE FOR EACH FUNCTION               #
    ###################################################
    text = ''
    for func in funcs:
        name, header = func
        namep = name + '()'
        line = prefix + namep.ljust(maxlen + 3) + '> ' + header + '\n'
        text += line

    return text
    
def text_for_script(script, prefix, headerpos):
    """Return a str (one line exactly)."""
    filename, header = script
    name = filename[:-3]
    text = (prefix + name).ljust(headerpos - 1)
    text = text + ' > ' + header + '\n'
    return text
  
def text_for_file(script, filename, prefix, last=False):
    """Return a str (one line for script and each function)."""
    sprefix = prefix + '+-{S} '
    if last:
        fprefix = prefix + '   +-[F] '
    else:
        fprefix = prefix + '|  +-[F] '

    ###################################################
    # SCRIPT HEADER                                   #
    ###################################################
    ftext = text_for_funcs_in_script(filename, fprefix)
    if '>' in ftext:
        headerpos = ftext.index('>')
    else:
        headerpos = 0
    stext = text_for_script(script, sprefix, headerpos)
    text = stext + ftext

    return text

def text_for_directory(path, prefix='', last=False):
    """Return a str."""
    if path[-1] == '/':
        path = path[:-1]

    scripts = scripts_in_dir(path)
    dirs = dirs_in_dir(path)
    

    ###################################################
    # DIRECTORY LINE                                  #
    ###################################################
    text = prefix + '+-<D> ' + path + '/\n'

    if last:
        dprefix = '   '
    else:
        dprefix = '|  '

    sprefix = prefix + dprefix 
    ###################################################
    # SCRIPTS                                         #
    ###################################################
    N = len(scripts)
    n = 0
    for script in scripts:
        n += 1
        if n < N:
            last = False
        else:
            last = not(any(dirs))

        name, header = script

        filename = path + '/' + name
        lines = text_for_file(script, filename, sprefix, last=last)

        text += lines

        if not last:
            line = sprefix + '|\n'
        else:
            line = ''
        text += line

    ###################################################
    # DIRECTORIES                                     #
    ###################################################
    N = len(dirs)
    n = 0
    for name in dirs:
        n += 1
        last = n == N
        subdir = path + '/' + name
        lines = text_for_directory(path=subdir, prefix=dprefix, last=last)
        text += lines

        if not last:
            line = sprefix + '|\n'
        else:
            line = ''
        text += line

    return text

def print_text(text):
    """Print color-coded text on screen."""
    colors = [
        ['<D>', '\n', '01;31'],   # directory
        ['{S}', '\n',  '01;33'],   # script
        ['[F]', '>',  '01;34'],   # function
#        [' >',  '\n', '10;39'],   # header
        ]
    ncol = '00;31'

    CSI = '\x1B['

    lines = text.split('\n')
    newtext = ''
    N = len(lines)
    for n in range(N):
        for color in colors:
            line = lines[n]
            start, end, col = color

            if not start in line:
                continue

            L = len(start)
            ibeg = line.index(start)

            if end in line[ibeg:]:
                iend = line.index(end, ibeg)
            else:
                iend = None

            if iend is not None:
                first  = line[:ibeg]
                middle = line[ibeg:iend]
                last   = line[iend:]
            else:
                first  = line[:ibeg]
                middle = line[ibeg:]
                last   = ''

            ins1 = CSI + col + 'm'
            ins2 = CSI + '0m'
            lines[n] = first + ins1 + middle + ins2 + last

        if any(lines[n]):
            newtext += (lines[n] + '\n')
    
    print(newtext)
    return newtext


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) > 0:
        base = args[0]
    else:
        base = '.'

    if len(args) > 1:
        out = args[1]
    else:
        out = None
    

    if not os.path.isdir(base):
        print('Could not access directory: ' + base)

    else:
        text = text_for_directory(base, last=True)
        if not out:
            print_text(text)
        else:
            if os.path.isfile(out):
                os.remove(out)
            f = open(out, 'w')
            f.write(text)
            f.close()










    
    







