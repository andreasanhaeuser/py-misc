#!/usr/bin/python3
"""Use a simple text file as a record.

    The module provide very basic functionality:
    - add a record to the file
    - check whether a record exists in the file
    - remove the record from the file
"""

# if running in python2
from __future__ import print_function

# standard modules
import os

class Record(object):
    """A class that facilitates use of the basic functions.
    
        Parameters
        ----------
        filename : str
            record file
        create_multiple_appearances : bool, optional
            (default: False) add entry another time if it is already in the
            file.
        remove_multiple_appearances : bool, optional
            (default: True) remove all appearances of entry if it exists
            multiple times
    """
    def __init__(
            self, filename, create_multiple_appearances=False,
            remove_multiple_appearances=None,
            ):
        if not isinstance(filename, str):
            raise TypeError('filename must be str.')
        if not isinstance(create_multiple_appearances, bool):
            raise TypeError('create_multiple_appearances must be str.')
        if remove_multiple_appearances is None:
            remove_multiple_appearances = not create_multiple_appearances
        if not isinstance(remove_multiple_appearances, bool):
            raise TypeError('remove_multiple_appearances must be str.')

        self.filename = os.path.expanduser(filename)
        self.create_multiple_appearances = create_multiple_appearances
        self.remove_multiple_appearances = remove_multiple_appearances
        initialize(filename, overwrite=False)

    def __str__(self):
        return '%i entries in %s' % (self.count_all(), self.filename)

    def __repr__(self):
        return 'Record with file %s' % self.filename

    def __len__(self):
        return self.count_all()

    def __iter__(self):
        return self.get_all().__iter__()

    def __contains__(self, entry):
        return self.contains(entry)

    def add(self, entry):
        """Add entry to record and return self."""
        add(entry, self.filename, self.create_multiple_appearances)
        return self

    def append(self, entry):
        return self.add(entry)

    def contains(self, entry):
        """Return a bool."""
        return contains(entry, self.filename)

    def count(self, entry):
        """Return number of appearances in record as int."""
        return count_appearances(entry, self.filename)

    def count_all(self):
        """Return number of entries (may contain duplicates)."""
        return count_all_entries(self.filename)

    def count_unique_entries(self):
        raise NotImplementedError()

    def erase(self):
        """Delete the record file if existant."""
        if os.path.isfile(self.filename):
            os.remove(self.filename)
        return self

    def get_unique_entries(self):
        raise NotImplementedError()

    def get_all(self):
        # file does not exist --> False
        if not os.path.isfile(self.filename):
            return []

        # read file
        with open(self.filename, 'r') as fid:
            lines = fid.readlines()

        return [line.strip() for line in lines]

    def pop(self, index=-1):
        """Remove and return item at index (default last).

            Parameters
            ----------
            index : int, optional
                (default: -1)
            
            Returns
            -------
            str : element at index position

            Raises
            ------
            IndexError
                record is empty or index is out of range.
        """
        # file does not exist
        if not os.path.isfile(self.filename):
            raise IndexError('pop from empty record')

        # read file
        with open(self.filename, 'r') as fid:
            lines = fid.readlines()

        # file empty
        if not any(lines):
            raise IndexError('pop from empty record')

        # index out of range
        L = len(lines)
        if not -L <= index < L:
            raise IndexError('pop index out of range')

        # regular case
        entry = lines[index]
        self.remove(entry)
        return entry

    def remove(self, entry):
        """Remove entry from record and return self."""
        remove(entry, self.filename, self.remove_multiple_appearances)
        return self

    def sort(self):
        # file does not exist --> False
        if not os.path.isfile(self.filename):
            return []

        # read file
        with open(self.filename, 'r') as fid:
            lines = fid.readlines()

        # sort
        lines_sorted = sorted(lines)

        # write file
        with open(self.filename, 'w') as fid:
            fid.writelines(lines_sorted)

        return self


###################################################
# USER FUNCTIONS                                  #
###################################################
def add(entry, filename, create_multiple_appearances=False):
    """Append a line to the record file.

        Creates the record file autonomously if it does not exist.

        Parameters
        ----------
        entry : str
            record to be added
        filename : str
            record file
        create_multiple_appearances : bool, optional
            (default: False) add entry another time if it is already in the
            file.

        Returns
        -------
        added : bool
            True if `entry` was added, False otherwise
    """
    # create if not already existing
    initialize(filename, overwrite=False)

    if not create_multiple_appearances:
        if contains(entry, filename):
            return False

    # normalize line
    line = str(entry).strip()

    # add it to file
    with open(filename, 'a') as fid:
        print(line, file=fid)

    return True

def contains(entry, filename):
    """Return a bool.

        Parameters
        ----------
        entry : str
            search item
        filename : str
            record file

        Returns
        -------
        found : bool
            True if `entry` is in file, False otherwise
    """
    nline = find_first_appearance(entry, filename)
    if nline is None:
        return False
    else:
        return True

def remove(entry, filename, remove_multiple_appearances=True):
    """Remove record from file.

        No exception is thrown if the entry does not exists.

        Parameters
        ----------
        entry : str
            search item
        filename : str
            record file
        remove_multiple_appearances : bool, optional
            (default: True) remove all appearances of entry if it exists
            multiple times

        Returns
        -------
        found : bool
            True if `entry` was found, False otherwise.
    """
    found = False
    repeat = True

    while repeat:
        found_once = remove_once(entry, filename)

        if found_once:
            found = True

        if not remove_multiple_appearances:
            repeat = False
        elif not contains(entry, filename):
            repeat = False

    return found

###################################################
# HELPERS                                         #
###################################################
def remove_once(entry, filename):
    """Remove first appearance of record from file.

        No exception is thrown if the entry does not exists.

        Parameters
        ----------
        entry : str
            search item
        filename : str
            record file

        Returns
        -------
        found : bool
            True if `entry` was found, False otherwise.
    """
    nline = find_first_appearance(entry, filename)
    if nline is None:
        return False

    # read file
    with open(filename, 'r') as fid:
        lines = fid.readlines()

    # remove one line
    reduced_lines = lines[:nline] + lines[nline+1:]

    # re-write file
    with open(filename, 'w') as fid:
        fid.writelines(reduced_lines)

    return True

def find_first_appearance(entry, filename):
    """Return line number of first appearance.
        
        Python counting (0 == first line).

        Parameters
        ----------
        entry : str
            search item
        filename : str
            record file

        Returns
        -------
        nline : int or None
            int: line number if the entry was found
            None: file does not exist or entry not found
    """
    # file does not exist --> False
    if not os.path.isfile(filename):
        return None

    # read file
    with open(filename, 'r') as fid:
        lines = fid.readlines()

    # search contents
    hit = False
    entry = str(entry).strip()
    for nline, line in enumerate(lines):
        line = line.strip()
        if entry == line:
            hit = True
            break

    if hit:
        return nline
    else:
        return None

def count_all_entries(filename):
    """Return an int.
        
        Parameters
        ----------
        filename : str
            record file

        Returns
        -------
        count : int or None
            number of records (may be duplicate)
    """
    # file does not exist --> False
    if not os.path.isfile(filename):
        return 0

    # read file
    with open(filename, 'r') as fid:
        lines = fid.readlines()

    return len(lines)

def count_appearances(entry, filename):
    """Return an int.
        
        Parameters
        ----------
        entry : str
            search item
        filename : str
            record file

        Returns
        -------
        count : int
            number of appearances in file
    """
    # file does not exist --> False
    if not os.path.isfile(filename):
        return 0

    # read file
    with open(filename, 'r') as fid:
        lines = fid.readlines()

    # search contents
    entry = str(entry).strip()
    count = 0
    for nline, line in enumerate(lines):
        line = line.strip()
        if entry == line:
            count += 1

    return count

def initialize(filename, overwrite=True):
    """Create record file.

        Parameters
        ----------
        filename : str
            record file
        overwrite : bool, optional
            (default: True) overwrite if file alreay exists

        Returns
        -------
        None
    """
    # remove already existing file
    if os.path.isfile(filename) and overwrite:
        os.remove(filename)

    # create directory
    directory = os.path.realpath(os.path.dirname(filename))
    if not os.path.isdir(directory):
        os.makedirs(directory)

    # create file
    if not os.path.isfile(filename):
        with open(filename, 'w'):
            pass

    return None


###################################################
# TESTS                                           #
###################################################
def test_write_and_look_up():
    import datetime as dt
    filename = './some/strange/path/record_file.tmp.txt'
    entries = ['first', 2, dt.datetime.now(), ('a', 'tuple')]

    print('record file: %s' % filename)
    for entry in entries:
        print('*' * 30)

        # look up
        found = contains(entry, filename)
        print('Entry %s in record file: %s' % (entry, found))

        # add
        print('Add it.')
        add(entry, filename)

        # look up again
        found = contains(entry, filename)
        print('Entry %s in record file: %s' % (entry, found))

    print('*' * 30)
    for entry in entries:
        # look up
        found = contains(entry, filename)
        print('Entry %s in record file: %s' % (entry, found))

def test_class_basic():
    import datetime as dt
    filename = '.record_file.tmp.txt'
    entries = ['first', 2, dt.datetime.now(), ('a', 'tuple')]
    record = Record(filename)

    print('record file: %s' % filename)
    for entry in entries:
        print('*' * 30)

        # look up
        found = record.contains(entry)
        print('Entry %s in record file: %s' % (entry, found))

        # add
        print('Add it.')
        record.add(entry)

        # look up again
        found = record.contains(entry)
        print('Entry %s in record file: %s' % (entry, found))

    print('*' * 30)
    for entry in entries:
        # look up
        found = record.contains(entry)
        print('Entry %s in record file: %s' % (entry, found))

    record.erase()


if __name__ == '__main__':
    test_class_basic()



