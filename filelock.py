"""A file lock.
    
    Based on code by EF's repository [1], but with large modifications by AA.

    Authors
    -------
    EF : Evan Fosmark <me@evanfosmark.com>
    AA : Andreas Anhaeuser <andreas.anhaeuser@posteo.net>

    References
    ----------
    [1] https://github.com/dmfrey/FileLock/blob/master/setup.py
"""
import os
import time
from pathlib import Path
 
class FileLockException(Exception):
    pass

class FileAlreadyLockedException(FileLockException, OSError):
    pass
 
class FileLock(object):
    """A file lock for use in a `with` context."""
    def __init__(self, file_name, timeout=86400):
        """Prepare the file locker."""
        self.file_name = file_name
        self.timeout = timeout
        self.lockfile = self._get_lock_file_name()
        self._clean_if_timeout()
        self.created_lock = False

    def __enter__(self):
        """Called upon entering a with statement."""
        return self
 
    def __exit__(self, exc_type, exc_value, traceback):
        """Release the lock."""
        if self.created_lock:
            self.release()
 
    def __del__(self):
        """Make sure the lockfile isn't left lying around."""
        if not hasattr(self, 'created_lock'):
            return
        if self.created_lock:
            self.release()

    def _get_lock_file_name(self):
        """Return a str."""
        return '%s.lock' % self.file_name

    def _clean_if_timeout(self):
        """Remove lock file if older than `self.timeout`."""
        if not os.path.isfile(self.lockfile):
            return self

        mtime = os.path.getmtime(self.lockfile)
        now = time.time()
        age = now - mtime
        if age > self.timeout:
            os.remove(self.lockfile)
        
        return self

    def is_locked(self):
        """Return a bool."""
        return os.path.isfile(self.lockfile)

    def touch(self):
        Path(self.lockfile).touch()

    def attempt_to_lock(self):
        """Try to lock. Return True if successful, False otherwise.

            The method is immune to any OSError raised by trying to lock the
            file, but not to other exceptions.

            Returns
            -------
            success : bool
                True if:
                - is already lock by self
                - has not been locked by other process

                False if:
                - locked by other process
                - OSError while trying to lock
        """
        if self.is_locked():
            if self.created_lock:
                self.touch()
                return True
            else:
                return False

        # If this line is reached, it can be assumed to be unlocked
        try:
            self.lock()
            return True
        except OSError as exception:
            return False

    def lock(self):
        """Create lock file or raise exception if it already exists."""
        # create directory
        dirname = os.path.dirname(self.lockfile)
        if dirname == '':
            # current directory does always exist
            pass
        elif not os.path.isdir(dirname):
            os.makedirs(dirname)

        # make sure it's not already locked by someone else
        if self.is_locked() and not self.created_lock:
            raise FileAlreadyLockedException(
                    'Trying to lock a file that has already been locked'
                    + ' by another instance: %s' % self.file_name
                    )

        # create file
        with open(self.lockfile, 'w'):
            pass

        self.created_lock = True
        return self
 
    def release(self):
        """Release the lock by removing the lock file."""
        if not self.is_locked():
            return self

        if not self.created_lock:
            raise FileLockException(
                    'Trying to remove a lock that was created'
                    + ' by another instance: %s' % self.lockfile
                    )

        if os.path.isfile(self.lockfile):
            try:
                os.remove(self.lockfile)
            except FileNotFoundError:
                # This happens if file is removed by other process between
                # check and deletion.
                pass

        self.created_lock = False

        return self
