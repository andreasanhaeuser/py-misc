#!/usr/bin/python3
"""Testing suite for filelock.FileLock."""

# standard modules
import os
import unittest

# local filelock module
from misc.filelock import FileLock

class RegularUsage(unittest.TestCase):
    def setUp(self):
        print('-----')

    def test_in_context(self):
        print('Context management')
        filename = 'testfile_to_be_locked.tmp'
        lockfile = 'testfile_to_be_locked.tmp.lock'
        if os.path.isfile(lockfile):
            os.remove(lockfile)

        # check empty `with` context
        with FileLock(filename) as lock:
            pass
        self.assertFalse(os.path.isfile(lockfile))

        # check that lock is cleaned up automatically
        with FileLock(filename) as lock:
            self.assertFalse(os.path.isfile(lockfile))
            lock.lock()
            self.assertTrue(os.path.isfile(lockfile))
        self.assertFalse(os.path.isfile(lockfile))

        # check that lock is created and then cleaned up
        with FileLock(filename) as lock:
            self.assertFalse(os.path.isfile(lockfile))
            lock.lock()
            self.assertTrue(os.path.isfile(lockfile))

            # lock it again
            lock.lock()
            self.assertTrue(os.path.isfile(lockfile))

            lock.release()
            self.assertFalse(os.path.isfile(lockfile))

        self.assertFalse(os.path.isfile(lockfile))

        # try to access an already locked file
        lock_other = FileLock(filename)
        lock_other.lock()
        self.assertTrue(os.path.isfile(lockfile))
        self.assertTrue(lock.is_locked())

        # check that lock realises the conflict
        with FileLock(filename) as lock:
            print('Some output text inside the FileLock context.')
            self.assertTrue(os.path.isfile(lock.lockfile))
            self.assertTrue(os.path.isfile(lockfile))
            self.assertTrue(lock.is_locked())
            self.assertRaises(Exception, lock.lock)
       
        self.assertTrue(os.path.isfile(lockfile))
        del lock_other
        self.assertFalse(os.path.isfile(lockfile))


if __name__ == '__main__':
    unittest.main()
