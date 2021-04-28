#!/usr/bin/python3
"""Testing suite for chronometer.Chronometer."""

# standard
import unittest
from time import sleep

# module to be tested
from misc.date_time.timer import Timer

class BasicFunctionality(unittest.TestCase):
    def setUp(self):
        print('----')
        self.title = 'No title.'

    def test_start_stop(self):
        self.title = 'Start and stop'
        timer = Timer(self.title)
        timer.start().show()
        timer.stop()

    def test_set(self):
        self.title = 'Set'
        timer = Timer(self.title)
        timer.set(2).show()

    def test_reset(self):
        self.title = 'Reset'
        timer = Timer(self.title)
        timer.start().stop().show().reset().show()

    def tearDown(self):
        print(self.title)
        print('----')

class ContextManager(unittest.TestCase):
    def setUp(self):
        Timer.autostart = True
        Timer.autoshow = True
        print('-----')

    def test_in_context(self):
        with Timer() as timer:
            print('Rush through context.')
        print('Sleep outside context')
        sleep(1)
        timer.show()

    def test_exit_function(self):
        with Timer() as timer:
            print('Sleep inside context.')
            sleep(1)

        print('Rush through outside')
        timer.show()

    def tearDown(self):
        Timer.autostart = False
        Timer.autoshow = False

if __name__ == '__main__':
    unittest.main()
