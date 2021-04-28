#!/usr/bin/python3
"""Testing suite for string_converters."""

# standard
import unittest

# PyPI
import numpy as np

# module to be tested
from misc.text import string_converters as conv

class NumberConversion(unittest.TestCase):
    def setUp(self):
        self.words = ['300', 'nan', '-', 4]
        self.missing_words = ('nan', '-')

    def test_float(self):
        self.title = 'To float'
        self.dtype = float
        self.missing_value = np.nan
        self.expected = np.array([300., np.nan, np.nan, 4.])

    def test_int(self):
        self.title = 'To int'
        self.dtype = int
        self.missing_value = 0
        self.expected = np.array([300, 0, 0, 4])

    def test_scalar(self):
        self.words = '300'
        self.title = 'Scalar to int'
        self.dtype = int
        self.missing_value = 0
        self.expected = 300

    def tearDown(self):
        result = conv.string_to_number(
                self.words, self.dtype, missing_str=self.missing_words,
                missing_value=self.missing_value, 
                )
        print('*' * 40)
        print(self.title)
        print('Input:    %s' % repr(self.words))
        print('Type:     %s' % self.dtype)
        print('Expected: %s' % self.expected)
        print('Result:   %s' % result)


if __name__ == '__main__':
    unittest.main()
