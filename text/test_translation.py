#!/usr/bin/python3
"""Testing suite for translation module."""

# standard
import unittest

# to be tested
from misc.text.translation import translate

class Translation(unittest.TestCase):
    def setUp(self):
        self.filename = 'test_files/dictionary.txt'
        self.kwargs = {}

    def test_regular(self):
        self.term = 'Brot'
        self.expected = 'bread'

    def test_case_sensitive(self):
        self.term = 'brot'
        self.expected = 'bread'
        self.kwargs['ignore_case'] = True

    def test_list(self):
        self.term = 'Bank'
        self.expected = ['bench', 'bank']

    def test_commas(self):
        self.term = 'wieso, weshalb, warum'
        self.expected = ['why, why, why', 'Why']

    def tearDown(self):
        expected = self.expected
        result = translate(self.term, filename=self.filename, **self.kwargs)

        if not isinstance(result, list):
            self.assertEqual(result, expected)

        else:
            Nres = len(result)
            Nexp = len(expected)
            self.assertEqual(Nres, Nexp)
            for n in range(Nres):
                self.assertEqual(result[n], expected[n])


if __name__ == '__main__':
    unittest.main()
