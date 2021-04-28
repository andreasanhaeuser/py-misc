#!/usr/bin/python3
"""Testing suite for tables."""

# standard modules
import unittest
from collections import Iterable

# PyPI modules
import numpy as np

# to be tested
import misc.in_out.tables as tab

class NameList(unittest.TestCase):
    def setUp(self):
        pass

    def test_recursive(self):
        filename = 'test_files/namelist_child.txt'
        result = tab.read_namelist(filename, convert_to_number=True)
        expected = {
                'parameter_a' : 'default_by_parent_1_a',
                'parameter_b' : ['default_by_parent_1', 'b'],
                'parameter_c' : 'default_by_parent_1_c',
                'parameter_d' : 'default_by_parent_2_d',
                'parameter_e' : ['default_by_parent_2', 'e'],
                'parameter_f' : 'default_by_parent_2_f',
                'parameter_g' : [4, 3], 
                'parameter_h' : 'a string', 
                'parameter_i' : ['two', 'strings'], 
                'parameter_j' : 1000, 
                'parameter_k' : 'default_by_grand_parent_k',
                'parameter_l' : 'default_by_grand_parent_l',
                }

        self.title = 'Recursion'
        self.result = result
        self.expected = expected


    def tearDown(self):
        print('')
        print('*' * 40)
        print(self.title)
        print('')
        result = self.result
        expected = self.expected
        for key in expected:
            value = expected[key]
            if not isinstance(value, str) and isinstance(value, Iterable):
                for n, val in enumerate(value):
                    self.assertEqual(value[n], result[key][n])
            self.assertEqual(value, result[key])

        print('Loaded_files:')
        for fn in result['loaded_files']:
            print(fn)


class ColumnList(unittest.TestCase):
    def setUp(self):
        pass

    def test_table(self):
        filename = 'tables__test_file_column_list.txt'
        sep = '|'
        convert_to_number = True

        data = tab.read_column_list_with_headers(
                filename, sep=sep, convert_to_number=convert_to_number,
                missing_str='-',
                )

        for key in sorted(data):
            print('%s: %s' % (key.ljust(32), str(data[key])))

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
