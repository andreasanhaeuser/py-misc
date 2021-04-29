#!/usr/bin/python3

# standard
import unittest
import string

# to be tested
import iso

class PaperFormat(unittest.TestCase):
    def test_paper_format(self):
        print('Format, size (m, m)')
        for series in string.ascii_lowercase[:8]:
            formatter =                     '%s4, %1.3f, %1.3f'
            args = series.upper(), *iso.paper_size(series, 4)
            print(formatter % args)
        print('A4.5', iso.paper_size('A', 4.5))
        print('A(-4.5)', iso.paper_size('A', -4.5))


if __name__ == '__main__':
    unittest.main()
