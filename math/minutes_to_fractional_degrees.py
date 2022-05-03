#!/usr/bin/python3
"""Convert angles D MM SS to D.DDDDD """

import sys

args = sys.argv[1:]

assert 1 <= len(args) <= 3

deg = float(args[0])
m = 0
s = 0

if len(args) > 1:
    m = float(args[1])

if len(args) > 2:
    s = float(args[2])

out = deg + m / 60 + s / 3600
print(out)
