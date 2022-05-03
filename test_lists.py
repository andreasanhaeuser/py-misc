#!/usr/bin/python3

from collections import Iterable
# PyPI
import numpy as np

# misc
from misc import lists

rand = np.random.random(12)
iterable = rand > 0.5

subints = lists.contiguous_intervals(iterable)
print(iterable* 1.)
for si in subints:
    print(si*1.)
