#!/usr/bin/python
"""Performance info for loops.

    See docstring of Chronometer for more info.

    Author
    ------
    Andreas Anhaeuser (AA) <andreas.anhaeuser@posteo.net>
    Institute for Geophysics and Meteorology, University of Cologne, Germany
    Greenpeace East Asia

    History
    -------
    2015       (AA): Created
    2017-01-25 (AA): Changed layout
    2018-01-25 (AA): Updated docstring.
"""
from __future__ import print_function

# standard modules
import os
import sys
from collections import Iterable
from copy import deepcopy as copy
import datetime as dt
import inspect
import textwrap
import warnings
if sys.version_info.major < 3:
    import __builtin__ as builtin
else:
    import builtins as builtin

# PyPI modules
import numpy as np

# local modules
if __name__ == '__main__':
    import string_utils
    from timer import Timer
    import chronometer_utils as utils
else:
    from . import string_utils
    from .timer import Timer
    from . import chronometer_utils as utils

_builtin_print = copy(builtin.print)
_builtin_warning = warnings.showwarning


class Node(object):
    def __init__(self, iterable):
        self.timer = Timer()
        self.total_count = len(iterable)
        self.count = 0
        self.children = []

class Chronometer(object):
    def __init__(
            self, header='', info='', file=None, print_colors=None,
            item_name='loop', item_plural=None, report_level='normal',
            ):
        self.root = []
        self.pivot = [0]
        self.iterables = []
        self.iter_pos = []

    def iterate(self, iterable):
        if not isinstance(iterable, Iterable):
            raise TypeError('Expected Iterable, got %s' % type(iterable))
        next_node = self.get_next_node()
        if next_node is None:
            self.append_node(iterable)
            self.iter_pos.append(0)
            self.iterable = iterable

        active_node.append(node.iterable)

    def __iter__(self):
        return self

    def __next__(self):
        i = self.iter_pos[-1]
        iterable = self.iterables[-1]
        N = len(iterable)
        node = self.get_active_node()
        if i <= N:
            self.pivot.pop()
            self.iterables.pop()
            node.timer.stop()
            
            raise StopIteration

        item = iterable[i]
        self.iter_pos[-1] += 1
        return item
        

    def range(self, N):
        return self.iterate(range(N))

    def enumerate(self, iterable):
        pass

    def get_next_node(self):
        i = pivot[-1]
        active_node = self.get_active_node()
        try:
            next_node = active_node.children[i]
        except IndexError:
            next_node = None
        return next_node
    def get_active_node(self):
        node = self.root
        for index in pivot[:-1]:
            node = node.children[index]
        return node
