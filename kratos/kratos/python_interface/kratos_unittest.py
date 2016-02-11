from __future__ import print_function, absolute_import, division
from unittest import *


def assertEqualTolerance(first, second, tolerance):
    if first < second + tolerance and first > second - tolerance:
        return True
    return False
