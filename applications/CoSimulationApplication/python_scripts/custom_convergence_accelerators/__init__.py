# Adding the current directory to the path such that the modules can be imported with __import__ in the factory
from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7
import sys, os
sys.path.append(os.path.dirname(__file__))