# Adding the current directory to the path such that the modules can be imported with __import__ in the factory
from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7
import sys, os
current_dir_name = os.path.dirname(__file__)
sys.path.append(current_dir_name)

# Append all the sub directories containing the solver interfaces to the path
for dirName, subdirList, fileList in os.walk(current_dir_name):
    sys.path.append(dirName)