#
# @file __init__.py
#
# SParse Matrix package python module intialization
#
# @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
# @version 6.0.0
# @author Pierre Ramet
# @author Mathieu Faverge
# @author Louis Poirel
# @date 2017-05-04
#
"""
PySpm
=====

Provides
  1. A sparse matrix structure
  2. Mathematical operations over this sparse matrix structure
  3. Driver to read from different file formats and to convert from Scipy package

"""
import ctypes
import ctypes.util

# Load the SPM library
libspm_name = ctypes.util.find_library('spm')
if libspm_name == None:
    raise EnvironmentError("Could not find shared library: spm."
                           "The path to libpastix_spm.so should be in "
                           "$LIBRARY_PATH")

try:
    libspm = ctypes.cdll.LoadLibrary(libspm_name)
except:
    raise EnvironmentError("Could not load shared library: spm."
                           "The path to libpastix_spm.so should be in "
                           "$LD_LIBRARY_PATH or $DYLD_LIBRARY_PATH on MacOS");

from .enum   import *
from .spm    import *

__all__ = [ 'libspm' ]

#__all__.extend(enum.__all__)
#__all__.extend(spm.__all__)
