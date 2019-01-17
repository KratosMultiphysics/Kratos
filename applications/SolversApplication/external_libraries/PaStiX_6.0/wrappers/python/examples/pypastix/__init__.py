"""
PyPaStiX
========

 @file __init__.py

 PaStiX python module intialization

 @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.1
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Louis Poirel
 @date 2018-07-16

"""
import ctypes
import ctypes.util
import spm

# Load the PASTIX library
libpastix_name = ctypes.util.find_library('pastix')
if libpastix_name == None:
    raise EnvironmentError("Could not find shared library: pastix."
                           "The path to libpastix.so should be in "
                           "$LIBRARY_PATH")
try:
    libpastix = ctypes.cdll.LoadLibrary(libpastix_name)
except:
    raise EnvironmentError("Could not load shared library: pastix."
                           "The path to libpastix.so should be in "
                           "$LD_LIBRARY_PATH or $DYLD_LIBRARY_PATH on MacOS");

__all__ = [ 'libpastix' ]

from .enum   import *
from .pastix import *
from .solver import *
