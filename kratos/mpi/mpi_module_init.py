from __future__ import absolute_import
import sys

from KratosMPI import * # this should probably use the RTLD flags too (since this is related to loading OpenMPI itself)

if sys.platform.startswith('linux'):
    # Note: from Python 3.3 onwards, dll load flags are available from module os
    # from Python 3.6 onwards, module DLFCN no longer exists
    flags = sys.getdlopenflags()
    if sys.version_info >= (3,3):
        import os
        dll_load_flags = os.RTLD_NOW | os.RTLD_GLOBAL
    else:
        import DLFCN as dl
        dll_load_flags = dl.RTLD_NOW | dl.RTLD_GLOBAL
    sys.setdlopenflags(dll_load_flags)
    import mpipython
    sys.setdlopenflags(flags)
else:
    import mpipython

mpi = mpipython.GetMPIInterface()

Hello()
