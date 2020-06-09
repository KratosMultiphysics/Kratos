import sys
import os

if sys.platform.startswith('linux'):
    # see https://github.com/open-mpi/ompi/issues/3705 for details
    flags = sys.getdlopenflags()
    dll_load_flags = os.RTLD_NOW | os.RTLD_GLOBAL
    sys.setdlopenflags(dll_load_flags)

from KratosMPI import *

if sys.platform.startswith('linux'):
    # restore default system flags
    sys.setdlopenflags(flags)
