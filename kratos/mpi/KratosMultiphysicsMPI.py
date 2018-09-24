from KratosMPI import *
import sys
import atexit

Hello()
MPIInitialize(sys.argv)

atexit.register(MPIFinalize)