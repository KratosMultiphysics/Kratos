from KratosMPI import *
import sys
import atexit

Hello()

MPIInitialize(sys.argv)

def finalize():
    MPIFinalize()
atexit.register(finalize)