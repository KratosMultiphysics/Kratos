from KratosMPI import *
import sys
import atexit

Hello()
MPIInitialize()

atexit.register(MPIFinalize)