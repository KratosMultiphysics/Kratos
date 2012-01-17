import os.path
import sys
import inspect
import kratos_globals

# this adds the libs/ and applications/ folders to sys.path
import KratosLoader

# import core library (Kratos.so)
from Kratos import *

KratosGlobals = kratos_globals.KratosGlobals(Kernel(),inspect.stack()[1],KratosLoader.kratos_applications)
