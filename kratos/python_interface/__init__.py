from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os.path
import sys
import inspect
from . import kratos_globals

from .kratos_utilities import *

# this adds the libs/ and applications/ folders to sys.path
from . import KratosLoader

# import core library (Kratos.so)
from Kratos import *

KratosGlobals = kratos_globals.KratosGlobals(
    Kernel(), inspect.stack()[1], KratosLoader.kratos_applications)

# Initialize Kernel so that core variables have an assigned Key even if we are not importing applications
KratosGlobals.Kernel.Initialize()

# adding the scripts in "kratos/python_scripts" such that they are treated as a regular python-module
__path__.append(KratosLoader.kratos_scripts)

def CheckForPreviousImport():
    warn_msg  = '"CheckForPreviousImport" is not needed any more and can be safely removed\n'
    warn_msg += 'It does nothing any more'
    Logger.PrintWarning('DEPRECATION', warn_msg)

def CheckRegisteredApplications(*applications):
    warn_msg  = '"CheckRegisteredApplications" is not needed any more and can be safely removed\n'
    warn_msg += 'It does nothing any more'
    Logger.PrintWarning('DEPRECATION', warn_msg)
