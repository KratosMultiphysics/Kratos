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

def CheckForPreviousImport():

    first_caller = KratosGlobals.AuthorizedCaller[
        1]  # path to file where KratosMultiphysics (or submodules) was first imported
    my_caller = inspect.stack()[1][
        1]  # path to file where this function was called

    if my_caller == first_caller:
        msg = "\n***\n* KratosMultiphysics was imported from the first time from\n"
        msg += "* " + str(my_caller) + "\n"
        msg += "* This file defines auxiliary Python tools and assumes that the \n"
        msg += "* KratosMultiphysics module and all required applications will be\n"
        msg += "* first imported from a \"main\" Python script. Please import\n"
        msg += "* them before importing this file.\n*\n"
        msg += "* A commmon cause of this error is the use of the DEPRECATED\n"
        msg += "* applications_interface to import Kratos. To solve it please \n"
        msg += "* remove the old import system and COMPLETELY replace it with\n"
        msg += "* KratosMultiphysics modules in your main script\n"
        msg += "***\n"
        if KratosGlobals.ApplicationsInterfaceIsDeprecated:
            raise RuntimeError(msg)
        else:
            print("\n\n WARNING: It appears that you are still using applications_interface")
            print(" to import Kratos. This will soon be OBSOLETED.")
            print(" See the following message for details:\n")
            print(msg)
