from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os.path
from . import kratos_globals

# this adds the libs/ and applications/ folders to sys.path
from . import KratosLoader

# import core library (Kratos.so)
from Kratos import *

KratosGlobals = kratos_globals.KratosGlobalsImpl(
    Kernel(), KratosLoader.kratos_applications)

# adding the scripts in "kratos/python_scripts" such that they are treated as a regular python-module
__path__.append(KratosLoader.kratos_scripts)

def _ImportApplicationAsModule(application, application_name, application_folder, mod_path):
    Kernel = KratosGlobals.Kernel
    applications_root = KratosGlobals.ApplicationsRoot

    Logger.PrintInfo("", "Importing    " + application_name)

    # adding the scripts in "APP_NAME/python_scripts" such that they are treated as a regular python-module
    application_path = os.path.join(applications_root, application_folder)
    python_path = os.path.join(application_path, 'python_scripts')
    mod_path.append(python_path)

    # Add application to kernel
    Kernel.ImportApplication(application)


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

def CheckRegisteredApplications(*applications):
    for app in applications:
       if not KratosGlobals.Kernel.IsImported(app):
           import __main__
           raise Exception("Application "+ app + " was not imported in the main script ("+__main__.__file__+")")
