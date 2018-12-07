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
    warn_msg  = '"CheckForPreviousImport" is not needed any more and can be safely removed\n'
    warn_msg += 'It does nothing any more'
    Logger.PrintWarning('DEPRECATION', warn_msg)

def CheckRegisteredApplications(*applications):
    warn_msg  = '"CheckRegisteredApplications" is not needed any more and can be safely removed\n'
    warn_msg += 'It does nothing any more'
    Logger.PrintWarning('DEPRECATION', warn_msg)
