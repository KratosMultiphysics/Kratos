from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os.path
import sys
import KratosMultiphysics
from KratosMultiphysics import Logger


def ImportApplication(application, application_name, application_folder, caller, mod_path=None):
    Globals = KratosMultiphysics.KratosGlobals
    Kernel = Globals.Kernel
    applications_root = Globals.ApplicationsRoot

    Logger.PrintInfo("", "Importing  " + application_name)

    # Add python scrips folder to path
    application_path = os.path.join(applications_root, application_folder)
    python_path = os.path.join(application_path, 'python_scripts')
    sys.path.append(python_path)
    # Add constitutive laws python scrips folder to path
    constitutive_laws_path = os.path.join(python_path, 'constitutive_laws')
    sys.path.append(constitutive_laws_path)

    # adding the scripts in "APP_NAME/python_scripts" such that they are treated as a regular python-module
    if mod_path is not None: # optional for backwards compatibility
        mod_path.append(python_path)

    # Add application to kernel
    Kernel.ImportApplication(application)
