from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os.path
import sys
import KratosMultiphysics
from KratosMultiphysics import Logger


def ImportApplication(application, application_name, application_folder, caller, mod_path=None):
    KratosGlobals = KratosMultiphysics.KratosGlobals
    Kernel = KratosGlobals.Kernel
    applications_root = KratosGlobals.ApplicationsRoot

    Logger.PrintInfo("", "Importing  " + application_name)

    # Add python scrips folder to path
    application_path = os.path.join(applications_root, application_folder)
    python_path = os.path.join(application_path, 'python_scripts')
    sys.path.append(python_path)
    # Add constitutive laws python scrips folder to path
    constitutive_laws_path = os.path.join(python_path, 'constitutive_laws')
    sys.path.append(constitutive_laws_path)

    warn_msg  = '\nThe python-import-mechanism used for application "' + application_name
    warn_msg += '" is DEPRECATED!\n'
    warn_msg += 'Please check the following website for instuctions on how to update it:\n'
    warn_msg += 'https://github.com/KratosMultiphysics/Kratos/wiki/Applications-as-python-modules\n'
    warn_msg += 'The old mechanism will be removed on 01.10.2019!\n'
    Logger.PrintWarning('\n\x1b[1;31mDEPRECATION-WARNING\x1b[0m', warn_msg)

    # adding the scripts in "APP_NAME/python_scripts" such that they are treated as a regular python-module
    if mod_path is not None: # optional for backwards compatibility
        mod_path.append(python_path)

    # Add application to kernel
    Kernel.ImportApplication(application)
