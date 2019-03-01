from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os.path
import sys
from . import kratos_globals

class KratosPaths(object):
    kratos_install_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

    kratos_libs = os.path.join(kratos_install_path, "libs")
    kratos_applications = os.path.join(kratos_install_path, "applications")
    kratos_scripts = os.path.join(kratos_install_path, "kratos", "python_scripts")
    kratos_tests = os.path.join(kratos_install_path, "kratos", "tests")

# import core library (Kratos.so)
sys.path.append(KratosPaths.kratos_libs)
from Kratos import *

def __ModuleInitDetail():
    """
    Configure the parallel environment.
    This is done before initializing the kernel to ensure that MPI
    and the parallel DataCommunicator are initialized when the Kernel is built.
    It is defined as a function to avoid polluting the Kratos namespace with local variables.
    """
    if "--using-mpi" in sys.argv[1:]:
        try:
            import KratosMultiphysics.mpi
        except ModuleNotFoundError:
            msg = [
                "\nRequesting MPI support through the \"--using-mpi\" command line option, ",
                "\nbut the mpi environment could not be initialized. The most likely cause ",
                "\nfor this warning is that this Kratos installation was compiled without ",
                "\nMPI support."
            ]
            Logger.PrintWarning("KRATOS INITIALIZATION WARNING:", "".join(msg))

__ModuleInitDetail()

KratosGlobals = kratos_globals.KratosGlobalsImpl(
    Kernel(), KratosPaths.kratos_applications)

# adding the scripts in "kratos/python_scripts" such that they are treated as a regular python-module
__path__.append(KratosPaths.kratos_scripts)
# To be purely pythonic, the following line should be removed
# and all imports of files in python_scrips should be made relative to the KratosMultiphysics module.
sys.path.append(KratosPaths.kratos_scripts)

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

