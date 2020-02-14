from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
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
    mpi_detected = (                         # Probing the environment to see if this is an MPI run
        "OMPI_COMM_WORLD_SIZE" in os.environ # OpenMPI implementation detected
        or "PMI_SIZE" in os.environ          # Intel MPI detected
        or "MPI_LOCALNRANKS" in os.environ   # Recent mpich detected
    )
    mpi_requested = "--using-mpi" in sys.argv[1:] # Forcing MPI initialization through command-line flag

    using_mpi = False
    if mpi_detected or mpi_requested:
        from KratosMultiphysics.kratos_utilities import IsMPIAvailable
        if IsMPIAvailable():
            import KratosMultiphysics.mpi
            mpi.InitializeMPIParallelRun()
            using_mpi = True
        else:
            if mpi_requested:
                msg = [
                    "\nThe MPI module could not be initialized."
                    "\nRequesting MPI support through the \"--using-mpi\" command line option, ",
                    "\nbut the Kratos mpi module could not be initialized. The most likely cause ",
                    "\nfor this warning is that this Kratos installation was compiled without ",
                    "\nMPI support."
                ]
                raise Exception("".join(msg))
            else: # if mpi_detected
                msg = [
                    "\nThis appears to be an MPI run, but the Kratos mpi module could not be found."
                    "\nMPI was not initialized. Running in serial mode."
                ]
                Logger.PrintWarning("KRATOS INITIALIZATION WARNING:", "".join(msg))

    return kratos_globals.KratosGlobalsImpl(Kernel(using_mpi), KratosPaths.kratos_applications)

KratosGlobals = __ModuleInitDetail()

def _ImportApplicationAsModule(application, application_name, application_folder, mod_path):
    Kernel = KratosGlobals.Kernel
    Logger.PrintInfo("", "Importing    " + application_name)

    # Add application to kernel
    Kernel.ImportApplication(application)

def IsDistributedRun():
    return KratosGlobals.Kernel.IsDistributedRun()
