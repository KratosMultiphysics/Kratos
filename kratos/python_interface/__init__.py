import os
import re
import sys
from . import kratos_globals

if sys.version_info < (3, 5):
    raise Exception("Kratos only supports Python version 3.5 and above")

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

    # Detect if MPI is used
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

            def CustomExceptionHook(exc_type, exc_value, exc_traceback):
                """Custom exception hook that also prints the source rank where the exception was thrown."""
                if issubclass(exc_type, KeyboardInterrupt):
                    # call the default excepthook saved at __excepthook__
                    sys.__excepthook__(exc_type, exc_value, exc_traceback)
                    return

                from KratosMultiphysics import DataCommunicator
                exc_value = exc_type("(on rank {}): {}".format(DataCommunicator.GetDefault().Rank(), exc_value))
                sys.__excepthook__(exc_type, exc_value, exc_traceback)

            sys.excepthook = CustomExceptionHook

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

    # Try to detect kratos library version
    try:
        kre = re.compile('Kratos\.([^\d]+)(\d+).+')
        kratos_version_info = [(kre.match(f))[2] for f in os.listdir(KratosPaths.kratos_libs) if kre.match(f)][0]

        if sys.version_info.major != int(kratos_version_info[0]) and sys.version_info.minor != int(kratos_version_info[1]):
            print("Warning: Kratos is running with python {}.{} but was compiled with python {}.{}. Please ensure the versions match.".format(
                sys.version_info.major, sys.version_info.minor, 
                kratos_version_info[0], kratos_version_info[1]
            ))
    except:
        print("Warning: Could not determine python version used to build kratos.")

    return kratos_globals.KratosGlobalsImpl(Kernel(using_mpi), KratosPaths.kratos_applications)

KratosGlobals = __ModuleInitDetail()

# print the process id e.g. for attatching a debugger
if KratosGlobals.Kernel.BuildType() != "Release":
    Logger.PrintInfo("Process Id", os.getpid())

def _ImportApplicationAsModule(application, application_name, application_folder, mod_path):
    Kernel = KratosGlobals.Kernel
    Logger.PrintInfo("", "Importing    " + application_name)
    Logger.PrintWarning('DEPRECATION-Warning', 'For importing "{}": "_ImportApplicationAsModule" is deprecated, please use "_ImportApplication"'.format(application_name))

    # Add application to kernel
    Kernel.ImportApplication(application)

def _ImportApplication(application, application_name):
    Kernel = KratosGlobals.Kernel
    Logger.PrintInfo("", "Importing    " + application_name)

    # Add application to kernel
    Kernel.ImportApplication(application)

def IsDistributedRun():
    return KratosGlobals.Kernel.IsDistributedRun()
