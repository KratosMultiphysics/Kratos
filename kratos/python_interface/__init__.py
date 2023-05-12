import os
import sys

# This is a "dirty" fix to force python to keep loading shared libraries from
# the PATH in windows (See https://docs.python.org/3/library/os.html#os.add_dll_directory)
# THIS NEEDS TO BE EXECUTED BEFORE ANY DLL / DEPENDENCY IS LOADED.
if os.name == 'nt': # This means "Windows"
    for lib_path in os.environ["PATH"].split(os.pathsep):
        try:
            os.add_dll_directory(lib_path.replace("\\","/"))
        except Exception as e:
            # We need to capture possible invalid paths.
            pass

from . import kratos_globals
from . import python_registry
from . import python_registry_lists
from . import python_registry_utilities

if sys.version_info < (3, 8):
    raise Exception("Kratos only supports Python version 3.8 and above")

class KratosPaths(object):
    kratos_install_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

    kratos_libs = os.path.join(kratos_install_path, "libs")
    kratos_module_libs = os.path.join(os.path.abspath(os.path.dirname(__file__)), ".libs")
    kratos_applications = os.path.join(kratos_install_path, "applications")
    kratos_scripts = os.path.join(kratos_install_path, "kratos", "python_scripts")
    kratos_tests = os.path.join(kratos_install_path, "kratos", "tests")

# import core library (Kratos.so)
sys.path.append(KratosPaths.kratos_libs)
sys.path.append(KratosPaths.kratos_module_libs)
from Kratos import *

Kernel.RegisterPythonVersion()

def __getattr__(name):
    if name == "CppRegistry":
        err_msg = "c++ registry must be accessed through 'KratosMultiphysics.Registry'."
        raise Exception(err_msg)
    raise AttributeError(f"Module {__name__} has no attribute {name}.")

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

                from KratosMultiphysics import ParallelEnvironment
                exc_value = exc_type("(on rank {}): {}".format(ParallelEnvironment.GetDataCommunicator("World").Rank(), exc_value))
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

    return kratos_globals.KratosGlobalsImpl(Kernel(using_mpi), KratosPaths.kratos_applications)

KratosGlobals = __ModuleInitDetail()

# Create the Python global registry
# Note that this interfaces the c++ registry.
Registry = python_registry.PythonRegistry()
RegisterPrototype = python_registry.RegisterPrototype

# Remove CppRegistry from locals() in order to give preference to the exception thrown in __getattr__
# This is required since we cannot use properties as usual due to the fact that we have no instance of CppRegistry (it is a static variable in c++)
locals().pop("CppRegistry")

# Loop and register the Python registry lists
python_registry_utilities.RegisterAll("KratosMultiphysics", python_registry_lists)

# Detect kratos library version
python_version = KratosGlobals.Kernel.PythonVersion()
python_version = python_version.replace("Python","")
kratos_version_info = python_version.split(".")

if sys.version_info.major != int(kratos_version_info[0]) and sys.version_info.minor != int(kratos_version_info[1]):
    Logger.PrintWarning("Warning: Kratos is running with python {}.{} but was compiled with python {}.{}. Please ensure the versions match.".format(
        sys.version_info.major, sys.version_info.minor,
        kratos_version_info[0], kratos_version_info[1]
    ))

# print the process id e.g. for attaching a debugger
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


# iterating through the parameters is deprecated
# the following wraps the __iter__ method to issue a deprecation warning
list_deprecation_warnings = []
orig_iter = Parameters.__iter__
import inspect
def iter_wrapper(self):
    # get information where the function is called
    # this is necessary to issue the deprecation warning only
    # once per call location
    frame = inspect.stack()[1]
    filename = frame.filename
    line_number = frame.lineno
    tup = (filename, line_number)
    # issue deprecation warning only once, providing file name and line number
    if tup not in list_deprecation_warnings:
        list_deprecation_warnings.append(tup)
        print(f'Deprecated method called in "{frame.filename}" in line {frame.lineno}: Iterating through "Parameters" object is deprecated, please use the "values" method instead')
    return orig_iter(self)
Parameters.__iter__ = iter_wrapper
