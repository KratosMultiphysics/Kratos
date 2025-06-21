import os, sys, shutil, re
import KratosMultiphysics as KM

def IssueDeprecationWarning(label, *message):
    KM.Logger.PrintWarning('DEPRECATION-Warning; '+label, ' '.join(map(str,message)))

def DeleteFileIfExisting(file_name):
    """This function tries to delete a file
    It uses try/except to also work in MPI
    """
    try:
        os.remove(file_name)
    except:
        pass

def DeleteDirectoryIfExisting(directory_name):
    """This function tries to delete a folder
    It uses try/except to also work in MPI
    """
    try:
        shutil.rmtree(directory_name)
    except:
        pass

def DeleteFilesEndingWith(directory_name, *ends_with):
    """This function deletes all the files in a directory that
    with the provided argument
    """
    for file_name in os.listdir(directory_name):
        for arg in ends_with:
            if file_name.endswith(arg):
                DeleteFileIfExisting(os.path.join(directory_name, file_name))

def DeleteTimeFiles(directory_name):
    """This function deletes all *.time files in a directory
    """
    DeleteFilesEndingWith(directory_name, ".time")

def GetKratosMultiphysicsPath():
    """Returning the path to the KratosMultiphysics-module
    """
    return os.path.dirname(KM.__file__)

def GetListOfAvailableApplications():
    """Return a list of compiled applications available in the
    KratosMultiphysics module.
    """
    kratos_path = GetKratosMultiphysicsPath()

    apps = sorted([
        f.split('.')[0] for f in os.listdir(kratos_path) if re.match(r'.*Application*', f)
    ])

    return apps

def GetListOfImportedApplications(module_name="KratosMultiphysics"):
    """Returns the list of moldules that are or were imported at
    call time, including the main module_name module.
    """
    loaded_modules = [mod for mod in sys.modules.keys() if '__file__' in dir(sys.modules[mod]) and sys.modules[mod].__file__]
    kratos_modules = [mod for mod in loaded_modules if module_name in sys.modules[mod].__file__]

    kratos_apps = [mod for mod in kratos_modules if os.path.basename(sys.modules[mod].__file__) == "__init__.py"]

    return kratos_apps

def GetMapOfImportedApplications(module_name="KratosMultiphysics"):
    """Returns the tree of modules that are or were imported at
    call time, including the main "module_name" module.
    """
    import_map = {}
    imported_apps = GetListOfImportedApplications(module_name)

    for app in imported_apps:
        names = app.split(".")

        entry_point = import_map
        for part in names:
            if part not in entry_point:
                entry_point[part] = {}
            entry_point = entry_point[part]

    return import_map

def IsMPIAvailable():
    """Check if the KratosMPI module (the MPI core) is available.
    """
    kratos_path = GetKratosMultiphysicsPath()

    return "mpi" in os.listdir(kratos_path)

def CheckIfApplicationsAvailable(*application_names):
    """Returns whether the inquired applications are available
    """
    available_apps = GetListOfAvailableApplications()
    return all(app_name in available_apps for app_name in application_names)

def GetNotAvailableApplications(*application_names):
    """Returns a list of applications that are not available out of a provided list of applications
    """
    available_apps = GetListOfAvailableApplications()
    return list(filter(lambda app_name: app_name not in available_apps, application_names))

def GenerateVariableListFromInput(param):
    """ Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
    """
    return [ KM.KratosGlobals.GetVariable( var_name ) for var_name in param.GetStringArray() ]

def GenerateFlagsListFromInput(param):
    """ Retrieve flag name from input (a string) and request the corresponding C++ object to the kernel
    """
    return [ KM.KratosGlobals.GetFlag( flag_name ) for flag_name in param.GetStringArray() ]
