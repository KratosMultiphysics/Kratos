from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7


def import_solver(SolverSettings):
    """this function imports a solver named "solver_type" from SolverSettings
    solver_type is expected to be the FILENAME of the solver to be imported"""
    obj = __import__(SolverSettings.solver_type)
    return obj

def DeleteFileIfExisting(file_name):
    """This function tries to delete a file
    It uses try/except to also work in MPI
    """
    from os import remove
    try:
        remove(file_name)
    except:
        pass

def DeleteDirectoryIfExisting(directory_name):
    """This function tries to delete a folder
    It uses try/except to also work in MPI
    """
    from shutil import rmtree
    try:
        rmtree(directory_name)
    except:
        pass

def DeleteTimeFiles(directory_name):
    """This function deletes all *.time files in a directory
    """
    import os
    for file_name in os.listdir(directory_name):
        if file_name.endswith(".time"):
            DeleteFileIfExisting(os.path.join(directory_name, file_name))

def GetKratosMultiphysicsPath():
    """Returning the path to the KratosMultiphysics-module
    """
    import KratosMultiphysics, os
    return os.path.dirname(KratosMultiphysics.__file__)

def GetListOfAvailableApplications():
    """Return a list of compiled applications available in the
    KratosMultiphysics module.
    """
    kratos_path = GetKratosMultiphysicsPath()
    import os, re

    apps = [
        f.split('.')[0] for f in os.listdir(kratos_path) if re.match(r'.*Application*', f)
    ]

    return apps

def IsApplicationAvailable(application_name):
    """Returns whether an application is available
    """
    from KratosMultiphysics import Logger
    warn_msg  = '"IsApplicationAvailable" is deprecated, please use "CheckIfApplicationsAvailable" instead'
    Logger.PrintWarning("DEPRECATION", warn_msg)
    return CheckIfApplicationsAvailable(application_name)

def AreApplicationsAvailable(list_application_names):
    """Returns whether several applications are available
    """
    from KratosMultiphysics import Logger
    warn_msg  = '"AreApplicationsAvailable" is deprecated, please use "CheckIfApplicationsAvailable" instead'
    Logger.PrintWarning("DEPRECATION", warn_msg)
    return CheckIfApplicationsAvailable(*list_application_names)

def CheckIfApplicationsAvailable(*application_names):
    """Returns whether the inquired applications are available
    """
    available_apps = GetListOfAvailableApplications()

    for app_name in application_names:
        if app_name not in available_apps:
            return False

    return True
