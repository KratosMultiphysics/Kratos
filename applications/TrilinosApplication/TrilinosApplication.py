from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import Logger, _ImportApplication
import KratosMultiphysics.mpi # importing the MPI-Core, since the TrilinosApp directly links to it
try:
    from KratosTrilinosApplication import *
except BaseException as error:
    Logger.PrintInfo("TrilinosApplication", "The import of the TrilinosApplication failed with the following error:")
    Logger.PrintInfo("TrilinosApplication", error)
    Logger.PrintInfo("TrilinosApplication", "Please visit the following website for known issues:")
    Logger.PrintInfo("TrilinosApplication", "https://github.com/KratosMultiphysics/Kratos/tree/master/applications/TrilinosApplication#known-issues")
application = KratosTrilinosApplication()
application_name = "KratosTrilinosApplication"

_ImportApplication(application, application_name)
