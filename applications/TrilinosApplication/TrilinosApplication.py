from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import _ImportApplication
import KratosMultiphysics.mpi # importing the MPI-Core, since the TrilinosApp directly links to it
try:
    from KratosTrilinosApplication import *
except BaseException as error:
    print("The import of the TrilinosApplication failed with the following error:")
    print(error)
    print("Please visit the following website for known issues:")
    print("https://github.com/KratosMultiphysics/Kratos/tree/master/applications/TrilinosApplication#known-issues")
application = KratosTrilinosApplication()
application_name = "KratosTrilinosApplication"

_ImportApplication(application, application_name)
