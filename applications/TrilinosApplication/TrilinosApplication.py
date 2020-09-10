from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import _ImportApplication
import KratosMultiphysics.mpi # importing the MPI-Core, since the TrilinosApp directly links to it
from KratosTrilinosApplication import *
application = KratosTrilinosApplication()
application_name = "KratosTrilinosApplication"

_ImportApplication(application, application_name)
