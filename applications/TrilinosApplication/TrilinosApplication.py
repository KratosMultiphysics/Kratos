from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics as KM
import KratosMultiphysics.mpi # importing the MPI-Core, since the TrilinosApp directly links to it
from KratosTrilinosApplication import *
application = KratosTrilinosApplication()
application_name = "KratosTrilinosApplication"
application_folder = "TrilinosApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)
