from KratosMultiphysics import _ImportApplication
import KratosMultiphysics.mpi # importing the MPI-Core, since the MetisApp directly links to it
from KratosMetisApplication import *
application = KratosMetisApplication()
application_name = "KratosMetisApplication"

_ImportApplication(application, application_name)
