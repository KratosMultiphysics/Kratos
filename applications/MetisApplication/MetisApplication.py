print("1")
from KratosMultiphysics import _ImportApplication
print("2")
import KratosMultiphysics.mpi # importing the MPI-Core, since the MetisApp directly links to it
print("3")
from KratosMetisApplication import *
print("4")
application = KratosMetisApplication()
print("5")
application_name = "KratosMetisApplication"
print("6")

_ImportApplication(application, application_name)
print("7")
