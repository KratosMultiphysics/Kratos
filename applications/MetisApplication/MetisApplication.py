print("1", flush=True)
from KratosMultiphysics import _ImportApplication
print("2", flush=True)
import KratosMultiphysics.mpi # importing the MPI-Core, since the MetisApp directly links to it
print("3", flush=True)
from KratosMetisApplication import *
print("4", flush=True)
application = KratosMetisApplication()
print(application)
print("5", flush=True)
application_name = "KratosMetisApplication"
print("6", flush=True)

_ImportApplication(application, application_name)
print("7", flush=True)
