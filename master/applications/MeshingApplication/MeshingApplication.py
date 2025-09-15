# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosMultiphysics import IsDistributedRun
if IsDistributedRun():
    import KratosMultiphysics.mpi # importing the MPI-Core
from KratosMeshingApplication import *
application = KratosMeshingApplication()
application_name = "KratosMeshingApplication"

_ImportApplication(application, application_name)

from KratosMultiphysics import python_registry_utilities
from . import python_registry_lists
python_registry_utilities.RegisterAll("KratosMultiphysics.MeshingApplication", python_registry_lists)