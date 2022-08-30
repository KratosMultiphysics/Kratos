from KratosMultiphysics import _ImportApplication
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
import KratosMultiphysics.mpi # importing the MPI-Core, since the TrilinosApp directly links to it
if CheckIfApplicationsAvailable("MetisApplication"):
    import KratosMultiphysics.MetisApplication # importing to avoid problems with loading the metis library
from KratosTrilinosApplication import *
application = KratosTrilinosApplication()
application_name = "KratosTrilinosApplication"

_ImportApplication(application, application_name)
