from KratosMultiphysics import _ImportApplication
from KratosHDF5Application import *
application = KratosHDF5Application()
application_name = "KratosHDF5Application"

_ImportApplication(application, application_name)
