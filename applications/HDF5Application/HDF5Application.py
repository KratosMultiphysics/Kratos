from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import _ImportApplication
from KratosHDF5Application import *
application = KratosHDF5Application()
application_name = "KratosHDF5Application"

_ImportApplication(application, application_name)
