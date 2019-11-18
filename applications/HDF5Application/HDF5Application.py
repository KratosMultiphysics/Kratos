from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics as KM
from KratosHDF5Application import *
application = KratosHDF5Application()
application_name = "KratosHDF5Application"
application_folder = "HDF5Application"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)
