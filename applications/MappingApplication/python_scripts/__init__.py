# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
import KratosMultiphysics as KM
from KratosMappingApplication import *
application = KratosMappingApplication()
application_name = "KratosMappingApplication"
application_folder = "MappingApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)
