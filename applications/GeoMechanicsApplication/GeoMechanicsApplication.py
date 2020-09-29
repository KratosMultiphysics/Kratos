# # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# from __future__ import print_function, absolute_import, division

# # Application dependent names and paths
# import KratosMultiphysics as KM
# from KratosGeoMechanicsApplication import *
# application = KratosGeoMechanicsApplication()
# application_name = "KratosGeoMechanicsApplication"
# application_folder = "GeoMechanicsApplication"

# KM._ImportApplication(application, application_name, application_folder, __path__)


from KratosMultiphysics import _ImportApplication
from KratosGeoMechanicsApplication import *
application = KratosGeoMechanicsApplication()
application_name = "KratosGeoMechanicsApplication"

_ImportApplication(application, application_name)
