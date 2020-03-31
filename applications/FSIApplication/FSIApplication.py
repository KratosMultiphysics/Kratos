from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics as KM
from KratosFSIApplication import *
application = KratosFSIApplication()
application_name = "KratosFSIApplication"

KM._ImportApplication(application, application_name)
