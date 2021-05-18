from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import _ImportApplication
from KratosFSIApplication import *
application = KratosFSIApplication()
application_name = "KratosFSIApplication"

_ImportApplication(application, application_name)
