from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosDEMApplication import *

from KratosMultiphysics import _ImportApplication
application = KratosDEMApplication()
application_name = "KratosDEMApplication"

_ImportApplication(application, application_name)
