from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics as KM
from KratosMetisApplication import *
application = KratosMetisApplication()
application_name = "KratosMetisApplication"

KM._ImportApplication(application, application_name)
