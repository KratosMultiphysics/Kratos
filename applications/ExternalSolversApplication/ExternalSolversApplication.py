from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics as KM
from KratosExternalSolversApplication import *
application = KratosExternalSolversApplication()
application_name = "KratosExternalSolversApplication"
application_folder = "ExternalSolversApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)
