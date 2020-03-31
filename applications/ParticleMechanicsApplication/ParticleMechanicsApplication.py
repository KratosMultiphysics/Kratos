from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics as KM
from KratosParticleMechanicsApplication import *
application = KratosParticleMechanicsApplication()
application_name = "KratosParticleMechanicsApplication"

KM._ImportApplication(application, application_name)
