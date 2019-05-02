#
#   KRATOS _______
#         / ____(_)___ ____  ____
#        / __/ / / __ `/ _ \/ __ \
#       / /___/ / /_/ /  __/ / / /
#      /_____/_/\__, /\___/_/ /_/ SolversApplication
#              /____/
#
#   Author: Thomas Oberbichler
#

import KratosMultiphysics as KM
from KratosEigenSolversApplication import *
application = KratosEigenSolversApplication()
application_name = "KratosEigenSolversApplication"
application_folder = "EigenSolversApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)
