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

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosEigenSolversApplication import *
application = KratosEigenSolversApplication()
application_name = "LinearSolversApplication"

_ImportApplication(application, application_name)
