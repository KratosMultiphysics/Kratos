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

# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosEigenSolversApplication import *
application = KratosEigenSolversApplication()
application_name = "EigenSolversApplication"

_ImportApplication(application, application_name)
