# KRATOS  _     _                       ____        _
#        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
#        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
#        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
#        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
#
#   Author: Thomas Oberbichler
#

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosLinearSolversApplication import *
application = KratosLinearSolversApplication()
application_name = "LinearSolversApplication"

_ImportApplication(application, application_name)
