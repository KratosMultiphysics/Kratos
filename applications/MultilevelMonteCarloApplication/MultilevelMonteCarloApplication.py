# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
from KratosMultiphysics import *
from KratosMultilevelMonteCarloApplication import *
application = KratosMultilevelMonteCarloApplication()
application_name = "KratosMultilevelMonteCarloApplication"
application_folder = "MultilevelMonteCarloApplication"

KratosMultiphysics._ImportApplicationAsModule(application, application_name, application_folder, __path__)