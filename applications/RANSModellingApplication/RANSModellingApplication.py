# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import sys

# Application dependent names and paths
import KratosMultiphysics as KM
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
if (CheckIfApplicationsAvailable("TrilinosApplication")):
    from KratosMultiphysics.TrilinosApplication import *
from KratosRANSModellingApplication import *
application = KratosRANSModellingApplication()
application_name = "KratosRANSModellingApplication"
application_folder = "RANSModellingApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, sys.path)
