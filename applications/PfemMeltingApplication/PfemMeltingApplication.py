# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# import the applications the PfemMelting depends on
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.ConvectionDiffusionApplication
import KratosMultiphysics.MeshingApplication

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosPfemMeltingApplication import *
application = KratosPfemMeltingApplication()
application_name = "KratosPfemMeltingApplication"

_ImportApplication(application, application_name)
