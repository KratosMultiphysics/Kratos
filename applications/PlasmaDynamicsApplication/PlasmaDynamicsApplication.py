from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import _ImportApplication

# import the applications the PlasmaDynamicsApplication depends on
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.ConvectionDiffusionApplication
import KratosMultiphysics.SwimmingDEMApplication
import KratosMultiphysics.PlasmaDynamicsApplication
from KratosPlasmaDynamicsApplication import *

from KratosMultiphysics.PlasmaDynamicsApplication.plasma_dynamics_analysis import PlasmaDynamicsAnalysis
application = KratosPlasmaDynamicsApplication()
application_name = "KratosPlasmaDynamicsApplication"

_ImportApplication(application, application_name)
