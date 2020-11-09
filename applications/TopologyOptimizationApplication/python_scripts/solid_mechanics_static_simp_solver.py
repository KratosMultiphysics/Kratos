# ==============================================================================
#  TopologyOptimizationApplication
#
#  License:         BSD License
#                   license: TopologyOptimizationApplication/license.txt
#
#  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.TopologyOptimizationApplication as KratosTopology

# Import the mechanical static solver base class
import solid_mechanics_static_solver

# ==============================================================================
def CreateSolver(main_model_part, custom_settings):
    return StaticMechanicalSIMPSolver(main_model_part, custom_settings)

# ==============================================================================
class StaticMechanicalSIMPSolver(solid_mechanics_static_solver.StaticMechanicalSolver):

    #### Specific internal functions ####

    def _GetSolutionScheme(self, analysis_type, component_wise, compute_contact_forces):

        if(analysis_type == "Linear"):
            mechanical_scheme = KratosTopology.ResidualBasedIncrementalUpdateStaticSIMPScheme()

        elif(analysis_type == "Non-Linear" ):
            raise TypeError("Non-Linear analysis_type not yet supported!")

        return mechanical_scheme
# ==============================================================================
