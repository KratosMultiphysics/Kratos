# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
#                   Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------
# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# ==============================================================================
def CreateSolver(model_part, opt_settings):

    # Dictionary to store solvers of all response functions defined in the optimization settings
    solver = {}

    # Collect all responses
    specified_responses = {}
    for response_id in opt_settings.objectives:
        specified_responses[response_id] = opt_settings.objectives[response_id]

    if not specified_responses:
        raise ValueError("No objective function specified!")

    for response_id in opt_settings.constraints:
        specified_responses[response_id] = opt_settings.constraints[response_id]

    # Creat response function solver according to specified settings and add relevant variables
    in_active_responses = True
    if "strain_energy" in specified_responses.keys():
        in_active_responses = False
        model_part.AddNodalSolutionStepVariable(STRAIN_ENERGY_SHAPE_GRADIENT)
        solver["strain_energy"] = StrainEnergyResponseFunction(model_part, specified_responses["strain_energy"])
    if "mass" in specified_responses.keys():
        in_active_responses = False
        model_part.AddNodalSolutionStepVariable(MASS_SHAPE_GRADIENT)
        solver["mass"] = MassResponseFunction(model_part, specified_responses["mass"])        
    if in_active_responses:
        raise ValueError("Specified response function not implemented. Implemented response functions are: \"strain_energy\", \"mass\"")

    return solver

# ==============================================================================
