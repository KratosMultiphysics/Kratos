from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import * #TODO move variable to struct_mech
import structural_response

def CreateResponseFunction(response_id, response_settings, model_part):
    response_type = response_settings["response_type"].GetString()

    if response_type == "strain_energy":
        model_part.AddNodalSolutionStepVariable(STRAIN_ENERGY_SHAPE_GRADIENT)
        response_function_utility = StrainEnergyResponseFunction( model_part, response_settings )
        return structural_response.SimpleResponseFunctionWrapper(response_id, response_settings, response_function_utility, model_part)

    elif response_type == "mass":
        model_part.AddNodalSolutionStepVariable(MASS_SHAPE_GRADIENT)
        response_function_utility = MassResponseFunction( model_part, response_settings )
        return structural_response.MassResponseFunctionWrapper(response_id, response_settings, response_function_utility)

    elif response_type == "eigenfrequency":
        model_part.AddNodalSolutionStepVariable(EIGENFREQUENCY_SHAPE_GRADIENT)
        if not response_settings.Has("weighting_method") or response_settings["weighting_method"].GetString() == "none":
            response_function_utility = EigenfrequencyResponseFunction( model_part, response_settings )
        elif response_settings["weighting_method"].GetString() == "linear_scaling":
            response_function_utility = EigenfrequencyResponseFunctionLinScal( model_part, response_settings )
        else:
            raise NameError("The following weighting_method is not valid for eigenfrequency response: " + response_settings["weighting_method"].GetString())
        return structural_response.SimpleResponseFunction(response_id, response_settings, response_function_utility, model_part)

    elif response_type == "adjoint_nodal_displacement":
        return structural_response.AdjointResponseFunction(response_id, response_settings, model_part)

    else:
        raise NameError("The type of the following response function is not specified: " + response_id)