from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
import structural_response

def CreateResponseFunction(response_id, response_settings, model_part):
    response_type = response_settings["response_type"].GetString()

    if response_type == "strain_energy":
        response_function_utility = StructuralMechanicsApplication.StrainEnergyResponseFunctionUtility( model_part, response_settings )
        return structural_response.SimpleResponseFunctionWrapper(response_id, response_settings, response_function_utility, model_part)

    elif response_type == "mass":
        response_function_utility = StructuralMechanicsApplication.MassResponseFunctionUtility( model_part, response_settings )
        return structural_response.MassResponseFunctionWrapper(response_id, response_settings, response_function_utility, model_part)

    elif response_type == "eigenfrequency":
        if not response_settings.Has("weighting_method") or response_settings["weighting_method"].GetString() == "none":
            response_function_utility = StructuralMechanicsApplication.EigenfrequencyResponseFunctionUtility( model_part, response_settings )
        elif response_settings["weighting_method"].GetString() == "linear_scaling":
            response_function_utility = StructuralMechanicsApplication.EigenfrequencyResponseFunctionLinScalUtility( model_part, response_settings )
        else:
            raise NameError("The following weighting_method is not valid for eigenfrequency response: " + response_settings["weighting_method"].GetString() +
                            ".\nAvailable weighting methods are: 'none', 'linear_scaling'. Default: 'none'")
        return structural_response.SimpleResponseFunctionWrapper(response_id, response_settings, response_function_utility, model_part)

    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: 'mass', 'strain_energy', 'eigenfrequency'." )