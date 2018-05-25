from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
import structural_response

def CreateResponseFunction(response_id, response_settings, model_part):
    response_type = response_settings["response_type"].GetString()

    if response_type == "strain_energy":
        return structural_response.StrainEnergyResponseFunction(response_id, response_settings, model_part)

    elif response_type == "mass":
        return structural_response.MassResponseFunction(response_id, response_settings, model_part)

    elif response_type == "eigenfrequency":
        return structural_response.EigenFrequencyResponseFunction(response_id, response_settings, model_part)

    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: 'mass', 'strain_energy', 'eigenfrequency'." )