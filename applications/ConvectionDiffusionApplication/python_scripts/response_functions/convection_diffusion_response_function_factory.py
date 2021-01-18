from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics.ConvectionDiffusionApplication.response_functions import convection_diffusion_response

def CreateResponseFunction(response_id, response_settings, model):
    response_type = response_settings["response_type"].GetString()

    if response_type == "point_temperature":
        return convection_diffusion_response.AdjointResponseFunction(response_id, response_settings, model)

    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: 'point_temperature'." )
