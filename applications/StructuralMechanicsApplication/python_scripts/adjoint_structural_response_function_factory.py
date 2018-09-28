from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
import adjoint_structural_response

def CreateAdjointResponseFunction(response_id, response_settings, model):
    response_type = response_settings["response_type"].GetString()

    if response_type == "adjoint_nodal_displacement":
        return adjoint_structural_response.AdjointResponseFunction(response_id, response_settings, model)

    elif response_type == "adjoint_linear_strain_energy":
        return adjoint_structural_response.AdjointResponseFunction(response_id, response_settings, model)

    elif response_type == "adjoint_local_stress":
        return adjoint_structural_response.AdjointResponseFunction(response_id, response_settings, model)

    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: 'adjoint_nodal_displacement', 'adjoint_linear_strain_energy', 'adjoint_local_stress'." )