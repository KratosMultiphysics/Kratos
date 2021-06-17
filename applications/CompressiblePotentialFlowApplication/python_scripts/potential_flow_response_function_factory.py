
from KratosMultiphysics.CompressiblePotentialFlowApplication import potential_flow_response

def CreateResponseFunction(response_id, response_settings, model):
    response_type = response_settings["response_type"].GetString()

    if response_type == "adjoint_lift_potential_jump":
        return potential_flow_response.AdjointResponseFunction(response_id, response_settings, model)
    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: 'adjoint_lift_potential_jump'." )
