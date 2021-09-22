
from KratosMultiphysics.CompressiblePotentialFlowApplication import potential_flow_response

def CreateResponseFunction(response_id, response_settings, model):
    response_type = response_settings["response_type"].GetString()

    if response_type == "adjoint_lift_potential_jump":
        return potential_flow_response.AdjointResponseFunction(response_id, response_settings, model)
    elif response_type == "embedded_adjoint_lift_potential_jump":
        return potential_flow_response.EmbeddedAdjointResponseFunction(response_id, response_settings, model)
    elif response_type == "angle_of_attack":
        return potential_flow_response.AngleOfAttackResponseFunction(response_id, response_settings, model)
    elif response_type == "chord_length":
        return potential_flow_response.ChordLengthResponseFunction(response_id, response_settings, model)
    elif response_type == "perimeter":
        return potential_flow_response.PerimeterResponseFunction(response_id, response_settings, model)
    elif response_type == "adjoint_lift_potential_jump_with_perimeter":
        return potential_flow_response.LiftPerimeterResponseFunction(response_id, response_settings, model)
    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: 'adjoint_lift_potential_jump', 'angle_of_attack', 'chord_length'," +
                        ",'perimeter' and 'adjoint_lift_potential_jump_with_perimeter'." )
