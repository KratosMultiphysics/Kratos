
from KratosMultiphysics.CompressiblePotentialFlowApplication import potential_flow_response
from KratosMultiphysics.CompressiblePotentialFlowApplication import airfoil_2d_responses

def CreateResponseFunction(response_id, response_settings, model):
    response_type = response_settings["response_type"].GetString()

    if response_type == "adjoint_lift_potential_jump":
        return potential_flow_response.AdjointResponseFunction(response_id, response_settings, model)
    elif response_type == "angle_of_attack":
        return airfoil_2d_responses.AngleOfAttackResponseFunction(response_id, response_settings, model)
    elif response_type == "chord_length":
        return airfoil_2d_responses.ChordLengthResponseFunction(response_id, response_settings, model)
    elif response_type == "perimeter":
        return airfoil_2d_responses.PerimeterResponseFunction(response_id, response_settings, model)
    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: 'adjoint_lift_potential_jump', 'angle_of_attack', 'chord_length'," +
                        " and 'perimeter'." )
