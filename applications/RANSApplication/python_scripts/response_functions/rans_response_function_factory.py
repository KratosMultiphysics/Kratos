# importing the Kratos Library
from KratosMultiphysics.RANSApplication.response_functions.lift_to_drag_response_function import LiftToDragResponseFunction
from KratosMultiphysics.RANSApplication.response_functions.average_location_deviation_response_function import AverageLocationDeviationResponseFunction

def CreateResponseFunction(response_id, response_settings, model):
    response_type = response_settings["response_type"].GetString()

    if response_type == "lift_to_drag_response":
        return LiftToDragResponseFunction(response_id, response_settings, model)
    elif response_type == "average_location_deviation_response":
        return AverageLocationDeviationResponseFunction(response_id, response_settings, model)
    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: " +
                        "\n\t 'lift_to_drag_response'" +
                        "\n\t 'average_location_deviation_response'." )
