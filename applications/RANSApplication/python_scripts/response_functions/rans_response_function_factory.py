# importing the Kratos Library
from KratosMultiphysics.RANSApplication.response_functions.lift_to_drag import LiftToDrag
from KratosMultiphysics.RANSApplication.response_functions.drag import Drag
from KratosMultiphysics.RANSApplication.response_functions.vortex_shedding_frequency import VortexSheddingFrequency

def CreateResponseFunction(response_id, response_settings, model):
    response_type = response_settings["response_type"].GetString()

    if response_type == "lift_to_drag":
        return LiftToDrag(response_id, response_settings, model)
    elif response_type == "drag":
        return Drag(response_id, response_settings, model)
    elif response_type == "vortex_shedding_frequency":
        return VortexSheddingFrequency(response_id, response_settings, model)
    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: " +
                        "\n\t 'lift_to_drag'" +
                        "\n\t 'drag'" +
                        "\n\t 'vortex_shedding_frequency'")
