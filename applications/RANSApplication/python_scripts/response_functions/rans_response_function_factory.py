# importing the Kratos Library
from KratosMultiphysics.RANSApplication.response_functions.lift_to_drag import LiftToDrag
from KratosMultiphysics.RANSApplication.response_functions.drag import Drag
from KratosMultiphysics.RANSApplication.response_functions.drag_frequency_max_amplitude import DragFrequencyMaxAmplitude
from KratosMultiphysics.RANSApplication.response_functions.transient_drag_steady_adjoint import TransientDragSteadyAdjoint
from KratosMultiphysics.RANSApplication.response_functions.domain_integrated_3d_vector_magnitude_square_power_mean import DomainIntegrated3DVectorMagnitudeSquarePowerMean


def CreateResponseFunction(response_id, response_settings, model):
    response_type = response_settings["response_type"].GetString()

    if response_type == "lift_to_drag":
        return LiftToDrag(response_id, response_settings, model)
    elif response_type == "drag":
        return Drag(response_id, response_settings, model)
    elif response_type == "drag_frequency_max_amplitude":
        return DragFrequencyMaxAmplitude(response_id, response_settings, model)
    elif response_type == "transient_drag_steady_adjoint":
        return TransientDragSteadyAdjoint(response_id, response_settings, model)
    elif response_type == "domain_integrated_3d_vector_magnitude_square_power_mean":
        return DomainIntegrated3DVectorMagnitudeSquarePowerMean(response_id, response_settings, model)
    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: " +
                        "\n\t 'lift_to_drag'" +
                        "\n\t 'drag'" +
                        "\n\t 'drag_frequency_max_amplitude'" +
                        "\n\t 'transient_drag_steady_adjoint'" +
                        "\n\t 'domain_integrated_3d_vector_magnitude_square_power_mean'")
