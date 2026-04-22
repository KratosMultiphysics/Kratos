# Importing the Kratos Library
import KratosMultiphysics

# Import new process
from KratosMultiphysics.FluidDynamicsApplication import compute_flow_forces_and_moments_process


def Factory(settings, model):
    if type(settings) != KratosMultiphysics.Parameters:
        raise Exception(
            "Expected input shall be a Parameters object, encapsulating a json string")

    return ComputeDragProcess(model, settings["Parameters"])


class ComputeDragProcess(
        compute_flow_forces_and_moments_process.ComputeFlowForcesAndMomentsProcess):
    """
    DEPRECATED (backward compatibility):

    Legacy wrapper around ComputeFlowForcesAndMomentsProcess.
    This class exists only to keep old input files working.
    All functionality is implemented in the base class.
    """

    def __init__(self, model, params):
        # Deprecation warning
        KratosMultiphysics.Logger.PrintWarning(
            "ComputeDragProcess",
            'DEPRECATED: "ComputeDragProcess" is deprecated and will be removed in a future release. '
            'Please use "ComputeFlowForcesAndMomentsProcess" instead.')

        # Initialize base class
        super().__init__(model, params)
