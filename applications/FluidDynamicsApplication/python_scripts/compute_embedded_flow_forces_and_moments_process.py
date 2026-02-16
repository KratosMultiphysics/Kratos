# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.compute_flow_forces_and_moments_process import ComputeFlowForcesAndMomentsProcess

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeEmbeddedFlowForcesAndMomentsProcess(model, settings["Parameters"])

class ComputeEmbeddedFlowForcesAndMomentsProcess(ComputeFlowForcesAndMomentsProcess):
    """
    The specific implementation for the output of embedded drag forces
    over obstacles in fluid dynamics problems.
    """
    def _GetFileHeader(self):
        header  = '# Embedded flow force and moment for model part ' + self.params["model_part_name"].GetString() + '\n'
        header += '# Time Fx Fy Fz Mx My Mz\n'
        return header

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo("ComputeEmbeddedFlowForcesAndMomentsProcess","EMBEDDED FLOW FORCES AND MOMENTS RESULTS:")
        KratosMultiphysics.Logger.PrintInfo("ComputeEmbeddedFlowForcesAndMomentsProcess","Current time: " + result_msg)

    def _GetCorrespondingFlowForcesAndMoments(self):
        return KratosCFD.FlowForcesAndMomentsUtilities().CalculateEmbeddedFlowForcesAndMoments(self.model_part, self.reference_point)