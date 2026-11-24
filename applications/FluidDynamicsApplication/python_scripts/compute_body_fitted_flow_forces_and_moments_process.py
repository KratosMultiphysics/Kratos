# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.compute_flow_forces_and_moments_process import ComputeFlowForcesAndMomentsProcess

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeBodyFittedFlowForcesAndMomentsProcess(model, settings["Parameters"])

class ComputeBodyFittedFlowForcesAndMomentsProcess(ComputeFlowForcesAndMomentsProcess):
    """
    The specific implementation for the output of body fitted flow forces and moments
    over obstacles in fluid dynamics problems.
    """
    def _GetFileHeader(self):
        header  = '# Body fitted flow force and moment for model part '
        header += self.params["model_part_name"].GetString() + '\n'

        if self.write_flow_moments_output_file:
            header += '# Time Fx Fy Fz Mx My Mz\n'
        else:
            header += '# Time Fx Fy Fz\n'

        return header

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo("ComputeBodyFittedFlowForcesAndMomentsProcess","BODY FITTED FLOW FORCES AND MOMENTS RESULTS:")
        KratosMultiphysics.Logger.PrintInfo("ComputeBodyFittedFlowForcesAndMomentsProcess","Current time: " + result_msg)

    def _GetCorrespondingFlowForcesAndMoments(self):
        return KratosCFD.FlowForcesAndMomentsUtilities().CalculateBodyFittedFlowForcesAndMoments(self.model_part, self.reference_point)