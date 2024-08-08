# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.compute_drag_process import ComputeDragProcess

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeBodyFittedMomentProcess(model, settings["Parameters"])

class ComputeBodyFittedMomentProcess(ComputeDragProcess):
    """
    The specific implementation for the output of body fitted drag forces
    over obstacles in fluid dynamics problems.
    """
    def __init__(self, model, params: KratosMultiphysics.Parameters):
        self.moment_point = params["moment_point"].GetVector()
        params.RemoveValue("moment_point")
        super().__init__(model, params)

    def _GetFileHeader(self):
        header  = '# Body fitted moment for model part ' + self.params["model_part_name"].GetString() + ' for moment point at ' + str(self.moment_point) + '\n'
        header += '# Time Mx My Mz\n'
        return header

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo("ComputeBodyFittedMomentProcess","BODY FITTED MOMENT RESULTS:")
        KratosMultiphysics.Logger.PrintInfo("ComputeBodyFittedMomentProcess","Current time: " + result_msg)

    def _GetCorrespondingDragForce(self):
        moment = KratosMultiphysics.Array3([0.0, 0.0, 0.0])
        node: KratosMultiphysics.Node
        for node in self.model_part.Nodes:
            reaction = node.GetSolutionStepValue(KratosMultiphysics.REACTION)
            r = KratosMultiphysics.Array3([node.X, node.Y, node.Z]) - self.moment_point
            moment[0] += r[1] * reaction[2] - r[2] * reaction[1]
            moment[1] += r[2] * reaction[0] - r[0] * reaction[2]
            moment[2] += r[0] * reaction[1] - r[1] * reaction[0]
        return moment * -1.0