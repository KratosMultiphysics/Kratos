# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from compute_drag_process import ComputeDragProcess

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeTimeAveragedBodyFittedDragProcess(model, settings["Parameters"])

class ComputeTimeAveragedBodyFittedDragProcess(ComputeDragProcess):
    """
    The specific implementation for the output of body fitted drag forces
    over obstacles in fluid dynamics problems.
    """
    def _GetFileHeader(self):
        header  = '# Time average body fitted drag for model part ' + self.params["model_part_name"].GetString() + '\n'
        header += '# Time Fx Fy Fz\n'
        return header

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo("ComputeTimeAveragedBodyFittedDragProcess","TIME AVERAFED BODY FITTED DRAG RESULTS:")
        KratosMultiphysics.Logger.PrintInfo("ComputeTimeAveragedBodyFittedDragProcess","Current time: " + result_msg)

    def _GetCorrespondingDragForce(self):
        return KratosCFD.DragUtilities().CalculateTimeAveragedBodyFittedDrag(self.model_part)