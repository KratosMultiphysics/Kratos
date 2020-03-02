# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.compute_drag_and_moment_process import ComputeDragAndMomentProcess

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeEmbeddedDragAndMomentProcess(model, settings["Parameters"])

class ComputeEmbeddedDragAndMomentProcess(ComputeDragAndMomentProcess):
    """
    The specific implementation for the output of body fitted drag forces
    over obstacles in fluid dynamics problems.
    """
    def _GetFileHeader(self):
        header  = '# Embedded drag and moment for model part ' + self.params["model_part_name"].GetString() + '\n'
        header += '# Time Fx Fy Fz Mx My Mz\n'
        return header

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo("ComputeEmbeddedDragAndMomentProcess","BODY FITTED DRAG RESULTS:")
        KratosMultiphysics.Logger.PrintInfo("ComputeEmbeddedDragAndMomentProcess","Current time: " + result_msg)

    def _GetCorrespondingDragForceAndMoment(self):
        return KratosCFD.DragAndMomentUtilities().CalculateEmbeddedDragAndMoment(self.model_part,self.reference_point)