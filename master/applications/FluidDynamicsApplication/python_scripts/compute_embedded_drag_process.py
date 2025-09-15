# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.compute_drag_process import ComputeDragProcess

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeEmbeddedDragProcess(model, settings["Parameters"])

class ComputeEmbeddedDragProcess(ComputeDragProcess):
    """
    The specific implementation for the output of embedded drag forces
    over obstacles in fluid dynamics problems.
    """
    def _GetFileHeader(self):
        header  = '# Embedded drag for model part ' + self.params["model_part_name"].GetString() + '\n'
        header += '# Time Fx Fy Fz\n'
        return header

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo("ComputeEmbeddedDragProcess","EMBEDDED DRAG RESULTS:")
        KratosMultiphysics.Logger.PrintInfo("ComputeEmbeddedDragProcess","Current time: " + result_msg)

    def _GetCorrespondingDragForce(self):
        return KratosCFD.DragUtilities().CalculateEmbeddedDrag(self.model_part)