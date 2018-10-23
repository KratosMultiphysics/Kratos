# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from compute_drag_process import ComputeDragProcess

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeBodyFittedDragProcess(model, settings["Parameters"])

class ComputeBodyFittedDragProcess(ComputeDragProcess):
    """
    The specific implementation for the output of body fitted drag forces
    over obstacles in fluid dynamics problems.
    """
    def _GetFileHeader(self):
        header  = '# Body fitted drag for model part ' + self.params["model_part_name"].GetString() + '\n'
        header += '# Time Fx Fy Fz\n'
        return header

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo("ComputeBodyFittedDragProcess","BODY FITTED DRAG RESULTS:")
        KratosMultiphysics.Logger.PrintInfo("ComputeBodyFittedDragProcess","Current time: " + result_msg)

    def _GetCorrespondingDragForce(self):
        return KratosCFD.DragUtilities().CalculateBodyFittedDrag(self.model_part)