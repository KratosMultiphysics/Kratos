# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.compute_base_moment_process import ComputeBaseMomentProcess

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeEmbeddedBaseMomentProcess(model, settings["Parameters"])

class ComputeEmbeddedBaseMomentProcess(ComputeBaseMomentProcess):
    """
    The specific implementation for the output of embedded fitted drag forces
    over obstacles in fluid dynamics problems.
    """
    def _GetFileHeader(self):
        header  = '# Embedded base moment for model part ' + self.params["model_part_name"].GetString() + '\n'
        header += '# Time Mx My Mz \n'
        return header

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo("ComputeEmbeddedBaseMomentProcess","EMBEDDED MOMENT RESULTS:")
        KratosMultiphysics.Logger.PrintInfo("ComputeEmbeddedBaseMomentProcess","Current time: " + result_msg)

    def _GetCorrespondingBaseMoment(self):
        return KratosCFD.DragUtilities().CalculateEmbeddedBaseMoment(self.model_part,self.reference_point)
