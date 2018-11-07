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

    return ComputeBodyShearForceProcess(model, settings["Parameters"])

# computing shear force on given structure using the drag implementation
class ComputeBodyShearForceProcess(ComputeDragProcess):
    """
    The specific implementation for the output of body fitted drag forces
    over obstacles in fluid dynamics problems.
    """
    def _GetFileHeader(self):
        header  = '# Shear Force for the model part ' + self.params["model_part_name"].GetString() + '\n'
        header += '# Time Shear Force(N) \n'
        return header

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo("ComputeBodyShearForceProcess","BODY SHEAR RESULTS:")
        KratosMultiphysics.Logger.PrintInfo("ComputeBodyShearForceProcess","Current time: " + result_msg)

    def _GetCorrespondingDragForce(self):
        return KratosCFD.DragUtilities().CalculateBodyNormalForce(self.model_part) #its shear force, not pressure
    
    def ExecuteFinalizeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        drag_force = 0
        if((current_time >= self.interval[0]) and  (current_time < self.interval[1])):
            # Compute the drag force
            drag_force = self._GetCorrespondingDragForce()

            # Write the drag force values
            if (self.model_part.GetCommunicator().MyPID() == 0):
                if (self.print_drag_to_screen):
                    result_msg = str(current_time) + " x-drag: " + format(drag_force,self.format)
                    self._PrintToScreen(result_msg)

                if (self.write_drag_output_file):
                    self.output_file.write(str(current_time)+" "+format(drag_force,self.format)+"\n")
