# Import Python libraries
import sys
import time
import numpy as np

# IMport Kratos
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis


class FluidDynamicsAnalysisMC(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters):
        super().__init__(model,project_parameters)
        self.drag_force_vector = np.zeros([0,4])
        self.default_time_step = self.project_parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()
        # set model part of interest
        self.interest_model_part = "FluidModelPart.NoSlip3D_structure"

    def ModifyInitialProperties(self):
        """
        function changing process settings
        input:  self: an instance of the class
        """
        super().ModifyInitialProperties()
        for aux_process in self.project_parameters["processes"]["auxiliar_process_list"]:
            if aux_process["python_module"].GetString() == "temporal_statistics_process":
                aux_process["Parameters"]["statistics_start_point_control_value"].SetDouble(self.project_parameters["problem_data"]["burnin_time"].GetDouble())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        # compute drag force
        drag_force_vector = KratosMultiphysics.FluidDynamicsApplication.DragUtilities().CalculateBodyFittedDrag(self.model.GetModelPart(self.interest_model_part))
        drag_force = [self.time,drag_force_vector[0],drag_force_vector[1],drag_force_vector[2]]
        self.drag_force_vector = np.vstack((self.drag_force_vector,drag_force))
        # store current force x for updating the time power sums
        self.current_drag_force_x = drag_force_vector[0]
        # set larger-smaller time step
        if (self.time >= self.project_parameters["problem_data"]["burnin_time"].GetDouble()): # burning time
            self.project_parameters["solver_settings"]["time_stepping"]["time_step"].SetDouble(self.default_time_step)
        else:
            self.project_parameters["solver_settings"]["time_stepping"]["time_step"].SetDouble(1.0*self.default_time_step)
            # self.project_parameters["solver_settings"]["time_stepping"]["time_step"].SetDouble(2.5*self.default_time_step)

    def Finalize(self):
        super().Finalize()
        burnin_time = self.project_parameters["problem_data"]["burnin_time"].GetDouble()
        drag_force_x_post_burnin = [self.drag_force_vector[i,1] for i in range (1,len(self.drag_force_vector[:,0])) if (self.drag_force_vector[i-1,0] >= burnin_time)] # not even check time step 0
        self.mean_drag_force_x = np.mean(drag_force_x_post_burnin)

if __name__ == "__main__":

    if len(sys.argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = sys.argv[1]
    else: # using default name
        parameter_file_name = "problem_settings/ProjectParametersCAARC_MLMC_steadyInlet.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = FluidDynamicsAnalysisMC(model,parameters)
    simulation.Run()
