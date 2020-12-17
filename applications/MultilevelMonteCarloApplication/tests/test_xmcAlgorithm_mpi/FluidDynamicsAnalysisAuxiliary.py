# Import Python libraries
import numpy as np
import time

# Import Kratos
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.ExaquteSandboxApplication
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
import KratosMultiphysics.kratos_utilities as kratos_utilities


class FluidDynamicsAnalysisAuxiliary(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters):
        super().__init__(model,project_parameters)
        self.drag_force_vector = np.zeros([0,4])
        # set model part of interest
        self.interest_model_part = "MainModelPart.NoSlip2D_No_Slip_Auto1"

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
        # store current force x and moment z for updating the time power sums
        self.current_force_x = drag_force_vector[0]

    def Finalize(self):
        super().Finalize()
        burnin_time = self.project_parameters["problem_data"]["burnin_time"].GetDouble()
        force_x_no_washout = [self.drag_force_vector[i,1] for i in range (1,len(self.drag_force_vector[:,0])) if (self.drag_force_vector[i-1,0] >= burnin_time)] # not even check time step 0
        self.mean_force_x = np.mean(force_x_no_washout)
        print("[INFO] Final averaged drag value", self.mean_force_x)

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python fluid_dynamics_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fluid_dynamics_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "problem_settings/ProjectParametersRectangularCylinder2D_Fractional.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()

    ini_time = time.time()
    simulation = FluidDynamicsAnalysisProblemZero(model, parameters)
    simulation.Run()
    print("[TIMER] Total analysis time:", time.time()-ini_time)
