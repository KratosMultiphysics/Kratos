# Importing Kratos
import KratosMultiphysics

# Importing the base class
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis 

class StructuralMechanicsAdjointDynamicAnalysis(StructuralMechanicsAnalysis):

    def __init__(self, model, project_parameters):
        solver_settings = project_parameters["solver_settings"]

        # Making sure that time step is negative
        if solver_settings["time_stepping"].Has("time_step"):
            if not solver_settings["time_stepping"]["time_step"].GetDouble() < 0:
                raise Exception("StructuralMechanicsAdjointDynamicAnalysis: " + '"time_step" in adjoint problem has to be negative!')

        # Setting start and end time
        if not project_parameters["problem_data"].Has("start_time"):
            project_parameters["problem_data"].AddEmptyValue("start_time")
            project_parameters["problem_data"]["start_time"].SetDouble( \
                            project_parameters["problem_data"]["start_step"].GetDouble() \
                                )    

        if not project_parameters["problem_data"].Has("end_time"):
            project_parameters["problem_data"].AddEmptyValue("end_time")
            project_parameters["problem_data"]["end_time"].SetDouble( \
                            project_parameters["problem_data"]["start_step"].GetDouble() + \
                            project_parameters["problem_data"]["nsteps"].GetInt()*solver_settings["time_stepping"]["time_step"].GetDouble()
                        )
        
        super().__init__(model, project_parameters)

    def KeepAdvancingSolutionLoop(self):
        """Note that the adjoint problem is solved in reverse time"""
        return self.time > self.end_time
    