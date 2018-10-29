from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
try:
    import KratosMultiphysics.ExternalSolversApplication
except ImportError:
    pass
from analysis_stage import AnalysisStage
import os
class PotentialFlowAnalysis(AnalysisStage):
        def __init__(self,model,parameters):
            super(PotentialFlowAnalysis,self).__init__(model,parameters)
            boundary_processes=self.project_parameters["processes"]["boundary_conditions_process_list"]
            case=self.project_parameters["solver_settings"]["problem_type"].GetString()
            defined=False         
            for i in range(0,boundary_processes.size()):
                python_module=boundary_processes[i]["python_module"].GetString()
                if python_module == "initialize_geometry":
                    defined=True
                    geometry_parameter=boundary_processes[i]["Parameters"]["geometry_parameter"].GetDouble()
                    skin_model_part_name=boundary_processes[i]["Parameters"]["skin_model_part_name"].GetString()
                    self.problem_name=case+"_"+self.project_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()+"_"+skin_model_part_name+"_"+str(geometry_parameter)
                if python_module == "compute_lift_level_set_process": 
                    self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"].AddEmptyValue("problem_name")   
                    self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"]["problem_name"].SetString(self.problem_name)
                if defined and python_module == 'define_wake_process':
                    self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"].AddEmptyValue("geometry_parameter")   
                    self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"]["geometry_parameter"].SetDouble(geometry_parameter)
                if not defined and python_module == "compute_lift_process":
                    self.problem_name=case+"_"+self.project_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()
                    self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"].AddEmptyValue("problem_name")   
                    self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"]["problem_name"].SetString(self.problem_name)
       
            if not os.path.exists('./gid_output/'):
                os.makedirs('./gid_output/')
            self.project_parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString('./gid_output/'+self.problem_name)
            
        def _CreateSolver(self):
            import potential_flow_solver 
            return potential_flow_solver.CreateSolver(self.model,self.project_parameters)
        def _GetOrderOfProcessesInitialization(self):
            return ["gravity",
                    "initial_conditions_process_list",
                    "boundary_conditions_process_list",
                    "auxiliar_process_list"]

        def _GetSimulationName(self):
            return "Potential Flow Analysis"
if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = PotentialFlowAnalysis(model,parameters)
    simulation.Run()