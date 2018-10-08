from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
try:
    import KratosMultiphysics.ExternalSolversApplication
except ImportError:
    pass
from analysis_stage import AnalysisStage

class PotentialFlowAnalysis(AnalysisStage):
        def __init__(self,model,parameters):
            super(PotentialFlowAnalysis,self).__init__(model,parameters)
            boundary_processes=self.project_parameters["processes"]["boundary_conditions_process_list"]
            defined=False            
            for i in range(0,boundary_processes.size()):
                python_module=boundary_processes[i]["python_module"].GetString()
                if python_module == "initialize_geometry":
                    defined=True
                    geometry_parameter=boundary_processes[i]["Parameters"]["geometry_parameter"].GetDouble()
                    initial_case=boundary_processes[i]["Parameters"]["initial_case"].GetString()
                    self.problem_name=self.project_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()+"_"+initial_case+"_"+str(geometry_parameter)
                if python_module == "compute_lift_level_set_process": 
                    self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"].AddEmptyValue("problem_name")   
                    self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"]["problem_name"].SetString(self.problem_name)
                if defined and python_module == 'define_wake_process':
                    self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"].AddEmptyValue("geometry_parameter")   
                    self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"]["geometry_parameter"].SetDouble(geometry_parameter)
                if not defined and python_module == "compute_lift_process":
                    self.problem_name=self.project_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()
                    self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"].AddEmptyValue("problem_name")   
                    self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"]["problem_name"].SetString(self.problem_name)
            
        def _CreateSolver(self):
            import potential_flow_solver 
            return potential_flow_solver.CreateSolver(self.model,self.project_parameters)
        def _GetOrderOfProcessesInitialization(self):
            return ["gravity",
                    "initial_conditions_process_list",
                    "boundary_conditions_process_list",
                    "auxiliar_process_list"]
        def _SetUpGiDOutput(self):
            '''Initialize a GiD output instance'''
            if self.parallel_type == "OpenMP":
                from gid_output_process import GiDOutputProcess as OutputProcess
            elif self.parallel_type == "MPI":
                from gid_output_process_mpi import GiDOutputProcessMPI as OutputProcess

            output = OutputProcess(self._GetSolver().GetComputingModelPart(),
                                    './gid_output/'+self.problem_name,
                                    self.project_parameters["output_configuration"])
            return output
        def RunSolutionLoop(self):
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()
            self.gid_output= self._SetUpGiDOutput()
            self.gid_output.ExecuteInitialize()
            self.gid_output.ExecuteBeforeSolutionLoop()
            self.gid_output.ExecuteInitializeSolutionStep()
            self.gid_output.PrintOutput()
            self.gid_output.ExecuteFinalize()
        # def _GetListOfOutputProcesses(self):
        #     self._list_of_output_processes=super(PotentialFlowAnalysis,self)._GetListOfOutputProcesses()
        #     if self.project_parameters.Has("output_configuration"):                
        #         gid_output= self._SetUpGiDOutput()
        #         self._list_of_output_processes += [gid_output,]
        #     return self._list_of_output_processes

        def _CreateProcesses(self, parameter_name, initialization_order):
            list_of_processes = super(PotentialFlowAnalysis, self)._CreateProcesses(parameter_name, initialization_order)
            # if self.project_parameters.Has("output_configuration"):
                
            #     self.gid_output= self._SetUpGiDOutput()
            #     list_of_processes += [self.gid_output,]


            return list_of_processes
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