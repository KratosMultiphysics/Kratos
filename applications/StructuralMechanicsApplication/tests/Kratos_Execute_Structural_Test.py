from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

import os
import process_factory

class Kratos_Execute_Test:

    def __init__(self, ProjectParameters):

        self.ProjectParameters = ProjectParameters

        self.main_model_part = KratosMultiphysics.ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.ProjectParameters["problem_data"]["domain_size"].GetInt())

        self.Model = {self.ProjectParameters["problem_data"]["model_part_name"].GetString(): self.main_model_part}

        # Construct the solver (main setting methods are located in the solver_module)
        import python_solvers_wrapper_structural
        self.solver = python_solvers_wrapper_structural.CreateSolver(self.main_model_part, self.ProjectParameters)

        # Add variables (always before importing the model part) (it must be integrated in the ImportModelPart)
        # If we integrate it in the model part we cannot use combined solvers
        self.solver.AddVariables()

        # Read model_part (note: the buffer_size is set here) (restart can be read here)
        self.solver.ImportModelPart()

        # Add dofs (always after importing the model part) (it must be integrated in the ImportModelPart)
        # If we integrate it in the model part we cannot use combined solvers
        self.solver.AddDofs()

        # Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
        # #Get the list of the submodel part in the object Model
        for i in range(self.ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
            part_name = self.ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
            self.Model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})

        # Obtain the list of the processes to be applied
        self.list_of_processes = process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(self.ProjectParameters["constraints_process_list"])
        self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(self.ProjectParameters["loads_process_list"])
        if (ProjectParameters.Has("list_other_processes") == True): 
            self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(self.ProjectParameters["list_other_processes"])
        if (ProjectParameters.Has("json_check_process") == True): 
            self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(self.ProjectParameters["json_check_process"]) 
        if (ProjectParameters.Has("check_analytic_results_process") == True): 
            self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(self.ProjectParameters["check_analytic_results_process"]) 
        if (ProjectParameters.Has("json_output_process") == True): 
            self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(self.ProjectParameters["json_output_process"]) 

        for process in self.list_of_processes:
            process.ExecuteInitialize()

        # ### START SOLUTION ####

        self.computing_model_part = self.solver.GetComputingModelPart()

        # ### Output settings start ####
        self.problem_path = os.getcwd()
        self.problem_name = self.ProjectParameters["problem_data"]["problem_name"].GetString()

        # ### Output settings start ####
        self.output_post = ProjectParameters.Has("output_configuration")
        if (self.output_post == True):
            from gid_output_process import GiDOutputProcess
            output_settings = ProjectParameters["output_configuration"]
            self.gid_output = GiDOutputProcess(self.computing_model_part,
                                               self.problem_name,
                                               output_settings)
            self.gid_output.ExecuteInitialize()
            
        # Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
        self.solver.Initialize()
        self.solver.SetEchoLevel(0) # Avoid to print anything 
        
        if (self.output_post == True):
            self.gid_output.ExecuteBeforeSolutionLoop()

    def Solve(self):
        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        # #Stepping and time settings (get from process info or solving info)
        # Delta time
        delta_time = self.ProjectParameters["problem_data"]["time_step"].GetDouble()
        # Start step
        self.main_model_part.ProcessInfo[KratosMultiphysics.TIME_STEPS] = 0
        # Start time
        time = self.ProjectParameters["problem_data"]["start_time"].GetDouble()
        # End time
        end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()

        # Solving the problem (time integration)
        while(time <= end_time):
            time = time + delta_time
            self.main_model_part.ProcessInfo[KratosMultiphysics.TIME_STEPS] += 1
            self.main_model_part.CloneTimeStep(time)

            for process in self.list_of_processes:
                process.ExecuteInitializeSolutionStep()
                
            if (self.output_post == True):
                self.gid_output.ExecuteInitializeSolutionStep()
                        
            self.solver.Solve()
            
            if (self.output_post == True):
                self.gid_output.ExecuteFinalizeSolutionStep()

            for process in self.list_of_processes:
                process.ExecuteFinalizeSolutionStep()

            for process in self.list_of_processes:
                process.ExecuteBeforeOutputStep()

            for process in self.list_of_processes:
                process.ExecuteAfterOutputStep()

            if (self.output_post == True):
                if self.gid_output.IsOutputStep():
                    self.gid_output.PrintOutput()

        if (self.output_post == True):
            self.gid_output.ExecuteFinalize()

        for process in self.list_of_processes:
            process.ExecuteFinalize()
