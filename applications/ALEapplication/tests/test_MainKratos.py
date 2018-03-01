from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

import process_factory

class MainKratos:

    def __init__(self, ProjectParameters):

        self.ProjectParameters = ProjectParameters

        self.parallel_type = self.ProjectParameters["problem_data"]["parallel_type"].GetString()

        if (self.parallel_type == "MPI"):
            import KratosMultiphysics.mpi as KratosMPI
            import KratosMultiphysics.TrilinosApplication as KratosTrilinos

        self.main_model_part_name = ProjectParameters["problem_data"]["model_part_name"].GetString()
        self.main_model_part = ModelPart(self.main_model_part_name)
        self.main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, self.ProjectParameters["problem_data"]["domain_size"].GetInt())


        import python_solvers_wrapper_mesh_motion
        self.solver = python_solvers_wrapper_mesh_motion.CreateSolver(self.main_model_part, self.ProjectParameters)
        self.solver.AddVariables()
        self.solver.ImportModelPart()
        self.solver.AddDofs()

##---------------------
        #if (self.parallel_type == "OpenMP"):
        #    from gid_output_process import GiDOutputProcess
        #    self.gid_output = GiDOutputProcess(self.solver.GetComputingModelPart(),
        #                                       self.ProjectParameters["problem_data"]["problem_name"].GetString(),
        #                                       self.ProjectParameters["output_configuration"])
        #elif (self.parallel_type == "MPI"):
        #    from gid_output_process_mpi import GiDOutputProcessMPI
        #    self.gid_output = GiDOutputProcessMPI(self.solver.GetComputingModelPart(),
        #                                          self.ProjectParameters["problem_data"]["problem_name"].GetString(),
        #                                          self.ProjectParameters["output_configuration"])
        #self.gid_output.ExecuteInitialize()
##---------------------

        self.Model = Model()
        self.Model.AddModelPart(self.main_model_part)

        self.list_of_processes = process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(self.ProjectParameters["boundary_conditions_process_list"])
        if (ProjectParameters.Has("list_other_processes") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(self.ProjectParameters["list_other_processes"])
        if (ProjectParameters.Has("json_check_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(self.ProjectParameters["json_check_process"])
        if (ProjectParameters.Has("check_analytic_results_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(self.ProjectParameters["check_analytic_results_process"])
        if (ProjectParameters.Has("json_output_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(self.ProjectParameters["json_output_process"])
        #print(self.list_of_processes)
        for process in self.list_of_processes:
            process.ExecuteInitialize()

        self.solver.Initialize()


    def Solve(self):
        start_time = self.ProjectParameters["problem_data"]["start_time"].GetDouble()
        end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()

        time = start_time
        step = 0
        out = 0.0
##---------------------
        #self.gid_output.ExecuteBeforeSolutionLoop()
##---------------------
        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()


        ## Stepping and time settings
        delta_time = self.ProjectParameters["problem_data"]["time_step"].GetDouble()
        start_time = self.ProjectParameters["problem_data"]["start_time"].GetDouble()
        end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()

        time = start_time
        #self.main_model_part.ProcessInfo[STEP] = 0

        while(time <= end_time):
            step += 1
            time = time + delta_time
            self.main_model_part.CloneTimeStep(time)

            #if (self.parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0):
            #    print("")
            #    print("STEP = ", step)
            #    print("TIME = ", time)

            for process in self.list_of_processes:
                process.ExecuteInitializeSolutionStep()
##---------------------
            #self.gid_output.ExecuteInitializeSolutionStep()
##---------------------

            self.solver.Solve()

            for process in self.list_of_processes:
                process.ExecuteFinalizeSolutionStep()
##---------------------
            #self.gid_output.ExecuteFinalizeSolutionStep()
##---------------------
            for process in self.list_of_processes:
                process.ExecuteBeforeOutputStep()
##---------------------
            #if self.gid_output.IsOutputStep():
            #    self.gid_output.PrintOutput()
##---------------------
            for process in self.list_of_processes:
                process.ExecuteAfterOutputStep()

            out = out + delta_time

        for process in self.list_of_processes:
            process.ExecuteFinalize()
##---------------------
        #self.gid_output.ExecuteFinalize()
##---------------------
if __name__ == '__main__':
    raise RuntimeError("This script should only be called from a test file.")
