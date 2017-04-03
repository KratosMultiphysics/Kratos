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

        self.main_model_part = ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
        self.main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, self.ProjectParameters["problem_data"]["domain_size"].GetInt())

        self.Model = {self.ProjectParameters["problem_data"]["model_part_name"].GetString() : self.main_model_part}

        self.solver_module = __import__(self.ProjectParameters["solver_settings"]["solver_type"].GetString())
        self.solver = self.solver_module.CreateSolver(self.main_model_part, self.ProjectParameters["solver_settings"])
        self.solver.AddVariables()
        self.solver.ImportModelPart()
        self.solver.AddDofs()

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

        for i in range(self.ProjectParameters["solver_settings"]["skin_parts"].size()):
            skin_part_name = self.ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
            self.Model.update({skin_part_name: self.main_model_part.GetSubModelPart(skin_part_name)})

        for i in range(self.ProjectParameters["solver_settings"]["no_skin_parts"].size()):
            no_skin_part_name = self.ProjectParameters["solver_settings"]["no_skin_parts"][i].GetString()
            self.Model.update({no_skin_part_name: self.main_model_part.GetSubModelPart(no_skin_part_name)})

        for i in range(self.ProjectParameters["initial_conditions_process_list"].size()):
            initial_cond_part_name = self.ProjectParameters["initial_conditions_process_list"][i]["Parameters"]["model_part_name"].GetString()
            self.Model.update({initial_cond_part_name: self.main_model_part.GetSubModelPart(initial_cond_part_name)})

        for i in range(self.ProjectParameters["gravity"].size()):
            gravity_part_name = self.ProjectParameters["gravity"][i]["Parameters"]["model_part_name"].GetString()
            self.Model.update({gravity_part_name: self.main_model_part.GetSubModelPart(gravity_part_name)})

        self.list_of_processes = process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses( self.ProjectParameters["gravity"] )
        self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses( self.ProjectParameters["initial_conditions_process_list"] )
        self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses( self.ProjectParameters["boundary_conditions_process_list"] )
        self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses( self.ProjectParameters["auxiliar_process_list"] )

        for process in self.list_of_processes:
            process.ExecuteInitialize()

        self.solver.Initialize()

        self.fluid_model_part = self.solver.GetComputingModelPart()


    def Solve(self):
        start_time = self.ProjectParameters["problem_data"]["start_time"].GetDouble()
        end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()

        time = start_time
        step = 0
        out = 0.0

        #self.gid_output.ExecuteBeforeSolutionLoop()

        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        while(time <= end_time):

            Dt = self.solver.ComputeDeltaTime()
            step += 1
            time = time + Dt
            self.main_model_part.CloneTimeStep(time)

            #if (self.parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0):
            #    print("")
            #    print("STEP = ", step)
            #    print("TIME = ", time)

            for process in self.list_of_processes:
                process.ExecuteInitializeSolutionStep()

            #self.gid_output.ExecuteInitializeSolutionStep()

            if(step >= 3):
                self.solver.SolveMeshMotion()

            for process in self.list_of_processes:
                process.ExecuteFinalizeSolutionStep()

            #self.gid_output.ExecuteFinalizeSolutionStep()

            for process in self.list_of_processes:
                process.ExecuteBeforeOutputStep()

            #if self.gid_output.IsOutputStep():
            #    self.gid_output.PrintOutput()

            for process in self.list_of_processes:
                process.ExecuteAfterOutputStep()

            out = out + Dt

        for process in self.list_of_processes:
            process.ExecuteFinalize()

        #self.gid_output.ExecuteFinalize()

if __name__ == '__main__':
    raise RuntimeError("This script should only be called from a test file.")
