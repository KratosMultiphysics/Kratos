from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.AdjointFluidApplication import *

import process_factory

class MainKratos:

    def __init__(self, ProjectParameters):

        self.ProjectParameters = ProjectParameters

        self.execute_solve = self.ProjectParameters["test_settings"]["execute_solve"].GetBool()

        self.main_model_part = ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
        self.main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, self.ProjectParameters["problem_data"]["domain_size"].GetInt())

        self.Model = {self.ProjectParameters["problem_data"]["model_part_name"].GetString() : self.main_model_part}

        ## Solver construction
        solver_module = __import__(self.ProjectParameters["solver_settings"]["solver_type"].GetString())
        self.solver = solver_module.CreateSolver(self.main_model_part, self.ProjectParameters["solver_settings"])
        self.solver.AddVariables()

        ## Read the model - note that SetBufferSize is done here
        self.solver.ImportModelPart()

        self.solver.AddDofs()

        #from gid_output_process import GiDOutputProcess
        #self.gid_output = GiDOutputProcess(self.solver.GetComputingModelPart(),
        #self.ProjectParameters["problem_data"]["problem_name"].GetString(),
        #self.ProjectParameters["output_configuration"])
        #self.gid_output.ExecuteInitialize()

        ## Get the list of the skin submodel parts in the object Model
        for i in range(self.ProjectParameters["solver_settings"]["skin_parts"].size()):
            skin_part_name = self.ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
            self.Model.update({skin_part_name: self.main_model_part.GetSubModelPart(skin_part_name)})

        ## Get the list of the initial conditions submodel parts in the object Model
        for i in range(self.ProjectParameters["initial_conditions_process_list"].size()):
            initial_cond_part_name = self.ProjectParameters["initial_conditions_process_list"][i]["Parameters"]["model_part_name"].GetString()
            self.Model.update({initial_cond_part_name: self.main_model_part.GetSubModelPart(initial_cond_part_name)})
    
        ## Get the gravity submodel part in the object Model
        for i in range(self.ProjectParameters["gravity"].size()):   
            gravity_part_name = self.ProjectParameters["gravity"][i]["Parameters"]["model_part_name"].GetString()
            self.Model.update({gravity_part_name: self.main_model_part.GetSubModelPart(gravity_part_name)})

        ## Obtain the list of the processes to be applied
        self.list_of_processes = process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses( self.ProjectParameters["initial_conditions_process_list"] )
        self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses( self.ProjectParameters["boundary_conditions_process_list"] )
        self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses( self.ProjectParameters["gravity"] )
        if (ProjectParameters.Has("list_other_processes") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(self.ProjectParameters["list_other_processes"])

        ##here all of the allocation of the strategies etc is done
        self.solver.Initialize()
        self.solver.SetEchoLevel(0)

        ## Processes initialization
        for process in self.list_of_processes:
            process.ExecuteInitialize()

        #TODO: think if there is a better way to do this
        self.fluid_model_part = self.solver.GetComputingModelPart()

    def Solve(self):
        ## Stepping and time settings
        Dt = self.ProjectParameters["problem_data"]["time_step"].GetDouble()
        nsteps = self.ProjectParameters["problem_data"]["nsteps"].GetInt()

        time = self.ProjectParameters["problem_data"]["start_step"].GetDouble()

        #self.gid_output.ExecuteBeforeSolutionLoop()

        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        self.main_model_part.CloneTimeStep(time)
        for step in range(1,nsteps+1):

            time = time + Dt
            self.main_model_part.CloneTimeStep(time)

            for process in self.list_of_processes:
                process.ExecuteInitializeSolutionStep()
            
            #self.gid_output.ExecuteInitializeSolutionStep()
                
            if self.execute_solve:
                self.solver.Solve()
        
            for process in self.list_of_processes:
                process.ExecuteFinalizeSolutionStep()
            
            #self.gid_output.ExecuteFinalizeSolutionStep()

            #if self.gid_output.IsOutputStep():
            #    self.gid_output.PrintOutput()
            
        for process in self.list_of_processes:
            process.ExecuteFinalize()
        
        #self.gid_output.ExecuteFinalize()

if __name__ == '__main__':
    raise RuntimeError("This script should only be called from a test file.")
