from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *

import json

import process_factory

class Kratos_Execute_Test:
  
  def __init__(self, ProjectParameters):
    self.ProjectParameters = ProjectParameters
    self.main_model_part = ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
    self.main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, self.ProjectParameters["problem_data"]["domain_size"].GetInt())

    ###TODO replace this "model" for real one once available in kratos core
    self.Model = {self.ProjectParameters["problem_data"]["model_part_name"].GetString() : self.main_model_part}

    #construct the solver (main setting methods are located in the solver_module)
    solver_module = __import__(self.ProjectParameters["solver_settings"]["solver_type"].GetString())
    self.solver = solver_module.CreateSolver(self.main_model_part, self.ProjectParameters["solver_settings"])

    #add variables (always before importing the model part) (it must be integrated in the ImportModelPart)
    # if we integrate it in the model part we cannot use combined solvers
    self.solver.AddVariables()
    ### Temporal 
    #self.main_model_part.AddNodalSolutionStepVariable(ALPHA_EAS)

    #read model_part (note: the buffer_size is set here) (restart can be read here)
    self.solver.ImportModelPart()

    #add dofs (always after importing the model part) (it must be integrated in the ImportModelPart)
    # if we integrate it in the model part we cannot use combined solvers
    self.solver.AddDofs()

    #build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
    ##TODO: replace MODEL for the Kratos one ASAP
    ##get the list of the submodel part in the object Model
    for i in range(self.ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
      part_name = self.ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
      self.Model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})

    #obtain the list of the processes to be applied
    self.list_of_processes = process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses( self.ProjectParameters["constraints_process_list"] )
    self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses( self.ProjectParameters["loads_process_list"] )
    self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses( self.ProjectParameters["list_other_processes"] )

    for process in self.list_of_processes:
      process.ExecuteInitialize()
    
    #### START SOLUTION ####

    #TODO: think if there is a better way to do this
    self.computing_model_part = self.solver.GetComputeModelPart()

    #### output settings start ####

    self.problem_path = os.getcwd()
    self.problem_name = self.ProjectParameters["problem_data"]["problem_name"].GetString()

    #### output settings start ####

    ## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
    self.solver.Initialize()

  def Solve(self):
    for process in self.list_of_processes:
      process.ExecuteBeforeSolutionLoop()
      
    ## Stepping and time settings (get from process info or solving info)
    #delta time
    delta_time = self.ProjectParameters["problem_data"]["time_step"].GetDouble()
    #start step
    step       = 0
    #start time
    time       = self.ProjectParameters["problem_data"]["start_time"].GetDouble()
    #end time
    end_time   = self.ProjectParameters["problem_data"]["end_time"].GetDouble()
    
    # solving the problem (time integration)
    while(time <= end_time):
      time = time + delta_time
      step = step + 1
      self.main_model_part.CloneTimeStep(time)
      
      for process in self.list_of_processes:
        process.ExecuteInitializeSolutionStep()
        
      self.solver.Solve()

      for process in self.list_of_processes:
        process.ExecuteFinalizeSolutionStep()

      for process in self.list_of_processes:
        process.ExecuteBeforeOutputStep()

      for process in self.list_of_processes:
        process.ExecuteAfterOutputStep()

      for process in self.list_of_processes:
        process.ExecuteFinalize()

  #def CheckResults(self, Results):
      
