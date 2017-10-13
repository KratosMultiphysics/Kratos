from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.IGAStructuralMechanicsApplication import *

import os
import process_factory
#import NonConformant_OneSideMap

class KratosExecuteIGATest:

	def __init__(self, ProjectParameters):
		self.model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
		self.model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

		
		###TODO replace this "model" for real one once available in kratos core
		self.Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : self.model_part}


		#construct the IGASolver (main setting methods are located in the solver_module)
		solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
		self.IGASolver = solver_module.CreateSolver(self.model_part, ProjectParameters["solver_settings"])

		import read_materials_process
		read_materials_process.Factory(ProjectParameters,self.Model)
		
		self.IGASolver.AddVariables()
		self.IGASolver.ImportModelPart(ProjectParameters)
		self.IGASolver.AddDofs()
		
		for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
			part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
			self.Model.update({part_name: self.model_part.GetSubModelPart(part_name)})
		
		import process_factory
		self.list_of_processes = process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses( ProjectParameters["constraints_process_list"] )
		self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses( ProjectParameters["loads_process_list"] )
		if (ProjectParameters.Has("list_other_processes") == True):
			self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(ProjectParameters["list_other_processes"])
		if (ProjectParameters.Has("json_check_process") == True):
			self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(ProjectParameters["json_check_process"])
		if (ProjectParameters.Has("json_output_process") == True):
			self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(ProjectParameters["json_output_process"])

		for process in self.list_of_processes:
			process.ExecuteInitialize()
			
		self.IGASolver.Initialize()
		
	def Solve(self):
	
		self.model_part.ProcessInfo[TIME_STEPS] = 1
		self.model_part.CloneTimeStep(1.0)
		
		for process in self.list_of_processes:
			process.ExecuteBeforeSolutionLoop()
			
		for process in self.list_of_processes:
			process.ExecuteInitializeSolutionStep()	
			
		self.IGASolver.Solve()
		
		for process in self.list_of_processes:
			process.ExecuteFinalizeSolutionStep()

		for process in self.list_of_processes:
			process.ExecuteBeforeOutputStep()

		for process in self.list_of_processes:
			process.ExecuteAfterOutputStep()

		for process in self.list_of_processes:
			process.ExecuteFinalize()