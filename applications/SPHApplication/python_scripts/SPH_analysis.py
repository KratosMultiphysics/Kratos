import json

import KratosMultiphysics as KM 
import KratosMultiphysics.SPHApplication as SPH

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.SPHApplication import python_solvers_wrapper_sph as sph_solvers

class SPHAnalysis(AnalysisStage):
	def __init__(self, model, project_parameters):

		solver_settings = project_parameters["solver_settings"]

		if solver_settings.Has("domain_size") and project_parameters["problem_data"].Has("domain_size"):
			raise Exception("SPHAnalysis: " + '"domain_size" defined both in "problem_data" and "solver_settings"!')

		if solver_settings.Has("model_part_name") and project_parameters["problem_data"].Has("model_part_name"):
			raise Exception("SPHAnalysis: " + '"model_part_name" defined both in problem_data" and "solver_settings"!')

		if solver_settings.Has("time_stepping") and project_parameters["problem_data"].Has("time_Step"):
			raise Exception("SPHAnalysis: " + '"time_stepping" defined both in "problem_data" and "solver_settings"!')

		if not solver_settings.Has("time_stepping"):
			raise Exception("SPHAnalysis: Using the old way to pass the time_step, this was removed!")

		if not solver_settings.Has("domain_size"):
			raise Exception("SPHAnalysis: Using the old way to pass the domain_size, this was removed!")

		try:
			with open("SPHProcesses.json", "r") as process_file:
				process_data = process_file.read()
				self.sph_process_parameters = KM.Parameters(process_data)
		except FileNotFoundError as exc:
			raise Exception("SPHProcesses.json not found") from exc

		solver_settings = project_parameters["solver_settings"]
        
		super().__init__(model, project_parameters)

	def _CreateSolver(self):
		"""Create the solver"""
		KM.Logger.PrintInfo("::[SPHAnalysis]:: ", "Creating SPH solver")
		return sph_solvers.CreateSolver(self.model, self.project_parameters)

	def _CreateProcesses(self, parameter_name, initialization_order):
		list_of_processes = super()._CreateProcesses(parameter_name, initialization_order)

		model_part = self.model["Structure"]

		if parameter_name != "processes":
			return list_of_processes

		# Adding SPH specific processes (to change parameters, see SPHProcesses.json)
		neighbours_search_process = SPH.NeighboursSearchProcess(model_part, self.sph_process_parameters["neighbours_search"])
		kernel_correction_process = SPH.ComputeKernelCorrectionProcess(model_part, self.sph_process_parameters["kernel_correction"])
		volume_process = SPH.ComputeVolumeProcess(self.model["Structure.Triangulation"], self.sph_process_parameters["volume"])

		custom_processes = [neighbours_search_process, kernel_correction_process, volume_process]

		return custom_processes + list_of_processes

	def _GetSimulationName(self):
		return "::[SPH Simulation]:: "

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()
    simulation = SPHAnalysis(model, parameters)
    simulation.Run()
