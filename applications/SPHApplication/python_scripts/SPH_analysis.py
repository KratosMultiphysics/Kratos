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

		params = KM.Parameters("""
		{
			"model_part_name" : "Structure",
			"coefficient"     : 0.5
		}
		""");
		neighbours_process = SPH.NeighboursSearchProcess(model_part, params)

		params = KM.Parameters("""
		{
			"model_part_name" : "Structure",
			"controls"     : true,
			"tol" : 1e-10
		}
		""");
		kernel_correction_process = SPH.ComputeKernelCorrectionProcess(model_part, params)
		volume_process = SPH.ComputeVolumeProcess(self.model["Structure.Triangulation"], KM.Parameters("{}"))

		## Il problema al momento è scrivere in modo più compatto questi input e assegnare quelli di default,
		## default non assegato perchè non viene richiamata la factory

		list_of_processes.insert(0, volume_process)
		list_of_processes.insert(1, neighbours_process)
		list_of_processes.insert(2, kernel_correction_process)

		return list_of_processes

	def _GetSimulationName(self):
		return "::[SPH Simulation]:: "

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()
    simulation = SPHAnalysis(model, parameters)
    simulation.Run()
