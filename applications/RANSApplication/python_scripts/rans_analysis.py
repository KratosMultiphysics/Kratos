from sys import argv

import KratosMultiphysics as Kratos

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.RANSApplication.coupled_rans_solver import CoupledRANSSolver

class RANSAnalysis(AnalysisStage):
    def __init__(self, model,parameters):
        """Main script for RANS fluid dynamics simulations using formulations via coupled_rans_solver

        Args:
            model (Kratos.Model): Model used in the simulation
            parameters (Kratos.Parameters): Parameters used in the simulation
        """
        super().__init__(model, parameters)

    def _CreateSolver(self):
        return CoupledRANSSolver(self.model, self.project_parameters["solver_settings"])

    def _GetOrderOfProcessesInitialization(self):
        # The list of processes will contain a list with each individual process already constructed (boundary conditions, initial conditions and gravity)
        # Note 1: gravity is constructed first. Outlet process might need its information.
        # Note 2: initial conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
        return ["gravity",
                "initial_conditions_process_list",
                "boundary_conditions_process_list",
                "auxiliar_process_list"]

    def _GetSimulationName(self):
        return "RANS Analysis"

    def KeepAdvancingSolutionLoop(self):
        """This function specifies the stopping criteria for breaking the solution loop
        It can be overridden by derived classes
        """
        return self.time < self.end_time and not self._GetSolver().IsConverged()

if __name__ == '__main__':
    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python rans_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python rans_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = RANSAnalysis(model,parameters)
    simulation.Run()
