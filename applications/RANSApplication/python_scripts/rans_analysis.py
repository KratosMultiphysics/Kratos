from sys import argv

import KratosMultiphysics as Kratos

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.RANSApplication.coupled_rans_solver import CoupledRANSSolver

class RANSAnalysis(FluidDynamicsAnalysis):
    '''Main script for fluid dynamics simulations using the navier_stokes family of python solvers.'''

    def __init__(self,model,parameters):
        super().__init__(model, parameters)

    def _CreateSolver(self):
        return CoupledRANSSolver(self.model, self.project_parameters["solver_settings"])

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
