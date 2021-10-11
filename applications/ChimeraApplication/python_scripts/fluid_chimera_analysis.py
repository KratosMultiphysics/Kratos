# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import absolute_import, division

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
import KratosMultiphysics.ChimeraApplication as KratosChimera
import KratosMultiphysics.ChimeraApplication.python_solvers_wrapper_fluid_chimera


class FluidChimeraAnalysis(FluidDynamicsAnalysis):
    '''
    Main script for fluid chimera simulations using the navier stokes family of python solvers.
    '''

    def __init__(self, model, parameters):
        super().__init__(model, parameters)

    def _CreateSolver(self):
        return KratosChimera.python_solvers_wrapper_fluid_chimera.CreateSolver(self.model, self.project_parameters)

    def _GetSimulationName(self):
        solver_type = self.project_parameters["solver_settings"]["solver_type"].GetString()
        if (solver_type == "monolithic_rans_chimera" or solver_type == "MonolithicRANSChimera"):
            return "RANS Chimera Analysis"
        else:
            return "Fluid Chimera Analysis"

    def KeepAdvancingSolutionLoop(self):
        """This function specfies the stopping criteria for breaking the solution loop.
        It can be overridden by derived class"""

        solver_type = self.project_parameters["solver_settings"]["solver_type"].GetString()
        if (solver_type == "monolithic_rans_chimera" or solver_type == "MonolithicRANSChimera"):
            # for RANSChimera: IsConverged() returns false if convergence occurs in ramp-up interval
            return (self.time < self.end_time) and (not self._GetSolver().IsConverged())
        else:
            # for FluidDynamicsChimera
            return (self.time < self.end_time)
            


if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg = 'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python fluid_chimera_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fluid_chimera_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2:  # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else:  # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name, 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = FluidChimeraAnalysis(model, parameters)
    simulation.Run()
