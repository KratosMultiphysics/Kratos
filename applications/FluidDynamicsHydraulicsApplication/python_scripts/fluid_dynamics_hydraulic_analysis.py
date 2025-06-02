
from sys import argv

import KratosMultiphysics as Kratos

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
import KratosMultiphysics.FluidDynamicsHydraulicsApplication as KratosHydraulics

from importlib import import_module


class FluidDynamicsHydraulicsAnalysis(FluidDynamicsAnalysis):
    '''Main script for fluid dynamics simulations using the navier_stokes family of python solvers.'''

    def _CreateSolver(self):
        solver_name = self.project_parameters["solver_settings"]["solver_type"].GetString()
        
        if solver_name == "two_fluid_hydraulic":
            solver_module_name = "navier_stokes_two_fluid_hydraulic_solver"
        elif solver_name == "two_fluid_hydraulic_fractional":
            solver_module_name = "navier_stokes_two_fluid_hydraulic_fractional_solver"

        module_full = 'KratosMultiphysics.FluidDynamicsHydraulicsApplication.' + solver_module_name
        solver = import_module(module_full).CreateSolver(self.model, self.project_parameters["solver_settings"])
        return solver
    
if __name__ == '__main__':
    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python fluid_dynamics_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fluid_dynamics_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = FluidDynamicsHydraulicsAnalysis(model, parameters)
    simulation.Run()
