from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
try:
    import KratosMultiphysics.ExternalSolversApplication
except ImportError:
    pass

from fluid_dynamics_analysis import FluidDynamicsAnalysis

class ALEFluidAnalysis(FluidDynamicsAnalysis):
    '''Main script for fluid dynamics simulations using the navier_stokes family of python solvers.'''

    def _CreateSolver(self):
        import ale_fluid_solver
        return ale_fluid_solver.CreateSolver(self.model, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(ALEFluidAnalysis, self)._CreateProcesses(parameter_name, initialization_order)
        ## TODO call base-analysis, no deprecated settings are allowed!!!

    def _GetSimulationName(self):
        return "ALE Fluid Analysis"

if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python ale_fluid_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python ale_fluid_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = ALEFluidAnalysis(model,parameters)
    simulation.Run()
