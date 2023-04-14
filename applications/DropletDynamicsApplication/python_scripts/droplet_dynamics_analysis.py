from sys import argv
from importlib import import_module

import KratosMultiphysics
# import KratosMultiphysics.FluidDynamicsApplication

from KratosMultiphysics.analysis_stage import AnalysisStage
# from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

class DropletDynamicsAnalysis(AnalysisStage): #FluidDynamicsAnalysis):

    # __init__ function can be defined specifically here to set deprecation warnings; using the template project parameters, there is no need for that

    def _CreateSolver(self):
        if not isinstance(self.model, KratosMultiphysics.Model):
            raise Exception("input is expected to be provided as a Kratos Model object")

        if not isinstance(self.project_parameters, KratosMultiphysics.Parameters):
            raise Exception("input is expected to be provided as a Kratos Parameters object")

        solver_settings = self.project_parameters["solver_settings"]

        module_full_name = 'KratosMultiphysics.DropletDynamicsApplication.droplet_dynamics_solver'
        solver = import_module(module_full_name).CreateSolver(self.model, solver_settings)

        return solver

    def _GetSimulationName(self):
        return "Droplet Dynamics Analysis"

    # def _CreateProcesses can be defined specifically here to give some hints or warnings regarding to the processes called within the project parameters

    def _GetOrderOfProcessesInitialization(self):
        # This is the duplicate of the same function defined in FluidDynamicsAnalysis. The aim is to provide essential data (e.g. gravity) first and prevent any unwanted modification to the BC
        return ["gravity",
                "initial_conditions_process_list",
                "boundary_conditions_process_list",
                "auxiliar_process_list"]

if __name__ == '__main__':
    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python droplet_dynamics_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python droplet_dynamics_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = DropletDynamicsAnalysis(model,parameters)
    simulation.Run()