# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

from KratosMultiphysics.analysis_stage_with_solver import AnalysisStageWithSolver

from importlib import import_module

class ShallowWaterAnalysis(AnalysisStageWithSolver):
    ''' Main script for shallow water simulations '''

    def _CreateSolver(self):
        python_module_name = "KratosMultiphysics.ShallowWaterApplication"
        full_module_name = python_module_name + "." + self.project_parameters["solver_settings"]["solver_type"].GetString()
        solver_module = import_module(full_module_name)
        solver = solver_module.CreateSolver(self.model, self.project_parameters["solver_settings"])
        return solver

    def _GetOrderOfProcessesInitialization(self):
        return ["topography_process_list",
                "initial_conditions_process_list",
                "boundary_conditions_process_list"]

    def _GetSimulationName(self):
        return "Shallow Water Analysis"

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("parameters_file_name", nargs="?", default="ProjectParameters.json")
    args = parser.parse_args()

    with open(args.parameters_file_name, 'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()
    ShallowWaterAnalysis(model, parameters).Run()
