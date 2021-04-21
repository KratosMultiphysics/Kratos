# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

from KratosMultiphysics.analysis_stage import AnalysisStage

from importlib import import_module

class ShallowWaterAnalysis(AnalysisStage):
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
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 shallow_water_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 shallow_water_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameter_file_name,'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()
    ShallowWaterAnalysis(model, parameters).Run()
