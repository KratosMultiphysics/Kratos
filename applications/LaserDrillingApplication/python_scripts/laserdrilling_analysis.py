
# Importing Kratos
import KratosMultiphysics

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

# Other imports
from importlib import import_module

class LaserDrillingAnalysis(AnalysisStage):
    """
    This class is the main-script of the LaserDrillingApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Making sure that older cases still work by properly initalizing the parameters
        solver_settings = project_parameters["solver_settings"]

        if not solver_settings.Has("domain_size"):
            KratosMultiphysics.Logger.PrintInfo("LaserDrillingAnalysis", "Using the old way to pass the domain_size, this will be removed!")
            solver_settings.AddEmptyValue("domain_size")
            solver_settings["domain_size"].SetInt(project_parameters["problem_data"]["domain_size"].GetInt())

        super(LaserDrillingAnalysis, self).__init__(model, project_parameters)

    def _CreateSolver(self):
        python_module_name = "KratosMultiphysics.LaserDrillingApplication"
        full_module_name = python_module_name + "." + self.project_parameters["solver_settings"]["solver_type"].GetString()
        solver_module = import_module(full_module_name)
        solver = solver_module.CreateSolver(self.model, self.project_parameters["solver_settings"])
        return solver

    def _GetSimulationName(self):
        return "::[Laser-Drilling Simulation]:: "

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 laserdrilling_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 laserdrilling_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = LaserDrillingAnalysis(model, parameters)
    simulation.Run()
