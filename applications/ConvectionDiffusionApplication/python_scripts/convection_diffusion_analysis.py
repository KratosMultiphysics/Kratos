
# Importing Kratos
import KratosMultiphysics
from KratosMultiphysics.ConvectionDiffusionApplication import python_solvers_wrapper_convection_diffusion as solver_wrapper

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

# Other imports
import sys

class ConvectionDiffusionAnalysis(AnalysisStage):
    """
    This class is the main-script of the ConvectionDiffusionApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Making sure that older cases still work by properly initializing the parameters
        solver_settings = project_parameters["solver_settings"]

        if not solver_settings.Has("domain_size"):
            KratosMultiphysics.Logger.PrintInfo("ConvectionDiffusionAnalysis", "Using the old way to pass the domain_size, this will be removed!")
            solver_settings.AddEmptyValue("domain_size")
            solver_settings["domain_size"].SetInt(project_parameters["problem_data"]["domain_size"].GetInt())

        super(ConvectionDiffusionAnalysis, self).__init__(model, project_parameters)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not already in the model) """
        ## Solver construction
        return solver_wrapper.CreateSolverByParameters(self.model, self.project_parameters["solver_settings"],self.project_parameters["problem_data"]["parallel_type"].GetString())

    def _GetSimulationName(self):
        return "::[Convection-Diffusion Simulation]:: "

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 convection_diffusion_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 convection_diffusion_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = ConvectionDiffusionAnalysis(model, parameters)
    simulation.Run()
