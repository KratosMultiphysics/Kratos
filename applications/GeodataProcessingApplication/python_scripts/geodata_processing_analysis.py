from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.GeodataProcessingApplication import geodata_processing_solver as geodata_solver

class GeodataProcessingAnalysis(AnalysisStage):
    """
    This class is the main-script of the GeodataProcessingApplication put in a class
    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        super(GeodataProcessingAnalysis, self).__init__(model, project_parameters)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not already in the model) """
        ## Solver construction
        return geodata_solver.CreateSolver(self.model, self.project_parameters)


if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 geodata_processing_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 geodata_processing_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = GeodataProcessingAnalysis(model, parameters)
    simulation.Run()