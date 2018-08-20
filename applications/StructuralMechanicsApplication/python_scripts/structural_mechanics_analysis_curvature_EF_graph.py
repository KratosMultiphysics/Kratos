from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Importing the base class
from structural_mechanics_analysis import StructuralMechanicsAnalysis

class StructuralMechanicsAnalysisCurvatureEFGraph(StructuralMechanicsAnalysis):
    """
    This class is the main-script of the StructuralMechanicsApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        #solver_settings = project_parameters["solver_settings"]
        self.time_array = [0.15,0.2,0.25]

        super(StructuralMechanicsAnalysisCurvatureEFGraph, self).__init__(model, project_parameters)

    def RunSolutionLoop(self):
        while self.time < self.end_time:
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 structural_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 structural_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysis(model, parameters)
    simulation.Run()
