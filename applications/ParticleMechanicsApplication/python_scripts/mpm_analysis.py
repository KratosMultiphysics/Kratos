from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos Core, Applications and Dependencies
import KratosMultiphysics

# Importing the solvers (if available)
try:
    import KratosMultiphysics.ExternalSolversApplication
    KratosMultiphysics.Logger.PrintInfo("ExternalSolversApplication", "succesfully imported")
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("ExternalSolversApplication", "not imported")

try:
    import KratosMultiphysics.StructuralMechanicsApplication
    KratosMultiphysics.Logger.PrintInfo("StructuralMechanicsApplication", "succesfully imported")
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("StructuralMechanicsApplication", "not imported")


# Importing the base class
from analysis_stage import AnalysisStage
import KratosMultiphysics.ParticleMechanicsApplication.python_solvers_wrapper_mpm as solvers_wrapper

class MpmAnalysis(AnalysisStage):
    """
    This class is the analysis script for MPM of the ParticleMechanicsApplication
    """

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        print('Create Solver')
        return solvers_wrapper.CreateSolver(self.model, self.project_parameters)

    def _GetSimulationName(self):
        return "::[Particle Mechanics Analysis]:: "

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 particle_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 particle_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = MpmAnalysis(model, parameters)
    simulation.Run()
