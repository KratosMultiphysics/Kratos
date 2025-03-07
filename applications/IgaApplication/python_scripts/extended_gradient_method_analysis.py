# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.IgaApplication 
from KratosMultiphysics.ConvectionDiffusionApplication import python_solvers_wrapper_convection_diffusion as solver_wrapper

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

# Other imports
import sys
import numpy as np
import scipy.sparse.linalg as spla
from shapely.geometry import LineString, Polygon
from scipy.interpolate import Rbf
import matplotlib
matplotlib.use('Agg')  # Use Agg for rendering plots to files
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from colorama import Fore, Style, init

class ExtendedGradientMethodAnalysis(AnalysisStage):
    """
    This class is the main-script of the ExtendedGradientConvectionDiffusion method put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Making sure that older cases still work by properly initalizing the parameters
        self.solver_settings = project_parameters["solver_settings"]
        
        self.model = model

        if not self.solver_settings.Has("domain_size"):
            KratosMultiphysics.Logger.PrintInfo("ExtendedGradientMethodPoissonProblemAnalysis", "Using the old way to pass the domain_size, this will be removed!")
            self.solver_settings.AddEmptyValue("domain_size")
            self.solver_settings["domain_size"].SetInt(project_parameters["problem_data"]["domain_size"].GetInt())

        super(ExtendedGradientMethodAnalysis, self).__init__(model, project_parameters)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        return solver_wrapper.CreateSolverByParameters(self.model, self.project_parameters["solver_settings"],self.project_parameters["problem_data"]["parallel_type"].GetString())

    def _GetSimulationName(self):
        return "::[Extended Gradient Method Poisson Problem IGA Simulation]:: "
    
    def RunSolutionLoop(self):
        # Get the unknown variable for the specific physical problem which is being solved
        unknown_variable_name = self.solver_settings["convection_diffusion_variables"]["unknown_variable"].GetString()
        unknown_variable = KratosMultiphysics.KratosGlobals.GetVariable(unknown_variable_name)

        # Write the output for the initial conditions
        self.OutputSolutionStep()

        while self.KeepAdvancingSolutionLoop():
            # Initialize the variables for the algorithm
            epsilon = 1e-3
            is_not_converged = True
            self.iteration_number = 0
            self.maximum_iterations = 20

            # Initialize the old solution field
            self.old_solution_field = np.zeros(self._GetSolver().GetComputingModelPart().NumberOfNodes())
            for element in self._GetSolver().GetComputingModelPart().Elements:
                for node in element.GetNodes():
                    self.old_solution_field[node.Id - 1] = node.GetSolutionStepValue(unknown_variable) + 1e-2

            self.time = self._AdvanceTime()

            while (is_not_converged == True and self.iteration_number < self.maximum_iterations):
                self._GetSolver().Predict()

                self.InitializeSolutionStep() 

                is_converged = self._GetSolver().SolveSolutionStep()

                solution_error = self.ComputeErrorAndUpdateOldSolutionField()

                print(f"{Fore.BLUE}Iteration # {self.iteration_number}{Style.RESET_ALL}")
                print(f"{Fore.RED}Error between iterations: {solution_error}{Style.RESET_ALL}")

                if (solution_error < epsilon or self.iteration_number == (self.maximum_iterations - 1)):
                    is_not_converged = False
                    print( f"{Fore.GREEN}*** CONVERGENCE ACHIEVED ***{Style.RESET_ALL}")

                self.iteration_number += 1

            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def ComputeErrorAndUpdateOldSolutionField(self):
        # Get the unknown variable for the specific physical problem which is being solved
        unknown_variable_name = self.solver_settings["convection_diffusion_variables"]["unknown_variable"].GetString()
        unknown_variable = KratosMultiphysics.KratosGlobals.GetVariable(unknown_variable_name)

        self.new_solution_field = np.zeros(self._GetSolver().GetComputingModelPart().NumberOfNodes())

        for element in self._GetSolver().GetComputingModelPart().Elements:
            for node in element.GetNodes():
                self.new_solution_field[node.Id - 1] = node.GetSolutionStepValue(unknown_variable)

        error = np.linalg.norm(self.new_solution_field - self.old_solution_field)/np.linalg.norm(self.old_solution_field)
        
        self.old_solution_field = self.new_solution_field

        return error
        

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
    simulation = ExtendedGradientMethodAnalysis(model, parameters)
    simulation.Run()
