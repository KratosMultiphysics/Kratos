# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.IgaApplication 
from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_structural as structural_solvers

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

# Other imports
import sys
import numpy as np
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')  # Use Agg for rendering plots to files
import matplotlib.pyplot as plt
from colorama import Fore, Style, init

# Ensure full printing of matrices
np.set_printoptions(threshold=np.inf, linewidth=200, suppress=True)

class ExtendedGradientMethodStructuralMechanicsAnalysis(AnalysisStage):
    """
    This class is the main-script of the ExtendedGradientMethodStructuralMechanicsAnalysis method put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Making sure that older cases still work by properly initalizing the parameters
        self.solver_settings = project_parameters["solver_settings"]
        
        self.model = model

        if not self.solver_settings.Has("domain_size"):
            KratosMultiphysics.Logger.PrintInfo("ExtendedGradientMethodStructuralMechanicsAnalysis", "Using the old way to pass the domain_size, this will be removed!")
            self.solver_settings.AddEmptyValue("domain_size")
            self.solver_settings["domain_size"].SetInt(project_parameters["problem_data"]["domain_size"].GetInt())

        super(ExtendedGradientMethodStructuralMechanicsAnalysis, self).__init__(model, project_parameters)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        return structural_solvers.CreateSolverByParameters(self.model, self.project_parameters["solver_settings"],self.project_parameters["problem_data"]["parallel_type"].GetString())

    def _GetSimulationName(self):
        return "::[Extended Gradient Method Structural Mechanics IGA Simulation]:: "
    
    def RunSolutionLoop(self):
        # Get the unknown variables for the displacement components
        displacement_variables = [
            KratosMultiphysics.KratosGlobals.GetVariable("DISPLACEMENT_X"),
            KratosMultiphysics.KratosGlobals.GetVariable("DISPLACEMENT_Y"),
            KratosMultiphysics.KratosGlobals.GetVariable("DISPLACEMENT_Z")
        ]

        # Write the output for the initial conditions
        self.OutputSolutionStep()

        while self.KeepAdvancingSolutionLoop():
            # Initialize the variables for the algorithm
            epsilon = 1e-2
            is_not_converged = True
            self.iteration_number = 0
            self.maximum_iterations = 3

            # Initialize the old solution field
            num_nodes = self._GetSolver().GetComputingModelPart().NumberOfNodes()
            self.old_solution_field = np.zeros((num_nodes, 3))  # Store all three components

            for element in self._GetSolver().GetComputingModelPart().Elements:
                for node in element.GetNodes():
                    node_id = node.Id - 1
                    for i, var in enumerate(displacement_variables):
                        self.old_solution_field[node_id, i] = node.GetSolutionStepValue(var) 
            
            self.time = self._AdvanceTime()

            while is_not_converged and self.iteration_number < self.maximum_iterations:
                self._GetSolver().Predict()
                self.InitializeSolutionStep()
                    
                is_converged = self._GetSolver().SolveSolutionStep()
                
                solution_error = self.ComputeErrorAndUpdateOldSolutionField()
                
                print(f"{Fore.BLUE}Iteration # {self.iteration_number}{Style.RESET_ALL}")
                print(f"{Fore.RED}Error between iterations: {solution_error}{Style.RESET_ALL}")
                
                if solution_error < epsilon or self.iteration_number == (self.maximum_iterations - 1):
                    is_not_converged = False
                    print(f"{Fore.GREEN}*** CONVERGENCE ACHIEVED ***{Style.RESET_ALL}")

                self.iteration_number += 1

            self.FinalizeSolutionStep()
            self.OutputSolutionStep()


    def ComputeErrorAndUpdateOldSolutionField(self):
        # Get the unknown variables for the displacement components
        displacement_variables = [
            KratosMultiphysics.KratosGlobals.GetVariable("DISPLACEMENT_X"),
            KratosMultiphysics.KratosGlobals.GetVariable("DISPLACEMENT_Y"),
            KratosMultiphysics.KratosGlobals.GetVariable("DISPLACEMENT_Z")
        ]
        
        num_nodes = self._GetSolver().GetComputingModelPart().NumberOfNodes()
        self.new_solution_field = np.zeros((num_nodes, 3))  # Store all three components
        
        for element in self._GetSolver().GetComputingModelPart().Elements:
            for node in element.GetNodes():
                node_id = node.Id - 1
                for i, var in enumerate(displacement_variables):
                    self.new_solution_field[node_id, i] = node.GetSolutionStepValue(var)

        # Compute norms of the old solution field per component
        old_norms = np.linalg.norm(self.old_solution_field, axis=0)

        # Identify valid components (where old_norm is not zero)
        valid_components = old_norms != 0  

        # Compute error only for valid components
        error_vector = np.zeros(3)  # Initialize error vector with zeros
        # error_vector[valid_components] = (
        #     np.linalg.norm(self.new_solution_field - self.old_solution_field, axis=0)[valid_components] / old_norms[valid_components]
        # )
        error_vector[valid_components] = (
            np.linalg.norm(self.new_solution_field - self.old_solution_field, axis=0)[valid_components])

        print(error_vector)

        # Compute the average error only over valid components
        if np.any(valid_components):  # Check if there is at least one valid component
            error = np.mean(error_vector[valid_components])
        else:
            error = 0.0  # If all components are zero, return zero error

        # Update the old solution field
        self.old_solution_field = self.new_solution_field.copy()

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
    simulation = ExtendedGradientMethodPoissonProblemAnalysis(model, parameters)
    simulation.Run()
