import os
import json
import shutil
import multiprocessing
import KratosMultiphysics
import KratosMultiphysics.ThermalDEMApplication as ThermalDEM
import KratosMultiphysics.ThermalDEMApplication.thermal_dem_analysis
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import Logger

result_print = True

#---------------------------------------------------------------------------------
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

#---------------------------------------------------------------------------------
class Solution(ThermalDEM.thermal_dem_analysis.ThermalDEMAnalysis, KratosUnittest.TestCase):
    def __init__(self, model, project_parameters, test_folder_name, solution_name, solution_file_path):
        self.test_folder_name = test_folder_name
        self.solution_name = solution_name
        ReadReferenceSolutionFile(self, solution_file_path)
        super().__init__(model, project_parameters)

    def GetMainPath(cls):
        return GetTestFolderPath(cls.test_folder_name)

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters['problem_name'].GetString())

    def Finalize(self):
        CheckResults(self)
        CleanFolders(self)
        super().Finalize()

#---------------------------------------------------------------------------------
def GetTestFolderPath(test_folder_name):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), test_folder_name)

#---------------------------------------------------------------------------------
def CreateAndRunStage(model, test_folder_name, parameters_file_name, solution_file_name, solution_name, number_of_threads):
    # Get file paths
    parameters_file_path = os.path.join(GetTestFolderPath(test_folder_name), parameters_file_name)
    solution_file_path   = os.path.join(GetTestFolderPath(test_folder_name), solution_file_name)

    # Set number of threads
    if "OMP_NUM_THREADS" in os.environ:
        initial_number_of_threads = os.environ['OMP_NUM_THREADS']
        KratosMultiphysics.ParallelUtilities.SetNumThreads(number_of_threads)
    else:
        Logger.PrintInfo("cores detected:", multiprocessing.cpu_count())
        KratosMultiphysics.ParallelUtilities.SetNumThreads(number_of_threads)

    # Read project parameters
    with open(parameters_file_path, 'r') as parameter_file:
        project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

    # Run test
    Solution(model, project_parameters, test_folder_name, solution_name, solution_file_path).Run()

    # Restore number of threads
    if "OMP_NUM_THREADS" in os.environ:
        KratosMultiphysics.ParallelUtilities.SetNumThreads(int(initial_number_of_threads))
    else:
        KratosMultiphysics.ParallelUtilities.SetNumThreads(multiprocessing.cpu_count())

#---------------------------------------------------------------------------------
def ReadReferenceSolutionFile(my_solution, solution_file_path):
    my_solution.result_types = None
    my_solution.result_nodes = None
    my_solution.result_tol   = None
    my_solution.result_refs  = None

    try:
        with open(solution_file_path, "r") as file:
            data = json.load(file)

        my_solution.result_types = data.get("result_types")
        my_solution.result_nodes = data.get("result_nodes")
        my_solution.result_tol   = data.get("tolerance")
        my_solution.result_refs  = data.get("solutions")

        if any(attr is None for attr in [my_solution.result_types, my_solution.result_nodes, my_solution.result_tol, my_solution.result_refs]):
            raise ValueError("One or more required fields are missing in the reference solution file.")

    except FileNotFoundError:
        print(f"Error: Reference solution file not found.")
    except json.JSONDecodeError:
        print(f"Error: Reference solution file is not a valid JSON file.")
    except Exception as e:
        print(f"An error occurred: {e}")

#---------------------------------------------------------------------------------
def CheckResults(my_solution):
    if my_solution.solution_name not in my_solution.result_refs:
        raise ValueError(f"Unknown solution name: {my_solution.solution_name}")

    refs = my_solution.result_refs[my_solution.solution_name]
    computed_results = {}

    if result_print:
        print(f"\nComputed results for {my_solution.solution_name}:")
        print("{")

    for node_id in my_solution.result_nodes:
        node = my_solution.spheres_model_part.GetNode(node_id)
        computed_values = []
        for result_type in my_solution.result_types:
            result_value = node.GetSolutionStepValue(eval(result_type))
            computed_values.append(result_value)
        computed_results[node_id] = computed_values
        if result_print:
            print(f"    {node_id}: {computed_values},")
    if result_print:
        print("}\n")

    for node_id, computed_values in computed_results.items():
        for i, result_value in enumerate(computed_values):
            my_solution.assertAlmostEqual(result_value, refs[str(node_id)][i], delta=my_solution.result_tol)

#---------------------------------------------------------------------------------
def CleanFolders(my_solution):
    my_solution.procedures.RemoveFoldersWithResults(my_solution.main_path, my_solution.problem_name, '')

    pycache_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '__pycache__')
    if os.path.exists(pycache_path) and os.path.isdir(pycache_path):
        shutil.rmtree(pycache_path)
