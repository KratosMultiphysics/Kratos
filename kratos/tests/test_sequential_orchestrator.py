import pathlib
import importlib

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest
import KratosMultiphysics.json_utilities
import KratosMultiphysics.kratos_utilities
from KratosMultiphysics.project import Project
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.analysis_stage import AnalysisStage

# Fake do-nothing solver class to test the orchestrator
class EmptySolver(PythonSolver):
    def __init__(self, model, settings):
        super().__init__(model, settings)

    @classmethod
    def GetDefaultParameters(cls):
        return KratosMultiphysics.Parameters("""{
            "echo_level" : 0,
            "model_part_name" : ""
        }""")

    def ImportModelPart(self):
        pass

    def GetComputingModelPart(self):
        return self.model.GetModelPart(self.settings["model_part_name"].GetString())

# Fake do-nothing analysis stage class to test the orchestrator
class EmptyAnalysisStage(AnalysisStage):
    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)

    def _CreateSolver(self):
        return EmptySolver(self.model, self.project_parameters["solver_settings"])

    def KeepAdvancingSolutionLoop(self):
        return False

# Fake do-nothing operation to test the orchestrator
class EmptyOperation(KratosMultiphysics.Operation):
    def __init__(self, model, project_parameters):
        super().__init__()

    def Execute(self):
        pass

# Register the fake classes to build the multistage simulation
# These will be removed from the registry when executing the test tearDown method
# Note that we check if the items have been already added in case this module is imported again by the factories creating the classes above
if not KratosMultiphysics.Registry.HasItem("Stages.KratosMultiphysics.EmptyAnalysisStage"):
    KratosMultiphysics.Registry.AddItem("Stages.KratosMultiphysics.EmptyAnalysisStage.ClassName", "EmptyAnalysisStage")
if not KratosMultiphysics.Registry.HasItem("Operations.KratosMultiphysics.EmptyOperation"):
    KratosMultiphysics.Registry.AddItem("Operations.KratosMultiphysics.EmptyOperation.ModuleName", "test_sequential_orchestrator")

class TestSequentialOrchestrator(KratosMultiphysics.KratosUnittest.TestCase):

    def test_sequential_orchestrator(self):
        # Parse sequential orchestrator settings
        project_parameters_filename = pathlib.Path("test_files/orchestrators_files/test_sequential_orchestrator.json")
        with open(project_parameters_filename, 'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

        # Construct the project instance with previous settings
        project = Project(project_parameters)

        # Instantiate the orchestrator and execute it
        orchestrator_reg_entry = KratosMultiphysics.Registry[project.GetSettings()["orchestrator"]["name"].GetString()]
        orchestrator_module = importlib.import_module(orchestrator_reg_entry["ModuleName"])
        orchestrator_class = getattr(orchestrator_module, orchestrator_reg_entry["ClassName"])
        orchestrator_instance = orchestrator_class(project)
        orchestrator_instance.Run()

        # Check that the checkpoint files exist
        # Note that content cannot be checked as some values depend on the compilation
        self.assertTrue(pathlib.Path.exists(pathlib.Path("test_files/orchestrators_files/checkpoints/stage_1")))
        self.assertTrue(pathlib.Path.exists(pathlib.Path("test_files/orchestrators_files/checkpoints/stage_2")))
        self.assertTrue(pathlib.Path.exists(pathlib.Path("test_files/orchestrators_files/checkpoints/test_sequential_orchestrator_stage_1.json")))
        self.assertTrue(pathlib.Path.exists(pathlib.Path("test_files/orchestrators_files/checkpoints/test_sequential_orchestrator_stage_2.json")))

    def tearDown(self):
        # Remove the checkpointing files
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting("test_files/orchestrators_files/checkpoints")

        # Remove the auxiliary EmptyAnalysisStage from the registry
        KratosMultiphysics.Registry.RemoveItem("Stages.KratosMultiphysics.EmptyAnalysisStage")
        KratosMultiphysics.Registry.RemoveItem("Operations.KratosMultiphysics.EmptyOperation")

if __name__ == '__main__':
    KratosMultiphysics.KratosUnittest.main()
