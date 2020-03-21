import os

# Import Kratos
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_analysis

import KratosMultiphysics.kratos_utilities as kratos_utils

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

# This utility will control the execution scope in case we need to access files or we depend
# on specific relative locations of the files.

# TODO: Should we move this to KratosUnittest?
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, the_type, value, traceback):
        os.chdir(self.currentPath)

class StructuralMechanicsRestartTestFactory(KratosUnittest.TestCase):
    """
    This class first runs a simulation and saves a restart file
    In the second step the created restart file is read and the simulation
    is started from it.

    The results are being compared in both cases

    The restart-related setup of the test is automatic as well as the comparison
    is automatic, all the user has to specify is the general settings, mdpa-
    and the materials-file
    """
    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # here we read the general parameters and add the load/save specific settings
            with open(self.file_name + "_parameters.json", 'r') as parameter_file:
                self.project_parameters_save = KratosMultiphysics.Parameters(parameter_file.read())

            # To avoid many prints
            if (self.project_parameters_save["problem_data"]["echo_level"].GetInt() == 0):
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

            # Set common settings
            self.time_step = 1.5
            self.start_time = 0.0
            self.end_time = 5.9

            self.project_parameters_save["solver_settings"]["time_stepping"]["time_step"].SetDouble(self.time_step)
            self.project_parameters_save["problem_data"]["start_time"].SetDouble(self.start_time)
            self.project_parameters_save["problem_data"]["end_time"].SetDouble(self.end_time)

            # Now clone the settings after the common settings are set
            self.project_parameters_load = self.project_parameters_save.Clone()

            # Adding the specific settings (minimal to test the default settings)
            save_restart_parameters = KratosMultiphysics.Parameters("""{
                "restart_processes" : [
                    {
                    "python_module"   : "save_restart_process",
                    "kratos_module"   : "KratosMultiphysics",
                    "process_name"    : "SaveRestartProcess",
                    "Parameters"            : {
                    }
                }]
             }""")

            save_restart_parameters["restart_processes"][0]["Parameters"].AddValue("model_part_name", self.project_parameters_save["solver_settings"]["model_part_name"])

            self.project_parameters_save.AddValue("output_processes", save_restart_parameters)

            load_mp_import_settings = self.project_parameters_load["solver_settings"]["model_import_settings"]
            load_mp_import_settings.AddEmptyValue("restart_load_file_label")
            load_mp_import_settings["input_type"].SetString("rest")
            load_mp_import_settings["restart_load_file_label"].SetString("3.0")

            # Correct the path
            restart_file_path = load_mp_import_settings["input_filename"].GetString()
            restart_file_path = load_mp_import_settings["input_filename"].SetString("Structure")

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model_save = KratosMultiphysics.Model()
            model_load = KratosMultiphysics.Model()
            structural_mechanics_analysis.StructuralMechanicsAnalysis(model_save, self.project_parameters_save).Run()
            structural_mechanics_analysis.StructuralMechanicsAnalysis(model_load, self.project_parameters_load).Run()

    def tearDown(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # remove the created restart files
            #raw_path, raw_file_name = os.path.split(self.file_name)
            raw_file_name = "Structure"
            folder_name = raw_file_name + "__restart_files"

            kratos_utils.DeleteDirectoryIfExisting(GetFilePath(folder_name))


class TestSmallDisplacement2D4N(StructuralMechanicsRestartTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_shear_qua"

class TestTotalLagrangian2D3N(StructuralMechanicsRestartTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_shear_tri"

class TestUpdatedLagrangian3D8N(StructuralMechanicsRestartTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_shear_hexa"


if __name__ == "__main__":
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSmallDisplacement2D4N]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTotalLagrangian2D3N]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestUpdatedLagrangian3D8N]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)
