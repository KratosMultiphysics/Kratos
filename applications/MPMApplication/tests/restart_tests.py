import os
import pathlib

# Import Kratos
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.MPMApplication import particle_mechanics_analysis

import KratosMultiphysics.kratos_utilities as kratos_utils

def GetFilePath(fileName):
    return str(pathlib.Path(__file__).absolute().parent / fileName)
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, the_type, value, traceback):
        os.chdir(self.currentPath)

# This utility will control the execution scope in case we need to access files or we depend
# on specific relative locations of the files.

class MPMRestartTestFactory(KratosUnittest.TestCase):
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
        with KratosUnittest.WorkFolderScope("", __file__):
            # here we read the general parameters and add the load/save specific settings
            with open(self.file_name + "_parameters.json", 'r') as parameter_file:
                self.project_parameters_save = KratosMultiphysics.Parameters(parameter_file.read())

            # To avoid many prints
            if (self.project_parameters_save["problem_data"]["echo_level"].GetInt() == 0):
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

            self.time_step = self.project_parameters_save["solver_settings"]["time_stepping"]["time_step"].GetDouble()
            self.start_time = self.project_parameters_save["problem_data"]["start_time"].GetDouble()
            self.end_time = self.project_parameters_save["problem_data"]["end_time"].GetDouble()

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

            model_part_name = self.project_parameters_save["solver_settings"]["model_part_name"]
            save_restart_parameters["restart_processes"][0]["Parameters"].AddValue("model_part_name", model_part_name)

            self.project_parameters_save.AddValue("output_processes", save_restart_parameters)

            load_mp_import_settings = self.project_parameters_load["solver_settings"]["model_import_settings"]
            load_mp_import_settings.AddEmptyValue("restart_load_file_label")
            load_mp_import_settings["input_type"].SetString("rest")
            load_mp_import_settings["restart_load_file_label"].SetString(str(self.start_time+self.time_step))

            # Correct the path
            load_mp_import_settings["input_filename"].SetString("MPM_Material")

    def test_execution(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope("", __file__):
            # remove the created restart files
            raw_file_name = "MPM_Material"
            folder_name = raw_file_name + "__restart_files"
            self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, GetFilePath(folder_name))

            model_save = KratosMultiphysics.Model()
            model_load = KratosMultiphysics.Model()
            particle_mechanics_analysis.ParticleMechanicsAnalysis(model_save, self.project_parameters_save).Run()
            particle_mechanics_analysis.ParticleMechanicsAnalysis(model_load, self.project_parameters_load).Run()


class MPMRestartTestDynamicCantilever2D(MPMRestartTestFactory):
    file_name = "beam_tests/dynamic_cantilever/dynamic_cantilever_consistent_mass_test"

class MPMRestartTestBeamStaticLineLoad2D(MPMRestartTestFactory):
    file_name = "beam_tests/cantilever_beam/static_line_load_2D_quad_test"


if __name__ == "__main__":
    suites = KratosUnittest.KratosSuites
    small_suite = suites['small'] # These tests are executed by the continuous integration tool
    small_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MPMRestartTestBeamStaticLineLoad2D]))
    small_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MPMRestartTestDynamicCantilever2D]))
    all_suite = suites['all']
    all_suite.addTests(small_suite)
    KratosUnittest.runTests(suites)
