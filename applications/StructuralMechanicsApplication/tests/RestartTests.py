import os

# Import Kratos
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import Kratos_Execute_Structural_Test as Execute_Test

import KratosMultiphysics.kratos_utilities as kratos_utils

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

            # Set common settings
            self.time_step = 1.5
            self.start_time = 0.0
            self.end_time = 5.9

            self.project_parameters_save["problem_data"]["time_step"].SetDouble(self.time_step)
            self.project_parameters_save["problem_data"]["start_time"].SetDouble(self.start_time)
            self.project_parameters_save["problem_data"]["end_time"].SetDouble(self.end_time)

            self.project_parameters_save["solver_settings"].RemoveValue("restart_settings") # to start clean

            # Now clone the settings after the common settings are set
            self.project_parameters_load = self.project_parameters_save.Clone()

            # Adding the specific settings (minimal to test the default settings)
            save_restart_parameters = KratosMultiphysics.Parameters("""{
                "save_restart" : true
            }""")
            load_restart_parameters = KratosMultiphysics.Parameters("""{
                "load_restart"              : true,
                "restart_load_file_label"   : "3.0"
            }""")

            self.project_parameters_save["solver_settings"].AddValue("restart_settings", save_restart_parameters)
            self.project_parameters_load["solver_settings"].AddValue("restart_settings", load_restart_parameters)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            test = Execute_Test.Kratos_Execute_Test(self.project_parameters_save)
            test.Solve()

            test = Execute_Test.Kratos_Execute_Test(self.project_parameters_load)
            test.Solve()

    def tearDown(self):
        # remove the created restart files
        raw_path, raw_file_name = os.path.split(self.file_name)
        folder_name = os.path.join(raw_path, raw_file_name + "__restart_files")

        t = self.start_time
        while(t < self.end_time): # using while bcs range does not support floats
            t +=  self.time_step
            file_path = os.path.join(folder_name, raw_file_name + "_" + str(t) + ".rest")
            kratos_utils.DeleteFileIfExisting(file_path)

        if os.path.isdir(folder_name):
            os.rmdir(folder_name)


class TestSmallDisplacement2D4N(StructuralMechanicsRestartTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_shear_qua"

class TestTotalLagrangian2D3N(StructuralMechanicsRestartTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_shear_tri"

class TestUpdatedLagrangian3D8N(StructuralMechanicsRestartTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_shear_hexa"
