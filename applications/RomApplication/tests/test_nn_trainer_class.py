import numpy as np
import importlib

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess

if kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication"):
    import KratosMultiphysics.StructuralMechanicsApplication
if kratos_utilities.CheckIfApplicationsAvailable("ConstitutiveLawsApplication"):
    import KratosMultiphysics.ConstitutiveLawsApplication

try:
    from KratosMultiphysics.RomApplication.rom_nn_trainer import RomNeuralNetworkTrainer
    have_tensorflow = True
except ImportError:
    have_tensorflow = False

from pathlib import Path

from KratosMultiphysics.RomApplication.rom_manager import RomManager

class SeededNN_RomManager(RomManager):

    def launch_test(self):
        try:
            self._LaunchTrainROM([None])
            self._LaunchFOM([None])
            error = self._LaunchTrainNeuralNetwork([None], [None])
        except Exception as e:
            print(f"An error occurred: {e}")
            return None  # or handle the error as needed
        return error

    def _LaunchTrainNeuralNetwork(self, mu_train, mu_validation):
        try:
            rom_nn_trainer = RomNeuralNetworkTrainer(self.general_rom_manager_parameters, mu_train, mu_validation, self.data_base)
            rom_nn_trainer.TrainNetwork(seed=324)  # Consistent solutions
            self.data_base.add_to_database("Neural_Network", mu_train, None)
            error = rom_nn_trainer.EvaluateNetwork()
        except Exception as e:
            print(f"Error during training or evaluation: {e}")
            return None  # or handle the error as needed
        return error



@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
@KratosUnittest.skipIfApplicationsNotAvailable("ConstitutiveLawsApplication")
@KratosUnittest.skipUnless(have_tensorflow,"Missing required python module: TensorFlow.")
class TestNeuralNetworkTrainerClass(KratosUnittest.TestCase):

    def setUp(self):
        self.reference_relative_error = 0.004697750827444179
        self.work_folder = "nn_trainer_class_test_files"

    def testNeuralNetworkTrainerClass(self):

        general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
                "rom_stages_to_train" : ["ROM"],
                "rom_stages_to_test" : [],
                "type_of_decoder" : "ann_enhanced",
                "paralellism" : null,
                "projection_strategy": "galerkin",
                "assembling_strategy": "global",
                "save_gid_output": false,
                "save_vtk_output": false,
                "output_name": "id",
                "ROM":{
                    "svd_truncation_tolerance": 1e-6,
                    "model_part_name": "Structure",
                    "nodal_unknowns": ["DISPLACEMENT_X","DISPLACEMENT_Y"],
                    "rom_basis_output_format": "numpy",
                    "rom_basis_output_name": "RomParameters",
                    "rom_basis_output_folder": "rom_data",
                    "snapshots_control_type": "step",
                    "snapshots_interval": 1,
                    "print_singular_values": true,
                    "ann_enhanced_settings":{
                        "modes":[1,2],
                        "layers_size":[2,2],
                        "batch_size":2,
                        "epochs":50,
                        "lr_strategy": {
                            "scheduler": "const",
                            "base_lr": 0.01,
                            "additional_params": []
                        }
                    }
                }
            }""")

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            rom_manager = SeededNN_RomManager(general_rom_manager_parameters = general_rom_manager_parameters)
            relative_error = rom_manager.launch_test()
            self.assertAlmostEqual(relative_error, self.reference_relative_error)

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            kratos_utilities.DeleteDirectoryIfExisting(Path('./rom_data'))


##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
