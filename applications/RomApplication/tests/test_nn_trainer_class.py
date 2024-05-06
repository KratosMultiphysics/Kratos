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


@KratosUnittest.skipUnless(have_tensorflow,"Missing required python module: TensorFlow.")
class TestNeuralNetworkTrainerClass(KratosUnittest.TestCase):

    def setUp(self):
        self.reference_relative_error = 0.0009752321235806322
        self.work_folder = "nn_trainer_class_test_files"

    def testNeuralNetworkTrainerClass(self):

        nn_training_parameters = KratosMultiphysics.Parameters("""{
            "ROM": {
                "ann_enhanced_settings":{
                    "saved_models_root_path": "rom_data/saved_nn_models/",
                    "training":{
                        "modes":[1,2],
                        "layers_size":[2,2],
                        "batch_size":2,
                        "epochs":50,
                        "lr_strategy": {
                            "scheduler": "const",         // "const", "steps", "sgdr"
                            "base_lr": 0.01,
                            "additional_params": []     // const:[], steps/sgdr:["min_lr", "reduction_factor","update_period"]
                        },
                        "database":{
                            "training_set": "rom_data/reference_fom_snapshots.npy",
                            "validation_set": "rom_data/reference_fom_snapshots_val.npy",
                            "phi_matrix": "rom_data/reference_RightBasisMatrix.npy",
                            "sigma_vector": "rom_data/reference_SingularValuesVector.npy"
                        },
                        "use_automatic_name": false,
                        "custom_name": "test_neural_network"
                    },
                    "online":{
                        "model_name": "test_neural_network"
                    }
                }
            }
        }""")

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):

            rom_neural_network_trainer = RomNeuralNetworkTrainer(nn_training_parameters)
            model_name = rom_neural_network_trainer.TrainNetwork(seed = 324)
            relative_error = rom_neural_network_trainer.EvaluateNetwork(model_name)

            self.assertAlmostEqual(relative_error, self.reference_relative_error)

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    @KratosUnittest.skipIfApplicationsNotAvailable("ConstitutiveLawsApplication")
    def testGenerateTrainingDatabases(self):

        parameters_filename = "ProjectParameters.json"

        rom_training_parameters = KratosMultiphysics.Parameters("""{
            "python_module" : "calculate_rom_basis_output_process",
            "kratos_module" : "KratosMultiphysics.RomApplication",
            "process_name"  : "CalculateRomBasisOutputProcess",
            "help"          : "This process should write the Rom basis",
            "Parameters"    :{
                "model_part_name": "Structure",
                "rom_manager" : false,
                "snapshots_control_type": "step",
                "snapshots_interval": 1.0,
                "nodal_unknowns":  ["DISPLACEMENT_X","DISPLACEMENT_Y"],
                "rom_basis_output_format": "numpy",
                "rom_basis_output_name": "RomParameters",
                "rom_basis_output_folder": "rom_data",
                "svd_truncation_tolerance": 1e-6,
                "print_singular_values": true
            }
        }""")

        # List containing the fake total times
        test_timesteps = [4,8,13,16]

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):

            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            parameters["output_processes"].AddEmptyArray("rom_output")
            parameters["output_processes"]["rom_output"].Append(rom_training_parameters)

            analysis_stage_module_name = parameters["analysis_stage"].GetString()
            analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
            analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))
            analysis_stage_module = importlib.import_module(analysis_stage_module_name)
            analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

            model = KratosMultiphysics.Model()
            simulation = analysis_stage_class(model,parameters)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    basis_output_process = process
            snapshots_matrix_train=basis_output_process._GetSnapshotsMatrix()

            snapshots_matrix_test = snapshots_matrix_train[:,test_timesteps].copy()
            snapshots_matrix_train = np.delete(snapshots_matrix_train, test_timesteps, 1)

            u, sigma = basis_output_process._ComputeSVD(snapshots_matrix_train)
            basis_output_process._PrintRomBasis(u, sigma)

            self.assertTrue(np.all(np.abs(snapshots_matrix_train - np.load('rom_data/reference_fom_snapshots.npy')<1e-10)))
            self.assertTrue(np.all(np.abs(snapshots_matrix_test - np.load('rom_data/reference_fom_snapshots_val.npy'))<1e-10))
            self.assertTrue(np.all(np.abs(np.load('rom_data/SingularValuesVector.npy') - np.load('rom_data/reference_SingularValuesVector.npy'))<1e-10))
            self.assertTrue(np.all(np.abs(np.load('rom_data/RightBasisMatrix.npy') - np.load('rom_data/reference_RightBasisMatrix.npy'))<1e-10))

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            rom_data_path = 'rom_data/'
            neural_network_models_path = rom_data_path+'saved_nn_models/test_neural_network/'
            kratos_utilities.DeleteDirectoryIfExisting(neural_network_models_path)
            kratos_utilities.DeleteFileIfExisting(rom_data_path+'RightBasisMatrix.npy')
            kratos_utilities.DeleteFileIfExisting(rom_data_path+'SingularValuesVector.npy')
            kratos_utilities.DeleteFileIfExisting(rom_data_path+'NodeIds.npy')
            kratos_utilities.DeleteFileIfExisting(rom_data_path+'RomParameters.json')

##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
