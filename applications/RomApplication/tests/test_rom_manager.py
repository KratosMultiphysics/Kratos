import os
import types
import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.RomApplication.rom_manager import RomManager
import json
from pathlib import Path


class TestRomManager(KratosUnittest.TestCase):

    def setUp(self):
        self.work_folder = "test_rom_manager"
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            original_rom_parameters = 'rom_data/RomParameters_test.json'
            with open(original_rom_parameters, 'r') as file:
                data = json.load(file)
            rom_parameters_to_erase = 'rom_data/RomParameters_test_to_erase.json'
            with open(rom_parameters_to_erase, 'w') as file:
                json.dump(data, file, indent=4)



    def test_initialization(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            rom_manager = RomManager()
            self.assertIsNotNone(rom_manager.general_rom_manager_parameters)
            self.assertIsInstance(rom_manager.general_rom_manager_parameters, KratosMultiphysics.Parameters)
            self.assertIsNotNone(rom_manager.project_parameters_name)
            self.assertEqual(rom_manager.project_parameters_name, "ProjectParameters.json")
            self.assertIsNotNone(rom_manager.CustomizeSimulation)
            self.assertTrue(callable(rom_manager.CustomizeSimulation))
            self.assertIsNotNone(rom_manager.UpdateProjectParameters)
            self.assertTrue(callable(rom_manager.UpdateProjectParameters))
            self.assertIsNotNone(rom_manager.UpdateMaterialParametersFile)
            self.assertTrue(callable(rom_manager.UpdateMaterialParametersFile))


    def test_setup_parameters(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            parameters = KratosMultiphysics.Parameters("""{
                "rom_stages_to_train" : ["HROM"],
                "projection_strategy": "lspg",
                "rom_stages_to_test" : ["HROM"],
                "ROM":{
                    "svd_truncation_tolerance": 1e-16 ,
                    "nodal_unknowns": ["VELOCITY_X", "VELOCITY_Y", "PRESSURE"]
                }
            }""")

            rom_manager_object = RomManager(general_rom_manager_parameters=parameters)
            self.assertListEqual(rom_manager_object.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray(), ["HROM"])
            self.assertListEqual(rom_manager_object.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray(), ["HROM"])
            self.assertEqual(rom_manager_object.general_rom_manager_parameters["projection_strategy"].GetString(), "lspg")
            self.assertEqual(rom_manager_object.general_rom_manager_parameters["ROM"]["svd_truncation_tolerance"].GetDouble(), 1e-16)
            self.assertListEqual(rom_manager_object.general_rom_manager_parameters["ROM"]["nodal_unknowns"].GetStringArray(), ["VELOCITY_X", "VELOCITY_Y", "PRESSURE"])


    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def test_orchestration(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):

            #We do not check for specific entries. The generated npy files depend heavily on the SVD, which contains randomization.
            # Here we are checking the generation of the files and their sizes, number of rom modes, number of selected hrom elements, etc.

            def UpdateProjectParameters(parameters, mu=None):
                """
                Customize ProjectParameters here for imposing different conditions to the simulations as needed
                """
                parameters["processes"]["loads_process_list"][0]["Parameters"]["value"].SetString("("+ str(mu[0]) + ")")
                parameters["processes"]["loads_process_list"][1]["Parameters"]["value"].SetString("("+ str(mu[1]) + ")")
                parameters["processes"]["loads_process_list"][2]["Parameters"]["value"].SetString("("+ str(mu[2]) + ")")

                return parameters

            parameters = KratosMultiphysics.Parameters("""{
                "projection_strategy": "galerkin",
                "ROM":{
                    "rom_basis_output_name": "RomParameters_test_to_erase"
                }
                }""")
            rom_manager_object = RomManager(general_rom_manager_parameters=parameters, UpdateProjectParameters = UpdateProjectParameters)
            train_mu = [[0.1, 0.1, 0.1]]
            rom_manager_object.Fit(train_mu)
            right_basis = 'RightBasisMatrix'
            self.assertEqual(np.load(f'rom_data/{right_basis}.npy').shape[1], np.load(f'ExpectedOutputs/galerkin/{right_basis}.npy').shape[1]) #number of right modes
            node_ids = 'NodeIds'
            self.assertListEqual(np.load(f'rom_data/{node_ids}.npy').tolist(),np.load(f'ExpectedOutputs/galerkin/{node_ids}.npy').tolist())
            elem_weights = 'HROM_ElementWeights'
            self.assertEqual(np.load(f'rom_data/{elem_weights}.npy').size, np.load(f'ExpectedOutputs/galerkin/{elem_weights}.npy').size)
            elem_ids = 'HROM_ElementIds'
            self.assertEqual(np.load(f'rom_data/{elem_ids}.npy').size,np.load(f'ExpectedOutputs/galerkin/{elem_ids}.npy').size)
            cond_weights = 'HROM_ConditionWeights'
            self.assertEqual(np.load(f'rom_data/{cond_weights}.npy').size,np.load(f'ExpectedOutputs/galerkin/{cond_weights}.npy').size)
            cond_ids = 'HROM_ConditionIds'
            self.assertEqual(np.load(f'rom_data/{cond_ids}.npy').size,np.load(f'ExpectedOutputs/galerkin/{cond_ids}.npy').size)
            self.clear_npys()

            # TODO: Uncomment this test section after reviewing LSPG HROM (it is not working properly).

            # parameters["projection_strategy"].SetString("lspg")
            # rom_manager_object.Fit(train_mu)
            # self.assertEqual(np.load(f'rom_data/{right_basis}.npy').shape[1], np.load(f'ExpectedOutputs/lspg/{right_basis}.npy').shape[1]) #number of right modes
            # self.assertListEqual(np.load(f'rom_data/{node_ids}.npy').tolist(),np.load(f'ExpectedOutputs/lspg/{node_ids}.npy').tolist())
            # self.assertEqual(np.load(f'rom_data/{elem_weights}.npy').size, np.load(f'ExpectedOutputs/lspg/{elem_weights}.npy').size)
            # self.assertEqual(np.load(f'rom_data/{elem_ids}.npy').size,np.load(f'ExpectedOutputs/lspg/{elem_ids}.npy').size)
            # self.assertEqual(np.load(f'rom_data/{cond_weights}.npy').size,np.load(f'ExpectedOutputs/lspg/{cond_weights}.npy').size)
            # self.assertEqual(np.load(f'rom_data/{cond_ids}.npy').size,np.load(f'ExpectedOutputs/lspg/{cond_ids}.npy').size)
            # self.clear_npys()

            parameters["projection_strategy"].SetString("petrov_galerkin")
            rom_manager_object.Fit(train_mu)
            self.assertEqual(np.load(f'rom_data/{right_basis}.npy').shape[1], np.load(f'ExpectedOutputs/petrov_galerkin/{right_basis}.npy').shape[1]) #number of right modes
            left_basis = 'LeftBasisMatrix'
            self.assertEqual(np.load(f'rom_data/{left_basis}.npy').shape[1], np.load(f'ExpectedOutputs/petrov_galerkin/{left_basis}.npy').shape[1]) #number of left modes
            self.assertListEqual(np.load(f'rom_data/{node_ids}.npy').tolist(),np.load(f'ExpectedOutputs/petrov_galerkin/{node_ids}.npy').tolist())
            self.assertEqual(np.load(f'rom_data/{elem_weights}.npy').size, np.load(f'ExpectedOutputs/petrov_galerkin/{elem_weights}.npy').size)
            self.assertEqual(np.load(f'rom_data/{elem_ids}.npy').size,np.load(f'ExpectedOutputs/petrov_galerkin/{elem_ids}.npy').size)
            self.assertEqual(np.load(f'rom_data/{cond_weights}.npy').size,np.load(f'ExpectedOutputs/petrov_galerkin/{cond_weights}.npy').size)
            self.assertEqual(np.load(f'rom_data/{cond_ids}.npy').size,np.load(f'ExpectedOutputs/petrov_galerkin/{cond_ids}.npy').size)




    def test_error_raising(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            #unavailable projection strategy
            parameters = KratosMultiphysics.Parameters("""{
                        "projection_strategy": "unavailable_strategy"
                    }""")
            rom_manager_object = RomManager(general_rom_manager_parameters=parameters)
            with self.assertRaises(Exception) as context:
                rom_manager_object.Fit()
            err_msg = 'Provided projection strategy unavailable_strategy is not supported. Available options are \'galerkin\', \'lspg\' and \'petrov_galerkin\'.'
            self.assertIn(err_msg, str(context.exception))

            #unavailable flag change
            parameters = KratosMultiphysics.Parameters("""{
                "ROM":{
                    "rom_basis_output_name": "RomParameters_test"
                }
            }""")
            rom_manager_object = RomManager(general_rom_manager_parameters=parameters)
            flag_change = 'unavailable_flag'
            with self.assertRaises(Exception) as context:
                rom_manager_object._ChangeRomFlags(flag_change)
            err_msg = f'Unknown flag "unavailable_flag" change for RomParameters.json'
            self.assertIn(err_msg, str(context.exception))



    def test_flag_changes(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            rom_parameters_to_erase = 'rom_data/RomParameters_test_to_erase.json'
            parameters = KratosMultiphysics.Parameters("""{
                "ROM":{
                    "rom_basis_output_name": "RomParameters_test_to_erase"
                }
            }""")
            rom_manager_object = RomManager(general_rom_manager_parameters=parameters)
            flag_change = 'GalerkinROM'
            rom_manager_object._ChangeRomFlags(flag_change)
            with open(rom_parameters_to_erase, 'r') as file:
                data = json.load(file)
                self.assertEqual(data['projection_strategy'], 'galerkin')
                self.assertEqual(data['train_hrom'], False)
                self.assertEqual(data['run_hrom'], False)
                # what to do with self._SetGalerkinBnSParameters()

            flag_change = 'trainHROMGalerkin'
            rom_manager_object._ChangeRomFlags(flag_change)
            with open(rom_parameters_to_erase, 'r') as file:
                data = json.load(file)
                self.assertEqual(data['projection_strategy'], 'galerkin')
                self.assertEqual(data['train_hrom'], True)
                self.assertEqual(data['run_hrom'], False)
                # what to do with self._SetGalerkinBnSParameters()

            flag_change = 'runHROMGalerkin'
            rom_manager_object._ChangeRomFlags(flag_change)
            with open(rom_parameters_to_erase, 'r') as file:
                data = json.load(file)
                self.assertEqual(data['projection_strategy'], 'galerkin')
                self.assertEqual(data['train_hrom'], False)
                self.assertEqual(data['run_hrom'], True)
                # what to do with self._SetGalerkinBnSParameters()

            flag_change = 'lspg'
            rom_manager_object._ChangeRomFlags(flag_change)
            with open(rom_parameters_to_erase, 'r') as file:
                data = json.load(file)
                self.assertEqual(data['projection_strategy'], 'lspg')
                self.assertEqual(data['train_hrom'], False)
                self.assertEqual(data['run_hrom'], False)
                # what to do with self._SetLSPGBnSParameters()

            flag_change = 'trainHROMLSPG'
            rom_manager_object._ChangeRomFlags(flag_change)
            with open(rom_parameters_to_erase, 'r') as file:
                data = json.load(file)
                self.assertEqual(data['projection_strategy'], 'lspg')
                self.assertEqual(data['train_hrom'], True)
                self.assertEqual(data['run_hrom'], False)
                # what to do with self._SetLSPGBnSParameters()

            flag_change = 'runHROMLSPG'
            rom_manager_object._ChangeRomFlags(flag_change)
            with open(rom_parameters_to_erase, 'r') as file:
                data = json.load(file)
                self.assertEqual(data['projection_strategy'], 'lspg')
                self.assertEqual(data['train_hrom'], False)
                self.assertEqual(data['run_hrom'], True)
                # what to do with self._SetLSPGBnSParameters()

            flag_change = 'TrainPG'
            rom_manager_object._ChangeRomFlags(flag_change)
            with open(rom_parameters_to_erase, 'r') as file:
                data = json.load(file)
                self.assertEqual(data['projection_strategy'], 'lspg')
                self.assertEqual(data['train_hrom'], False)
                self.assertEqual(data['run_hrom'], False)
                self.assertEqual(data["rom_settings"]['rom_bns_settings']['train_petrov_galerkin'], True)
                # what to do with self._SetLSPGBnSParameters()

            flag_change = 'PG'
            rom_manager_object._ChangeRomFlags(flag_change)
            with open(rom_parameters_to_erase, 'r') as file:
                data = json.load(file)
                self.assertEqual(data['projection_strategy'], 'petrov_galerkin')
                self.assertEqual(data['train_hrom'], False)
                self.assertEqual(data['run_hrom'], False)
                # what to do with self._SetPetrovGalerkinBnSParameters()

            flag_change = 'trainHROMPetrovGalerkin'
            rom_manager_object._ChangeRomFlags(flag_change)
            with open(rom_parameters_to_erase, 'r') as file:
                data = json.load(file)
                self.assertEqual(data['projection_strategy'], 'petrov_galerkin')
                self.assertEqual(data['train_hrom'], True)
                self.assertEqual(data['run_hrom'], False)
                # what to do with self._SetPetrovGalerkinBnSParameters()

            flag_change = 'runHROMPetrovGalerkin'
            rom_manager_object._ChangeRomFlags(flag_change)
            with open(rom_parameters_to_erase, 'r') as file:
                data = json.load(file)
                self.assertEqual(data['projection_strategy'], 'petrov_galerkin')
                self.assertEqual(data['train_hrom'], False)
                self.assertEqual(data['run_hrom'], True)
                # what to do with self._SetPetrovGalerkinBnSParameters()


    def clear_npys(self):
        #this should be called from the rom_data scope
        for file_name in os.listdir():
                if file_name.endswith(".npy"):
                    kratos_utilities.DeleteFileIfExisting(file_name)



    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            for file_name in os.listdir():
                if file_name.endswith(".time") or file_name.endswith("HROM.mdpa") or file_name.endswith("HROMVisualization.mdpa"):
                    kratos_utilities.DeleteFileIfExisting(file_name)
        with KratosUnittest.WorkFolderScope(self.work_folder+'/rom_data', __file__):
            kratos_utilities.DeleteDirectoryIfExisting(Path('./rom_database/'))
            for file_name in os.listdir():
                if file_name.endswith("test_to_erase.json"):
                    kratos_utilities.DeleteFileIfExisting(file_name)
            self.clear_npys()


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
