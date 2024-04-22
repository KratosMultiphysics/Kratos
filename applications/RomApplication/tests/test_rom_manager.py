import os
import types
import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
#from unittest.mock import patch, MagicMock
from KratosMultiphysics.RomApplication.rom_manager import RomManager







class TestRomManager(KratosUnittest.TestCase):
    def setUp(self):
        self.rom_manager = RomManager()

    def test_initialization(self):
        self.assertIsNotNone(self.rom_manager.general_rom_manager_parameters)

    def test_setup_parameters(self):
        parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train" : ["ROM","HROM"],
            "projection_strategy": "lspg",
            "rom_stages_to_test" : ["HROM"],
            "ROM":{
                "svd_truncation_tolerance": 1e-16 ,
                "nodal_unknowns": ["VELOCITY_X", "VELOCITY_Y", "PRESSURE"]
            }
        }""")

        rom_manager_object = RomManager(general_rom_manager_parameters=parameters)

        self.assertListEqual(rom_manager_object.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray(), ["ROM","HROM"])
        self.assertListEqual(rom_manager_object.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray(), ["HROM"])
        self.assertEqual(rom_manager_object.general_rom_manager_parameters["projection_strategy"].GetString(), "lspg")
        self.assertEqual(rom_manager_object.general_rom_manager_parameters["ROM"]["svd_truncation_tolerance"].GetDouble(), 1e-16)
        self.assertListEqual(rom_manager_object.general_rom_manager_parameters["ROM"]["nodal_unknowns"].GetStringArray(), ["VELOCITY_X", "VELOCITY_Y", "PRESSURE"])


        #add other test here


    def test_error_raising(self):
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
                "rom_basis_output_name": "RomParameters_test",
                "rom_basis_output_folder": "test_rom_manager"
            }
        }""")
        rom_manager_object = RomManager(general_rom_manager_parameters=parameters)
        flag_change = 'unavailable_flag'
        with self.assertRaises(Exception) as context:
            rom_manager_object._ChangeRomFlags(flag_change)
        err_msg = f'Unknown flag "unavailable_flag" change for RomParameters.json'
        self.assertIn(err_msg, str(context.exception))


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
