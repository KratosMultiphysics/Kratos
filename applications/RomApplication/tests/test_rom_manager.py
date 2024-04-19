import os
import types
import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
#from unittest.mock import patch, MagicMock
from KratosMultiphysics.RomApplication.rom_manager import RomManager

class TestRomManager(KratosUnittest.TestCase):
    def setUp(self):
        # Setup reusable assets here
        self.rom_manager = RomManager()

    def test_initialization(self):
        # Test initialization logic
        self.assertIsNotNone(self.rom_manager.general_rom_manager_parameters)

    # @patch('your_module.KratosMultiphysics.Parameters')
    # def test_load_configuration(self, mock_params):
    #     # Test loading configurations
    #     mock_params.return_value.read.return_value = '{"test": "value"}'
    #     result = self.rom_manager.load_configuration("dummy_path.json")
    #     self.assertEqual(result["test"], "value")

# Run tests
if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
