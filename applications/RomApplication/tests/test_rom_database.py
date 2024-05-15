import os
import numpy as np
import json
import sqlite3

try:
    import pandas as pd
except ModuleNotFoundError:
    pd = None

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from pathlib import Path
from KratosMultiphysics.RomApplication.rom_manager import RomManager
from KratosMultiphysics.RomApplication.rom_database import RomDatabase


class TestRomDatabase(KratosUnittest.TestCase):

    def setUp(self):
        self.work_folder = "test_rom_database"
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            rm  = RomManager()
            self.parameters = rm.general_rom_manager_parameters
            self.mu_names = ["alpha", "beta", "gamma"]
            self.database_name = "test_rom_database.db"

    def test_database_initialization(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            db = RomDatabase(self.parameters, self.mu_names, self.database_name)
            self.assertTrue(db.database_name.exists(), "Database file was not created.")

    @KratosUnittest.skipIf(pd == None, "this test requires pandas")
    def test_table_creation(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            db = RomDatabase(self.parameters, self.mu_names, self.database_name)
            conn = sqlite3.connect("rom_data/rom_database/"+self.database_name)
            query = "SELECT name FROM sqlite_master WHERE type='table';"
            tables = pd.read_sql_query(query, conn)
            expected_tables = [
                'FOM', 'ROM', 'HROM', 'RightBasis', 'SingularValues_Solution', 'LeftBasis', 'PetrovGalerkinSnapshots',
                'ResidualsProjected', 'SingularValues_Residuals', 'HROM_Elements', 'HROM_Weights', 'Neural_Network'
            ]
            self.assertEqual(sorted(tables['name']), sorted(expected_tables), "Not all expected tables were created.")
            conn.close()

    def test_add_and_retrieve_data(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            db = RomDatabase(self.parameters, self.mu_names, self.database_name)
            mu_test = [0.1, 0.01, 0.001]
            numpy_array = np.array([mu_test])
            db.add_to_database("FOM", mu_test, numpy_array)
            exists, _ = db.check_if_in_database("FOM", mu_test)
            self.assertTrue(exists, "Data was not added correctly to the database.")

    def test_file_operations(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            db = RomDatabase(self.parameters, self.mu_names, self.database_name)
            numpy_array = np.array([1, 2, 3])
            hash_mu, _ = db.get_hashed_mu_for_table("FOM", [1, 2, 3])
            db.save_as_npy(numpy_array, hash_mu)
            file_path = db.npys_directory / f"{hash_mu}.npy"
            self.assertTrue(file_path.exists(), "Numpy file was not saved correctly.")

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            kratos_utilities.DeleteDirectoryIfExisting(Path('./rom_data'))

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
