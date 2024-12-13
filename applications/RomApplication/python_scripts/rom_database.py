import sqlite3
import hashlib
import json
from pathlib import Path
import numpy as np
try:
    import pandas as pd
    from xlsxwriter import Workbook
    missing_pandas = False
    missing_xlsxwriter = False
except ImportError as e:
    if 'pandas' in str(e):
        missing_pandas = True
        missing_xlsxwriter = True
    elif 'xlsxwriter' in str(e):
        missing_xlsxwriter = True
        missing_pandas = False


class RomDatabase(object):

    def __init__(self, general_rom_manager_parameters, mu_names, database_name = "rom_database.db"):
        """
        Initialize the RomDatabase object.

        Args:
            general_rom_manager_parameters: Parameters for the ROM manager.
            mu_names: List of parameter names. A single list with the name of each parameter to be used (for reference) in the database, e.g. ["alpha", "permeability", "mach"].
                      self.mu_names are expected to match the position and the number of inputs in the lists mu_train, mu_test, mu_run
            database_name: Name of the SQLite database file.
        """
        self.general_rom_manager_parameters = general_rom_manager_parameters
        self.database_root_directory = Path(self.general_rom_manager_parameters["ROM"]["rom_basis_output_folder"].GetString() + "/rom_database"  )
        self.database_name = self.database_root_directory / database_name
        self.npys_directory = self.database_root_directory / "npy_files"
        self.xlsx_directory = self.database_root_directory / "xlsx_files"
        self.create_directories([self.database_root_directory, self.npys_directory, self.xlsx_directory])
        self.mu_names = mu_names
        self.set_up_or_get_data_base()


    def create_directories(self, directories):
        """
        Create directories if they do not exist.

        Args:
            directories: List of Path objects representing directories to create.
        """
        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)


    def set_up_or_get_data_base(self):
        """
        Set up the database by creating tables if they do not exist.
        """

        table_definitions = {
            "FOM": '''CREATE TABLE IF NOT EXISTS FOM
                        (id INTEGER PRIMARY KEY, parameters TEXT, file_name TEXT)''',
            "ROM": '''CREATE TABLE IF NOT EXISTS ROM
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL,  type_of_projection TEXT, type_of_decoder TEXT, using_non_converged_sols REAL, file_name TEXT)''',
            "HROM": '''CREATE TABLE IF NOT EXISTS HROM
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, tol_res REAL, type_of_projection TEXT, type_of_decoder TEXT, using_non_converged_sols REAL, file_name TEXT)''',
            "NonconvergedFOM": '''CREATE TABLE IF NOT EXISTS NonconvergedFOM
                        (id INTEGER PRIMARY KEY, parameters TEXT, file_name TEXT)''',
            "NonconvergedROM": '''CREATE TABLE IF NOT EXISTS NonconvergedROM
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, type_of_projection TEXT, type_of_decoder TEXT, using_non_converged_sols REAL, file_name TEXT)''',
            "NonconvergedHROM": '''CREATE TABLE IF NOT EXISTS NonconvergedHROM
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, tol_res REAL, type_of_projection TEXT, type_of_decoder TEXT, using_non_converged_sols REAL, file_name TEXT)''',
            "RightBasis": '''CREATE TABLE IF NOT EXISTS RightBasis
                        (id INTEGER PRIMARY KEY, tol_sol REAL, using_non_converged_sols INTEGER,  file_name TEXT)''',
            "SingularValues_Solution": '''CREATE TABLE IF NOT EXISTS SingularValues_Solution
                        (id INTEGER PRIMARY KEY, tol_sol REAL, using_non_converged_sols INTEGER, file_name TEXT)''',
            "LeftBasis": '''CREATE TABLE IF NOT EXISTS LeftBasis
                        (id INTEGER PRIMARY KEY, tol_sol REAL, type_of_projection TEXT, type_of_decoder TEXT, using_non_converged_sols REAL, basis_strategy TEXT, include_phi INTEGER, tol_pg REAL, solving_technique TEXT, monotonicity_preserving INTEGER , file_name TEXT)''',
            "PetrovGalerkinSnapshots": '''CREATE TABLE IF NOT EXISTS PetrovGalerkinSnapshots
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, type_of_projection TEXT, type_of_decoder TEXT, using_non_converged_sols REAL, basis_strategy TEXT, include_phi INTEGER, tol_pg REAL, solving_technique TEXT, monotonicity_preserving INTEGER , file_name TEXT)''',
            "ResidualsProjected": '''CREATE TABLE IF NOT EXISTS ResidualsProjected
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, tol_res REAL, type_of_projection TEXT, type_of_decoder TEXT, using_non_converged_sols REAL, file_name TEXT)''',
            "SingularValues_Residuals": '''CREATE TABLE IF NOT EXISTS SingularValues_Residuals
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, tol_res REAL, type_of_projection TEXT, type_of_decoder TEXT, using_non_converged_sols REAL, file_name TEXT)''',
            "HROM_Elements": '''CREATE TABLE IF NOT EXISTS HROM_Elements
                        (id INTEGER PRIMARY KEY, tol_sol REAL, tol_res REAL, type_of_projection TEXT, type_of_decoder TEXT, using_non_converged_sols REAL, file_name TEXT)''',
            "HROM_Weights": '''CREATE TABLE IF NOT EXISTS HROM_Weights
                        (id INTEGER PRIMARY KEY, tol_sol REAL, tol_res REAL, type_of_projection TEXT, type_of_decoder TEXT, using_non_converged_sols REAL, file_name TEXT)''',
            "Neural_Network": '''CREATE TABLE IF NOT EXISTS Neural_Network
                        (id INTEGER PRIMARY KEY, tol_sol REAL, type_of_projection TEXT, type_of_decoder TEXT, using_non_converged_sols REAL, modes TEXT, layers_size TEXT, batch_size INTEGER, epochs INTEGER, scheduler TEXT, base_lr REAL, additional_params TEXT, model_number INTEGER, NNgrad_regularisation_weight REAL, file_name TEXT)''',
            "QoI_FOM": '''CREATE TABLE IF NOT EXISTS QoI_FOM
                        (id INTEGER PRIMARY KEY, parameters TEXT, is_active INTEGER , file_name TEXT)''',
            "QoI_ROM": '''CREATE TABLE IF NOT EXISTS QoI_ROM
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL,  type_of_projection TEXT, type_of_decoder TEXT, using_non_converged_sols REAL, is_active INTEGER , file_name TEXT)''',
            "QoI_HROM": '''CREATE TABLE IF NOT EXISTS QoI_HROM
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, tol_res REAL, type_of_projection TEXT, type_of_decoder TEXT, using_non_converged_sols REAL, is_active INTEGER  , file_name TEXT)''',
            "RBF_Model": '''CREATE TABLE IF NOT EXISTS RBF_Model
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, type_of_projection TEXT, type_of_decoder TEXT, using_non_converged_sols REAL, kernel TEXT, epsilon REAL, neighbors INTEGER, file_name TEXT)'''
        }

        self.table_names = table_definitions.keys()
        try:
            with sqlite3.connect(self.database_name) as conn:
                cursor = conn.cursor()
                for table_name, table_sql in table_definitions.items():
                    cursor.execute(table_sql)
                    print(f"Table {table_name} created successfully.")
        except sqlite3.OperationalError as e:
            print(f"Error creating tables: {e}")



    def hash_parameters(self, *args):
        """
        Create a SHA-256 hash of the given parameters.

        Args:
            *args: Parameters to hash.

        Returns:
            str: SHA-256 hash of the parameters.
        """
        params_str = '_'.join(str(arg) if isinstance(arg, str) else f"{arg:.10f}" if isinstance(arg, float) else str(arg) for arg in args)
        return hashlib.sha256(params_str.encode()).hexdigest()

    def get_hashed_file_name_for_table(self, table_name, mu):
        """
        Generate a hashed file name based on the table name and parameters.

        Args:
            table_name: Name of the database table.
            mu: Parameters to include in the hash.

        Returns:
            tuple: Hashed file name and serialized parameters.
        """
        type_of_list = self.identify_list_type(mu)
        if type_of_list=="mu":
            serialized_mu = self.serialize_mu(self.make_mu_dictionary(mu))
        elif type_of_list=="complete_mu":
            serialized_mu = self.serialize_entire_mu_train(mu)
        elif type_of_list=="Empty list":
            serialized_mu = 'No mu provided, running single case as described by the ProjectParameters.json'
        else:
            err_msg = f'Error: {self.identify_list_type(mu)}'
            raise Exception(err_msg)
        tol_sol, tol_res, projection_type, decoder_type, pg_data1_str, pg_data2_bool, pg_data3_double, pg_data4_str, pg_data5_bool, nn_data6_str, nn_data7_str, nn_data8_int, nn_data9_int, nn_data10_str, nn_data11_double, nn_data12_str, nn_data13_int, nn_data14_double, non_converged_fom_14_bool = self.get_curret_params()
        if table_name == 'FOM':
            hash_mu = self.hash_parameters(serialized_mu, table_name)
        elif table_name == 'ROM':
            if decoder_type=="ann_enhanced":
                ann_params = self.hash_parameters(serialized_mu, tol_sol, projection_type, decoder_type, non_converged_fom_14_bool, nn_data6_str, nn_data7_str, nn_data8_int, nn_data9_int, nn_data10_str, nn_data11_double, nn_data12_str, nn_data13_int, nn_data14_double, table_name)
                hash_mu = self.hash_parameters(serialized_mu, tol_sol, projection_type, decoder_type, ann_params ,non_converged_fom_14_bool,table_name)
            else:
                hash_mu = self.hash_parameters(serialized_mu, tol_sol, projection_type, decoder_type, non_converged_fom_14_bool,table_name)
        elif table_name == 'HROM':
            if decoder_type=="ann_enhanced":
                ann_params = self.hash_parameters(serialized_mu, tol_sol, projection_type, decoder_type, non_converged_fom_14_bool, nn_data6_str, nn_data7_str, nn_data8_int, nn_data9_int, nn_data10_str, nn_data11_double, nn_data12_str, nn_data13_int, nn_data14_double, table_name)
                hash_mu = self.hash_parameters(serialized_mu, tol_sol, tol_res, projection_type, decoder_type, ann_params, non_converged_fom_14_bool, table_name)
            else:
                hash_mu = self.hash_parameters(serialized_mu, tol_sol, tol_res, projection_type, decoder_type, non_converged_fom_14_bool, table_name)
        elif table_name == 'NonconvergedFOM':
            hash_mu = self.hash_parameters(serialized_mu, table_name)
        elif table_name == 'NonconvergedROM':
            hash_mu = self.hash_parameters(serialized_mu, tol_sol, projection_type, decoder_type, non_converged_fom_14_bool,table_name)
        elif table_name == 'NonconvergedHROM':
            hash_mu = self.hash_parameters(serialized_mu, tol_sol, tol_res, projection_type, decoder_type, non_converged_fom_14_bool,table_name)
        elif table_name == 'ResidualsProjected':
            hash_mu = self.hash_parameters(serialized_mu, tol_sol, tol_res, projection_type, decoder_type, non_converged_fom_14_bool,table_name)
        elif table_name == 'SingularValues_Residuals':
            hash_mu = self.hash_parameters(serialized_mu, tol_sol, tol_res, projection_type, decoder_type, non_converged_fom_14_bool,table_name)
        elif table_name == 'PetrovGalerkinSnapshots':
            hash_mu = self.hash_parameters(serialized_mu, tol_sol, projection_type, decoder_type, non_converged_fom_14_bool, pg_data1_str,pg_data2_bool,pg_data3_double,pg_data4_str,pg_data5_bool,table_name)
        elif table_name == 'RightBasis':
            hash_mu = self.hash_parameters(serialized_mu, tol_sol, non_converged_fom_14_bool, table_name)
        elif table_name == 'SingularValues_Solution':
            hash_mu = self.hash_parameters(serialized_mu, tol_sol, non_converged_fom_14_bool, table_name)
        #TODO some other params might need to be added the "ann_enhanced" part to notice when they are generated with it. To be added in PR including online AnnEnhanced PROM
        elif table_name == 'LeftBasis':
            hash_mu = self.hash_parameters(serialized_mu, tol_sol, projection_type, decoder_type, non_converged_fom_14_bool, pg_data1_str,pg_data2_bool,pg_data3_double,pg_data4_str,pg_data5_bool,table_name)
        elif table_name == "HROM_Elements":
            hash_mu= self.hash_parameters(serialized_mu, tol_sol,tol_res,projection_type, decoder_type, non_converged_fom_14_bool,table_name)
        elif table_name == "HROM_Weights":
            hash_mu= self.hash_parameters(serialized_mu, tol_sol,tol_res,projection_type, decoder_type, non_converged_fom_14_bool,table_name)
        elif table_name == "Neural_Network":
            hash_mu = self.hash_parameters(serialized_mu, tol_sol, projection_type, decoder_type, non_converged_fom_14_bool, nn_data6_str, nn_data7_str, nn_data8_int, nn_data9_int, nn_data10_str, nn_data11_double, nn_data12_str, nn_data13_int, nn_data14_double, table_name)
        elif table_name == "QoI_FOM":
            hash_mu = self.hash_parameters(serialized_mu, table_name)
        elif table_name == "QoI_ROM":
            if decoder_type=="ann_enhanced":
                ann_params = self.hash_parameters(serialized_mu, tol_sol, projection_type, decoder_type, non_converged_fom_14_bool, nn_data6_str, nn_data7_str, nn_data8_int, nn_data9_int, nn_data10_str, nn_data11_double, nn_data12_str, nn_data13_int, nn_data14_double, table_name)
                hash_mu = self.hash_parameters(serialized_mu, tol_sol,projection_type,decoder_type, ann_params, non_converged_fom_14_bool, table_name)
            else:
                hash_mu = self.hash_parameters(serialized_mu, tol_sol,projection_type,decoder_type, non_converged_fom_14_bool, table_name)
        elif table_name == "QoI_HROM":
            if decoder_type=="ann_enhanced":
                ann_params = self.hash_parameters(serialized_mu, tol_sol, projection_type, decoder_type, non_converged_fom_14_bool, nn_data6_str, nn_data7_str, nn_data8_int, nn_data9_int, nn_data10_str, nn_data11_double, nn_data12_str, nn_data13_int, nn_data14_double, table_name)
                hash_mu = self.hash_parameters(serialized_mu, tol_sol, tol_res,projection_type,decoder_type, ann_params, non_converged_fom_14_bool,table_name)
            else:
                hash_mu = self.hash_parameters(serialized_mu, tol_sol, tol_res,projection_type,decoder_type, non_converged_fom_14_bool,table_name)
        else:
            err_msg = f'Error: table_name: {table_name} not available. Available options are: {", ".join(self.table_names)}'
            raise Exception(err_msg)

        return hash_mu, serialized_mu


    def get_curret_params(self):
        """
        Retrieve current parameters from general_rom_manager_parameters.

        Returns:
            tuple: Various parameters.
        """
        tol_sol = self.general_rom_manager_parameters["ROM"]["svd_truncation_tolerance"].GetDouble()
        tol_res =  self.general_rom_manager_parameters["HROM"]["element_selection_svd_truncation_tolerance"].GetDouble()
        projection_type = self.general_rom_manager_parameters["projection_strategy"].GetString()
        decoder_type = self.general_rom_manager_parameters["type_of_decoder"].GetString()
        pg_data1_str = self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"]["basis_strategy"].GetString()
        pg_data2_bool = self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"]["include_phi"].GetBool()
        pg_data3_double = self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"]["svd_truncation_tolerance"].GetDouble()
        pg_data4_str = self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"]["solving_technique"].GetString()
        pg_data5_bool = self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"]["monotonicity_preserving"].GetBool()
        nn_data6_str = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]["modes"].WriteJsonString()
        nn_data7_str = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]["layers_size"].WriteJsonString()
        nn_data8_int = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]["batch_size"].GetInt()
        nn_data9_int = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]["epochs"].GetInt()
        nn_data10_str = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]["lr_strategy"]["scheduler"].GetString()
        nn_data11_double = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]["lr_strategy"]["base_lr"].GetDouble()
        nn_data12_str = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]["lr_strategy"]["additional_params"].WriteJsonString()
        nn_data13_int = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]["online"]["model_number"].GetInt()
        nn_data14_double = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]["NN_gradient_regularisation_weight"].GetDouble()
        non_converged_fom_14_bool = self.general_rom_manager_parameters["ROM"]["use_non_converged_sols"].GetBool()

        return tol_sol, tol_res, projection_type, decoder_type, pg_data1_str, pg_data2_bool, pg_data3_double, pg_data4_str, pg_data5_bool, nn_data6_str, nn_data7_str, nn_data8_int, nn_data9_int, nn_data10_str, nn_data11_double, nn_data12_str, nn_data13_int, nn_data14_double, non_converged_fom_14_bool



    def check_if_in_database(self, table_name, mu):
        """
        Check if a given set of parameters is already in the specified table.

        Args:
            table_name: Name of the database table.
            mu: Parameters to check.

        Returns:
            tuple: Boolean indicating existence and file name.
        """
        file_name, _ = self.get_hashed_file_name_for_table(table_name, mu)
        with sqlite3.connect(self.database_name) as conn:
            cursor = conn.cursor()
            cursor.execute(f'SELECT COUNT(*) FROM {table_name} WHERE file_name = ?', (file_name,))
            count = cursor.fetchone()[0]
        return count > 0, file_name


    def add_to_database(self, table_name, mu, numpy_array):
        """
        Add a numpy array and its metadata to the specified table in the database.

        Args:
            table_name: Name of the database table.
            mu: Parameters associated with the data.
            numpy_array: Numpy array to store.
        """
        file_name, serialized_mu = self.get_hashed_file_name_for_table(table_name, mu)
        tol_sol, tol_res, projection_type, decoder_type, pg_data1_str, pg_data2_bool, pg_data3_double, pg_data4_str, pg_data5_bool, nn_data6_str, nn_data7_str, nn_data8_int, nn_data9_int, nn_data10_str, nn_data11_double, nn_data12_str, nn_data13_int, nn_data14_double, non_converged_fom_14_bool = self.get_curret_params()

        queries = {
            'FOM': 'INSERT INTO {table} (parameters, file_name) VALUES (?, ?)',
            'ROM': 'INSERT INTO {table} (parameters, tol_sol , type_of_projection, type_of_decoder, using_non_converged_sols, file_name) VALUES (?, ?, ?, ?, ?, ?)',
            'HROM': 'INSERT INTO {table} (parameters, tol_sol , tol_res , type_of_projection, type_of_decoder, using_non_converged_sols, file_name) VALUES (?, ?, ?, ?, ?, ?, ?)',
            'NonconvergedFOM': 'INSERT INTO {table} (parameters, file_name) VALUES (?, ?)',
            'NonconvergedROM': 'INSERT INTO {table} (parameters, tol_sol , type_of_projection, type_of_decoder, using_non_converged_sols, file_name) VALUES (?, ?, ?, ?, ?, ?)',
            'NonconvergedHROM': 'INSERT INTO {table} (parameters, tol_sol , tol_res , type_of_projection, type_of_decoder, using_non_converged_sols, file_name) VALUES (?, ?, ?, ?, ?, ?, ?)',
            'ResidualsProjected': 'INSERT INTO {table} (parameters, type_of_projection, tol_sol , tol_res , type_of_decoder, using_non_converged_sols, file_name) VALUES (?, ?, ?, ?, ?, ?, ?)',
            'SingularValues_Residuals': 'INSERT INTO {table} (parameters, type_of_projection, tol_sol , tol_res , type_of_decoder, using_non_converged_sols, file_name) VALUES (?, ?, ?, ?, ?, ?, ?)',
            'PetrovGalerkinSnapshots': 'INSERT INTO {table} (parameters, tol_sol , type_of_projection, type_of_decoder, using_non_converged_sols, basis_strategy, include_phi, tol_pg, solving_technique, monotonicity_preserving, file_name) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
            'RightBasis': 'INSERT INTO {table} (tol_sol, using_non_converged_sols, file_name) VALUES (?, ?, ?)',
            'SingularValues_Solution': 'INSERT INTO {table} (tol_sol, using_non_converged_sols, file_name) VALUES (?, ?, ?)',
            'LeftBasis': 'INSERT INTO {table} (tol_sol, type_of_projection, type_of_decoder, using_non_converged_sols, basis_strategy, include_phi, tol_pg, solving_technique, monotonicity_preserving, file_name) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
            'HROM_Elements': 'INSERT INTO {table} (tol_sol , tol_res , type_of_projection, type_of_decoder, using_non_converged_sols, file_name) VALUES (?, ?, ?, ?, ?, ?)',
            'HROM_Weights': 'INSERT INTO {table} (tol_sol , tol_res , type_of_projection, type_of_decoder, using_non_converged_sols, file_name) VALUES (?, ?, ?, ?, ?, ?)',
            'Neural_Network': 'INSERT INTO {table} (tol_sol , modes , layers_size, batch_size, epochs, scheduler, base_lr, additional_params, model_number, NNgrad_regularisation_weight, file_name) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
            'QoI_FOM': 'INSERT INTO {table} (parameters, file_name, is_active) VALUES (?, ?, ?)',
            'QoI_ROM': 'INSERT INTO {table} (parameters, tol_sol, type_of_projection, type_of_decoder, using_non_converged_sols, file_name, is_active) VALUES (?, ?, ?, ?, ?, ?, ?)',
            'QoI_HROM': 'INSERT INTO {table} (parameters, tol_sol, tol_res, type_of_projection, type_of_decoder, using_non_converged_sols, file_name, is_active) VALUES (?, ?, ?, ?, ?, ?, ?, ?)'
        }

        if table_name not in queries:
            err_msg = f'Error: table_name: {table_name} not available. Available options are: {", ".join(queries.keys())}'
            raise Exception(err_msg)

        query = queries[table_name].format(table=table_name)

        with sqlite3.connect(self.database_name) as conn:
            cursor = conn.cursor()

            if table_name in ['FOM', 'NonconvergedFOM']:
                cursor.execute(query, (serialized_mu, file_name))
            elif table_name in ['ROM', 'NonconvergedROM']:
                cursor.execute(query, (serialized_mu, tol_sol, projection_type, decoder_type, non_converged_fom_14_bool, file_name))
            elif table_name in ['HROM', 'NonconvergedHROM']:
                cursor.execute(query, (serialized_mu, tol_sol, tol_res, projection_type, decoder_type, non_converged_fom_14_bool, file_name))
            elif table_name in ['ResidualsProjected', 'SingularValues_Residuals']:
                cursor.execute(query, (serialized_mu, projection_type, tol_sol, tol_res, decoder_type, non_converged_fom_14_bool, file_name))
            elif table_name == 'PetrovGalerkinSnapshots':
                cursor.execute(query, (serialized_mu, tol_sol, projection_type, decoder_type, non_converged_fom_14_bool, pg_data1_str, pg_data2_bool, pg_data3_double, pg_data4_str, pg_data5_bool, file_name))
            elif table_name in ['RightBasis', 'SingularValues_Solution']:
                cursor.execute(query, (tol_sol, non_converged_fom_14_bool, file_name))
            elif table_name == 'LeftBasis':
                cursor.execute(query, (tol_sol, projection_type, decoder_type, non_converged_fom_14_bool, pg_data1_str, pg_data2_bool, pg_data3_double, pg_data4_str, pg_data5_bool, file_name))
            elif table_name in ['HROM_Elements', 'HROM_Weights']:
                cursor.execute(query, (tol_sol, tol_res, projection_type, decoder_type, non_converged_fom_14_bool, file_name))
            elif table_name == 'Neural_Network':
                cursor.execute(query, (tol_sol, nn_data6_str, nn_data7_str, nn_data8_int, nn_data9_int, nn_data10_str, nn_data11_double, nn_data12_str, nn_data13_int, nn_data14_double, file_name))
            elif table_name == 'QoI_FOM':
                if len(numpy_array) > 0:
                    cursor.execute(query, (serialized_mu, file_name, True))
                    for key, value in numpy_array.items():
                        cursor.execute(f"PRAGMA table_info({table_name})")
                        columns = [info[1] for info in cursor.fetchall()]
                        if key not in columns:
                            cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {key} TEXT")
                        cursor.execute(f'UPDATE {table_name} SET {key} = ? WHERE file_name = ?', (value, file_name))
                        self.save_as_npy(value, file_name+key)
                else:
                    cursor.execute(query, (serialized_mu, file_name, False))
            elif table_name == 'QoI_ROM':
                if len(numpy_array) > 0:
                    cursor.execute(query, (serialized_mu, tol_sol, projection_type, decoder_type, non_converged_fom_14_bool, file_name, True))
                    for key, value in numpy_array.items():
                        cursor.execute(f"PRAGMA table_info({table_name})")
                        columns = [info[1] for info in cursor.fetchall()]
                        if key not in columns:
                            cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {key} TEXT")
                        cursor.execute(f'UPDATE {table_name} SET {key} = ? WHERE file_name = ?', (value, file_name))
                        self.save_as_npy(value, file_name+key)
                else:
                    cursor.execute(query, (serialized_mu, tol_sol, projection_type, decoder_type, non_converged_fom_14_bool, file_name, False))
            elif table_name == 'QoI_HROM':
                if len(numpy_array) > 0:
                    cursor.execute(query, (serialized_mu, tol_sol, tol_res, projection_type, decoder_type, non_converged_fom_14_bool, file_name, True))
                    for key, value in numpy_array.items():
                        cursor.execute(f"PRAGMA table_info({table_name})")
                        columns = [info[1] for info in cursor.fetchall()]
                        if key not in columns:
                            cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {key} TEXT")
                        cursor.execute(f'UPDATE {table_name} SET {key} = ? WHERE file_name = ?', (value, file_name))
                        self.save_as_npy(value, file_name+key)
                else:
                    cursor.execute(query, (serialized_mu, tol_sol, tol_res, projection_type, decoder_type, non_converged_fom_14_bool, file_name, False))

        if table_name in ["RBF_Model", "Neural_Network", 'QoI_FOM' , 'QoI_ROM' ,'QoI_HROM']:
            pass
        else:
            self.save_as_npy(numpy_array, file_name)



    def identify_list_type(self, lst):
        """
        Identify the type of list provided.

        Args:
            lst: List to identify.

        Returns:
            str: Type of list.
        """
        if not lst:
            return "Empty list"
        if isinstance(lst[0], list):
            if any(isinstance(subitem, list) for item in lst for subitem in item):
                return "Invalid structure - nested lists found"
            return "complete_mu"
        else:
            if any(isinstance(item, list) for item in lst):
                return "Invalid structure - mixed types"
            return "mu"



    def serialize_mu(self, parameters):
        """
        Serialize a parameters dictionary to a JSON string.

        Args:
            parameters: Dictionary of parameters.

        Returns:
            str: JSON string of the parameters.
        """
        return json.dumps(parameters)


    def save_as_npy(self, result, hash_mu):
        """
        Save a numpy array to a .npy file.

        Args:
            result: Numpy array to save.
            hash_mu: Hashed file name.
        """
        file_path = self.npys_directory / f"{hash_mu}.npy"
        np.save(file_path, result)
        print(f"numpy saved to {file_path}")


    def make_mu_dictionary(self, mu):
        """
        Create a dictionary from mu parameters using provided names.

        Args:
            mu: List of parameters.

        Returns:
            dict: Dictionary of parameters.
        """
        if self.mu_names is None:
            self.mu_names = []
            for i in range(len(mu)):
                self.mu_names.append(f'generic_name_{i}')
        return dict(zip(self.mu_names , mu))


    def get_snapshots_matrix_from_database(self, mu_list, table_name='FOM', QoI = None):
        """
        Retrieve a matrix of snapshots from the database based on the given parameters.

        Args:
            mu_list: List of parameter sets.
            table_name: Name of the database table.
            QoI: Quantity of Interest, optional.

        Returns:
            np.array: Matrix of snapshots.
        """
        if mu_list == [None]: #this happens when no mu is passed, the simulation run is the one in the ProjectParameters.json
            mu_list_unique = mu_list
        else:
            unique_tuples = set(tuple(item) for item in mu_list)
            mu_list_unique = [list(item) for item in unique_tuples] #unique members in mu_list
        SnapshotsMatrix = []
        unavailable_cases = []
        with sqlite3.connect(self.database_name) as conn:
            cursor = conn.cursor()
            for mu in mu_list_unique:
                hash_mu, _ = self.get_hashed_file_name_for_table(table_name, mu)
                cursor.execute(f"SELECT file_name FROM {table_name} WHERE file_name = ?", (hash_mu,))
                result = cursor.fetchone()
                if result:
                    file_name = result[0]
                    if QoI is not None:
                        file_name += QoI
                    SnapshotsMatrix.append(self.get_single_numpy_from_database(file_name))
                else:
                    print(f"No entry found for hash {hash_mu}")
                    unavailable_cases.append(mu)

        if unavailable_cases:
            print(f"Retrieved snapshots matrix for {table_name} does not contain {len(unavailable_cases)} cases: {unavailable_cases}")

        return np.block(SnapshotsMatrix) if SnapshotsMatrix else None



    def serialize_entire_mu_train(self, mu_train):
        """
        Serializes the entire mu_train list of lists into a string.
        """
        serialized_mu_train = json.dumps(mu_train, sort_keys=True)
        return serialized_mu_train


    def get_single_numpy_from_database(self, hashed_name):
        """
        Load a numpy array from a .npy file.

        Args:
            hashed_name: Hashed file name.

        Returns:
            np.array: Loaded numpy array.
        """
        return np.load(self.npys_directory / (hashed_name + '.npy'))


    def generate_excel(self, full_tables = False, number_of_terms=5):
        """
        Generate an Excel file summarizing the database contents.

        Args:
            full_tables: Boolean indicating whether to export full tables or a summary.
            number_of_terms: Number of terms to include in the summary.
        """
        if missing_pandas == True:
            print('\x1b[1;31m[MISSING LIBRARY] \x1b[0m'," pandas library not installed. No excel file generated")
        if missing_xlsxwriter == True:
            print('\x1b[1;31m[MISSING LIBRARY] \x1b[0m'," xlsxwriter library not installed. No excel file generated")
        if missing_pandas==False and missing_xlsxwriter==False:
            with sqlite3.connect(self.database_name) as conn:
                query = "SELECT name FROM sqlite_master WHERE type='table';"
                tables = pd.read_sql_query(query, conn)
                if full_tables==False:
                    output_file_path = self.xlsx_directory / 'Summary.xlsx'
                elif full_tables==True:
                    output_file_path = self.xlsx_directory / 'DataBaseCompleteDump.xlsx'
                with pd.ExcelWriter(output_file_path, engine='xlsxwriter') as writer:
                    for table_name in tables['name']:
                        if full_tables==True:
                            df = pd.read_sql_query(f"SELECT * FROM {table_name}", conn)
                        else:
                            df = pd.read_sql_query(f"SELECT * FROM {table_name}  LIMIT {number_of_terms}", conn)
                        df.to_excel(writer, sheet_name=table_name, index=False)
                        print(f"Exported {table_name} to Excel sheet.")
                print(f"Data exported to Excel file {output_file_path} successfully.")

    def generate_database_summary(self):
        """
        Generate a summary Excel file with a limited number of terms.
        """
        number_of_terms = 5 #TODO user-defined
        self.generate_excel(full_tables=False, number_of_terms=number_of_terms)

    def dump_database_as_excel(self):
        """
        Export the entire database to an Excel file.
        """
        self.generate_excel(full_tables=True)