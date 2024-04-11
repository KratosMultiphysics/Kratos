import sqlite3
import hashlib
import json
from pathlib import Path
import numpy as np

class RomDatabase(object):

    def __init__(self, general_rom_manager_parameters, rom_training_parameters, hrom_training_parameters, mu_names):
        self.database_name = "rom_database_sqlite3.db" #TODO This should be made user-defined. Use the same for FOM ROM HROM HHROM??


        self.mu_names = mu_names
        # A single list with the name of each parameter to be used (for reference) in the database, e.g. ["alpha", "permeability", "mach"].
        # self.mu_names are expected to match the position and the number of inputs in the lists mu_train, mu_test, mu_run

        self.general_rom_manager_parameters = general_rom_manager_parameters
        self.rom_training_parameters = rom_training_parameters
        self.hrom_training_parameters = hrom_training_parameters
        self.set_up_or_get_data_base()



    def set_up_or_get_data_base(self):
        rom_output_folder_name = self.rom_training_parameters["Parameters"]["rom_basis_output_folder"].GetString()
        directory = Path(rom_output_folder_name)
        directory.mkdir(parents=True, exist_ok=True)
        self.database_name = directory / "rom_database.db"  # Adjust the database file name as needed
        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()

        # Updated table definitions including new tables
        table_definitions = {
            "FOM": '''CREATE TABLE IF NOT EXISTS FOM
                      (id INTEGER PRIMARY KEY, parameters TEXT, file_name TEXT)''',
            "ROM": '''CREATE TABLE IF NOT EXISTS ROM
                      (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, file_name TEXT)''',
            "HROM": '''CREATE TABLE IF NOT EXISTS HROM
                       (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, tol_res REAL, file_name TEXT)''',
            "HHROM": '''CREATE TABLE IF NOT EXISTS HHROM
                        (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, tol_res REAL, file_name TEXT)''',
            "LeftBasis": '''CREATE TABLE IF NOT EXISTS LeftBasis
                            (id INTEGER PRIMARY KEY, tol_sol REAL, file_name TEXT)''',
            "RightBasis": '''CREATE TABLE IF NOT EXISTS RightBasis
                             (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, file_name TEXT)''',
            "ResidualsProjected": '''CREATE TABLE IF NOT EXISTS ResidualsProjected
                            (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, tol_res REAL, type_of_projection TEXT, file_name TEXT)''',
            "Conditions": '''CREATE TABLE IF NOT EXISTS Conditions
                             (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, tol_res REAL, file_name TEXT)''',
            "ConditionsWeights": '''CREATE TABLE IF NOT EXISTS ConditionsWeights
                                    (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, tol_res REAL, file_name TEXT)''',
            "Elements": '''CREATE TABLE IF NOT EXISTS Elements
                           (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, tol_res REAL, file_name TEXT)''',
            "ElementsWeights": '''CREATE TABLE IF NOT EXISTS ElementsWeights
                                  (id INTEGER PRIMARY KEY, parameters TEXT, tol_sol REAL, tol_res REAL, file_name TEXT)'''
        }

        # Create each table using its definition
        for table_name, table_sql in table_definitions.items():
            try:
                cursor.execute(table_sql)
                conn.commit()
                print(f"Table {table_name} created successfully.")
            except sqlite3.OperationalError as e:
                print(f"Error creating table {table_name}: {e}")

        # Close the database connection
        conn.close()


    def hash_parameters(self, mu):
        # Convert parameters list to a string and encode
        mu_str = '_'.join(map(str, mu))
        return hashlib.sha256(mu_str.encode()).hexdigest()

    def hash_parameters_with_tol(self, mu, tol):
        # Convert the parameters list and the tolerance value to a string
        # Ensure the tolerance is converted to a string with consistent formatting
        mu_str = '_'.join(map(str, mu)) + f"_{tol:.10f}"  # Example with 10 decimal places
        # Generate the hash
        return hashlib.sha256(mu_str.encode()).hexdigest()


    def hash_parameters_with_tol_2_tols(self, mu, tol, tol2):
        # Convert the parameters list to a string and encode
        mu_str = '_'.join(map(str, mu)) + f"_{tol:.10f}_{tol2:.10f}"  # Convert both tol and tol2 to strings with consistent formatting
        # Generate the hash of the combined string including both tolerance values
        return hashlib.sha256(mu_str.encode()).hexdigest()

    def hash_parameters_with_tol_2_tols_and_string(self, mu, tol, tol2, projection):
        # Convert the parameters list to a string and encode
        mu_str = '_'.join(map(str, mu)) + f"_{tol:.10f}_{tol2:.10f}_" +  projection # Convert both tol and tol2 to strings with consistent formatting
        # Generate the hash of the combined string including both tolerance values
        return hashlib.sha256(mu_str.encode()).hexdigest()


    def check_if_mu_already_in_database(self, mu):
        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()
        hash_mu = self.hash_parameters(mu)
        cursor.execute('SELECT COUNT(*) FROM FOM WHERE file_name = ?', (hash_mu,))
        count = cursor.fetchone()[0]
        conn.close()
        return count > 0, hash_mu


    def check_if_rom_already_in_database(self, mu, tol_sol):
        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()
        hash_mu = self.hash_parameters_with_tol(mu, tol_sol)
        cursor.execute('SELECT COUNT(*) FROM ROM WHERE file_name = ?', (hash_mu,))
        count = cursor.fetchone()[0]
        conn.close()
        return count > 0, hash_mu


    def check_if_hrom_already_in_database(self, mu, tol_sol, tol_res):
        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()
        hash_mu = self.hash_parameters_with_tol_2_tols(mu, tol_sol, tol_res)
        cursor.execute('SELECT COUNT(*) FROM HROM WHERE file_name = ?', (hash_mu,))
        count = cursor.fetchone()[0]
        conn.close()
        return count > 0, hash_mu


    def check_if_res_already_in_database(self, mu, tol_sol, tol_res, projection_type):
        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()
        hash_mu = self.hash_parameters_with_tol_2_tols_and_string(mu, tol_sol, tol_res, projection_type)
        cursor.execute('SELECT COUNT(*) FROM ResidualsProjected WHERE file_name = ?', (hash_mu,))
        count = cursor.fetchone()[0]
        conn.close()
        return count > 0, hash_mu


    def check_if_hrom_elems_and_weights_already_in_database(self, *args):#TODO implement the check on the elements and weights
        return False


    def add_FOM_to_database(self, parameters, file_name):
        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()
        parameters_str = self.serialize_mu(parameters)
        cursor.execute('INSERT INTO FOM (parameters, file_name) VALUES (?, ?)',
                       (parameters_str, file_name))
        conn.commit()
        conn.close()

    def add_ROM_to_database(self, parameters, file_name, tol_sol):
        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()
        parameters_str = self.serialize_mu(parameters)
        cursor.execute('INSERT INTO ROM (parameters, tol_sol , file_name) VALUES (?, ?, ?)',
                       (parameters_str, tol_sol, file_name))
        conn.commit()
        conn.close()


    def add_HROM_to_database(self, parameters, file_name, tol_sol, tol_res):
        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()
        parameters_str = self.serialize_mu(parameters)
        cursor.execute('INSERT INTO HROM (parameters, tol_sol , tol_res , file_name) VALUES (?, ?, ?, ?)',
                       (parameters_str, tol_sol, tol_res, file_name))
        conn.commit()
        conn.close()


    def add_ResidualProjected_to_database(self, parameters, file_name, tol_sol, tol_res, type_of_projection):
        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()
        parameters_str = self.serialize_mu(parameters)
        cursor.execute('INSERT INTO ResidualsProjected (parameters, type_of_projection, tol_sol , tol_res , file_name) VALUES (?, ?, ?, ?, ?)',
                       (parameters_str, type_of_projection, tol_sol, tol_res, file_name))
        conn.commit()
        conn.close()

    def serialize_mu(self, parameters):
        return json.dumps(parameters)


    def save_as_npy(self, result, hash_mu):
        rom_output_folder_name = self.rom_training_parameters["Parameters"]["rom_basis_output_folder"].GetString()
        directory = Path(rom_output_folder_name) / 'rom_database'  #TODO hardcoded names here
        file_path = directory / f"{hash_mu}.npy"
        directory.mkdir(parents=True, exist_ok=True)
        np.save(file_path, result)
        return file_path


    def generate_database_summary_file(self):
        table_names = ["FOM", "ROM", "HROM", "HHROM","LeftBasis", "RightBasis", "ResidualsProjected","Conditions","ConditionsWeights", "Elements", "ElementsWeights"]
        rom_output_folder_name = self.rom_training_parameters["Parameters"]["rom_basis_output_folder"].GetString()
        directory = Path(rom_output_folder_name)
        directory.mkdir(parents=True, exist_ok=True)
        summary_path = directory / "rom_database_summary.dat"
        number_of_samples_to_include_in_summary = 5  # Adjust as needed

        with summary_path.open('w') as f:
            conn = sqlite3.connect(self.database_name)
            cursor = conn.cursor()

            for table_name in table_names:
                f.write(f"\nTable Name: {table_name}\n")
                f.write("-" * (10 + len(table_name)) + "\n")

                if table_name =="LeftBasis":

                    # Query the database for a limited number of entries from the current table
                    cursor.execute(f"SELECT tol_sol, file_name FROM {table_name} LIMIT ?", (number_of_samples_to_include_in_summary,))
                    rows = cursor.fetchall()

                    f.write(" | tol_sol | File Name (Hash)\n")
                    f.write("-" * 15 + "\n")
                    # Write each row to the file
                    for tol_sol, file_name in rows:
                        f.write(f"{tol_sol} | {file_name}\n")

                elif table_name =="ResidualsProjected":
                    # Query the database for a limited number of entries from the current table
                    cursor.execute(f"SELECT parameters, tol_sol, tol_res, type_of_projection, file_name FROM {table_name} LIMIT ?", (number_of_samples_to_include_in_summary,))
                    rows = cursor.fetchall()

                    if rows:
                        # Assuming the first row's parameters to figure out the headers
                        sample_params = json.loads(json.loads(rows[0][0]))
                        headers = list(sample_params.keys())
                        f.write(" | ".join(headers) + " |   tol_sol   |  tol_res  |   Projection type  |   File Name (Hash)\n")
                        f.write("-" * ((len(headers)+4) * 15) + "\n")
                        # Write each row to the file
                        for parameters, tol_sol, tol_res, type_of_projection, file_name in rows:
                            params_dict = json.loads(json.loads(parameters))
                            params_str = " | ".join(f"{params_dict.get(name, ''):.3f}" for name in headers)
                            f.write(f"{params_str} | {tol_sol} | {tol_res} | {type_of_projection} | {file_name}\n")
                    else:
                        f.write("No data available.\n")

                else:

                    # Query the database for a limited number of entries from the current table
                    cursor.execute(f"SELECT parameters, file_name FROM {table_name} LIMIT ?", (number_of_samples_to_include_in_summary,))
                    rows = cursor.fetchall()

                    if rows:
                        # Assuming the first row's parameters to figure out the headers
                        sample_params = json.loads(json.loads(rows[0][0]))
                        headers = list(sample_params.keys())
                        f.write(" | ".join(headers) + " | File Name (Hash)\n")
                        f.write("-" * (len(headers) * 15) + "\n")

                        # Write each row to the file
                        for parameters, file_name in rows:
                            params_dict = json.loads(json.loads(parameters))
                            params_str = " | ".join(f"{params_dict.get(name, ''):.3f}" for name in headers)
                            f.write(f"{params_str} | {file_name}\n")
                    else:
                        f.write("No data available.\n")

            print(f"Summary file generated at {summary_path}")

        conn.close()


    def FOM_make_mu_dictionary(self, mu):
        if self.mu_names is None:
            self.mu_names = []
            for i in range(len(mu)):
                self.mu_names.append(f'generic_name_{i}')
        return dict(zip(self.mu_names , mu))

    def LeftBasis_make_dictionary(self, mu):
        keys = ["hashed_name", "svd_tol_solution"]
        return dict(zip(keys, mu))


    def get_snapshots_matrix_from_database(self, mu_list, table_name='FOM'):
        unique_tuples = set(tuple(item) for item in mu_list)
        mu_list_unique = [list(item) for item in unique_tuples] #unique members in mu_lust
        rom_output_folder_name = self.rom_training_parameters["Parameters"]["rom_basis_output_folder"].GetString()
        directory = Path(rom_output_folder_name) / 'rom_database' #TODO hardcoded names here
        directory.mkdir(parents=True, exist_ok=True)
        SnapshotsMatrix = []
        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()
        tol_sol = self.rom_training_parameters["Parameters"]["svd_truncation_tolerance"].GetDouble()
        tol_res =  self.hrom_training_parameters["element_selection_svd_truncation_tolerance"].GetDouble()
        projection_type = self.general_rom_manager_parameters["projection_strategy"].GetString()
        for mu in mu_list_unique:
            serialized_mu = self.serialize_mu(self.FOM_make_mu_dictionary(mu))
            if table_name == 'FOM':
                hash_mu = self.hash_parameters(serialized_mu)
            elif table_name == 'ROM':
                hash_mu = self.hash_parameters_with_tol(serialized_mu, tol_sol)
            elif table_name == 'HROM':
                hash_mu = self.hash_parameters_with_tol_2_tols(serialized_mu, tol_sol, tol_res)
            elif table_name == 'ResidualsProjected':
                hash_mu = self.hash_parameters_with_tol_2_tols_and_string(serialized_mu, tol_sol, tol_res, projection_type)
            cursor.execute(f"SELECT file_name FROM {table_name} WHERE file_name = ?", (hash_mu,)) # Query the database for the file name using the hash
            result = cursor.fetchone()

            if result:
                file_name = result[0]
                file_path = directory / (file_name + '.npy')
                SnapshotsMatrix.append(np.load(file_path))
            else:
                print(f"No entry found for hash {hash_mu}")

        conn.close()

        return np.block(SnapshotsMatrix) if SnapshotsMatrix else None


    def check_if_basis_already_in_database(self, mu_train, tol_sol):
        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()
        serialized_mu_train = self.serialize_entire_mu_train(mu_train)
        hashed_mu_train = self.hash_mu_train(serialized_mu_train, tol_sol)
        # Include both file_name and tol_sol in the WHERE clause
        cursor.execute('SELECT COUNT(*) FROM LeftBasis WHERE file_name = ? AND tol_sol = ?', (hashed_mu_train, tol_sol))
        count = cursor.fetchone()[0]
        conn.close()
        return count > 0


    def add_Basis_to_database(self, mu_train, real_value):
        print(f"Attempting to add tol_sol with value: {real_value}")  # Debugging line

        serialized_mu_train = self.serialize_entire_mu_train(mu_train)
        hashed_mu_train = self.hash_mu_train(serialized_mu_train, real_value)
        file_path = self.save_as_npy(mu_train, hashed_mu_train)  # Save numpy array and get file path

        conn = sqlite3.connect(self.database_name)
        cursor = conn.cursor()

        # Debugging: Print the values before executing the SQL command
        print(f"Inserting into LeftBasis: tol_sol={real_value}, file_name={hashed_mu_train}")

        cursor.execute('INSERT INTO LeftBasis (tol_sol, file_name) VALUES (?, ?)',
                    (real_value, str(hashed_mu_train)))
        conn.commit()
        conn.close()


    def hash_mu_train(self, serialized_mu_train, real_value):
        """
        Generates a hash for the serialized mu_train, including a real value in the hash.
        """
        # Concatenate the serialized mu_train with the real_value (converted to string) for hashing
        # Ensure real_value is converted to a string with consistent formatting
        combined_str = f"{serialized_mu_train}_{real_value:.10f}"  # Example with 10 decimal places
        # Generate the hash of the combined string
        return hashlib.sha256(combined_str.encode()).hexdigest()


    def serialize_entire_mu_train(self, mu_train):
        """
        Serializes the entire mu_train list of lists into a string.
        """
        # Ensure consistent serialization by sorting sublists if needed
        # Here, we assume mu_train's structure is consistent without needing sorting
        serialized_mu_train = json.dumps(mu_train, sort_keys=True)
        return serialized_mu_train