import numpy as np
import pathlib
import json
import KratosMultiphysics


class RBF_ROM_Interface:
    def __init__(self, mu_train, data_base):
        """
        Initialize the RBF Enhanced ROM Interface.

        Args:
            mu_train: Training parameters.
            data_base: Reference to the database.
        """
        model_name, _ = data_base.get_hashed_file_name_for_table("RBF_Model", mu_train)
        model_path = pathlib.Path(data_base.database_root_directory / 'saved_models' / model_name)

        with open(pathlib.Path(model_path / 'train_config.json'), "r") as config_file:
            model_config = json.load(config_file)

        self.n_inf = int(model_config["modes"][0])
        self.n_sup = int(model_config["modes"][1])
        self.kernel = model_config.get("kernel", "gaussian")
        self.epsilon = model_config.get("epsilon", 1.0)

        # Load phi and singular values
        _, hash_basis = data_base.check_if_in_database("RightBasis", mu_train)
        self.phi = data_base.get_single_numpy_from_database(hash_basis)
        _, hash_sigma = data_base.check_if_in_database("SingularValues_Solution", mu_train)
        self.sigma = data_base.get_single_numpy_from_database(hash_sigma) / np.sqrt(len(mu_train))

        # Load reference snapshot
        self.ref_snapshot = np.zeros(self.phi.shape[0])

    def get_encode_function(self):
        """
        Generate the encoding function for RBF ROM.

        Returns:
            encode_function: A callable encoding function.
        """
        phi_inf = self.phi[:, :self.n_inf]
        ref_snapshot = self.ref_snapshot

        def encode_function(s):
            output_data = (s - ref_snapshot).copy()
            output_data = (phi_inf.T @ output_data).T
            return np.expand_dims(output_data, axis=0), None

        return encode_function

    def get_phi_matrices(self):
        """
        Get the phi matrices for reduced-order modeling.

        Returns:
            Tuple: phi_inf, phi_sup.
        """
        phi_inf = self.phi[:, :self.n_inf]
        phi_sup = self.phi[:, self.n_inf:self.n_sup]

        return KratosMultiphysics.Matrix(phi_inf), KratosMultiphysics.Matrix(phi_sup)

    def get_rbf_params(self):
        """
        Get the RBF parameters.

        Returns:
            Tuple: kernel, epsilon.
        """
        return self.kernel, self.epsilon

    def get_ref_snapshot(self):
        """
        Get the reference snapshot.

        Returns:
            ref_snapshot: Reference snapshot array.
        """
        return self.ref_snapshot
