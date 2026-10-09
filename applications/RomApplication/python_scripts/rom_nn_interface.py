import numpy as np
import pathlib
import json

import KratosMultiphysics


class NN_ROM_Interface():
    def __init__(self, mu_train, data_base):

        model_name, _ = data_base.get_hashed_file_name_for_table("Neural_Network", mu_train)
        self.model_path = pathlib.Path(data_base.database_root_directory / 'saved_nn_models' / model_name)

        _, hash_basis = data_base.check_if_in_database("RightBasis", mu_train)
        self.phi = data_base.get_single_numpy_from_database(hash_basis)

        # Coupled solvers have one network per sub-solver and no global one (see GetSubSolverInterface)
        config_path = self.model_path / 'train_config.json'
        if config_path.exists():
            with open(config_path, "r") as config_file:
                model_config = json.load(config_file)
            _, hash_sigma = data_base.check_if_in_database("SingularValues_Solution", mu_train)
            sigma = data_base.get_single_numpy_from_database(hash_sigma)/np.sqrt(len(mu_train))
            self._SetData(self.phi, sigma, model_config["modes"], self.model_path / 'model_weights.npy')

    def GetSubSolverInterface(self, sub_solver_name, rows):
        """Returns the interface of a coupled sub-solver: its own network and singular values (files prefixed with the
        sub-solver name, as saved by the RomNeuralNetworkTrainer) and the rows of the basis of its unknowns."""
        file_prefix = f"{sub_solver_name}_"
        with open(self.model_path / f'{file_prefix}train_config.json', "r") as config_file:
            model_config = json.load(config_file)
        interface = self.__class__.__new__(self.__class__)
        interface._SetData(
            self.phi[rows, :],
            np.load(self.model_path / f'{file_prefix}SingularValues.npy'),
            model_config["modes"],
            self.model_path / f'{file_prefix}model_weights.npy')
        return interface

    @classmethod
    def FromNumpyFiles(cls, rom_data_folder, modes, file_prefix="", rows=None):
        """Creates the interface without a RomManager database, e.g. for a network trained outside the RomManager.
        The folder must contain:
        - 'RightBasisMatrix.npy': the solution basis, with at least n_sup columns and the rows in the order used for training.
        - '<file_prefix>SingularValues.npy': the singular values used to scale the basis during training (the RomManager uses sigma/sqrt(number_of_training_cases)).
        - '<file_prefix>model_weights.npy': the list of weight matrices of the network (as saved by the RomManager, layers without bias).
        'modes' is [n_inf, n_sup].
        For a coupled sub-solver, 'file_prefix' is '<sub_solver_name>_' and 'rows' are the rows of the basis of its unknowns.
        """
        rom_data_folder = pathlib.Path(rom_data_folder)
        phi = np.load(rom_data_folder / 'RightBasisMatrix.npy')
        interface = cls.__new__(cls)
        interface._SetData(
            phi if rows is None else phi[rows, :],
            np.load(rom_data_folder / f'{file_prefix}SingularValues.npy'),
            modes,
            rom_data_folder / f'{file_prefix}model_weights.npy')
        return interface

    def _SetData(self, phi, sigma, modes, network_weights_path):
        self.n_inf = int(modes[0])
        self.n_sup = int(modes[1])
        self.network_weights_path = network_weights_path
        self.phi = phi
        self.sigma = sigma
        self.ref_snapshot = np.zeros(self.phi.shape[0])

    def get_encode_function(self):

        phi_inf = self.phi[:,:self.n_inf]
        ref_snapshot = self.ref_snapshot

        def encode_function(s):
            output_data=(s-ref_snapshot).copy()
            output_data = (phi_inf.T@output_data).T
            return np.expand_dims(output_data,axis=0), None
        
        return encode_function
    
    def get_phi_matrices(self):
        phi_inf = self.phi[:,:self.n_inf]
        phi_sup = self.phi[:,self.n_inf:self.n_sup]
        sigma_inf = self.sigma[:self.n_inf]
        sigma_sup = self.sigma[self.n_inf:self.n_sup]
        phi_inf_weighted = phi_inf @ np.diag(sigma_inf)
        sigma_inf_inv = np.linalg.inv(np.diag(sigma_inf))
        phi_sup_weighted = phi_sup @ np.diag(sigma_sup)

        return KratosMultiphysics.Matrix(phi_inf), KratosMultiphysics.Matrix(phi_sup_weighted), KratosMultiphysics.Matrix(sigma_inf_inv)
    
    def get_NN_layers(self):
        layers = np.load(self.network_weights_path, allow_pickle=True)
        for layer in layers:
            layer = KratosMultiphysics.Matrix(layer)
        return layers
    
    def get_ref_snapshot(self):
        return self.ref_snapshot