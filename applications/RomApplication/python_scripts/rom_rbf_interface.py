import numpy as np
import pathlib

import KratosMultiphysics


class RBF_ROM_Interface():
    def __init__(self, mu_train, data_base):

        model_name, _ = data_base.get_hashed_file_name_for_table("RBF", mu_train)
        self.model_path = pathlib.Path(data_base.database_root_directory / 'saved_rbf_models' / model_name)

        _, hash_basis = data_base.check_if_in_database("RightBasis", mu_train)
        self.phi = data_base.get_single_numpy_from_database(hash_basis)

        # Coupled solvers have one RBF model per sub-solver and no global one (see GetSubSolverInterface)
        if (self.model_path / "model_data.npz").exists():
            self._SetDataFromModelFile(self.phi, self.model_path / "model_data.npz")

    def GetSubSolverInterface(self, sub_solver_name, rows):
        """Returns the interface of a coupled sub-solver: its own RBF model (file prefixed with the sub-solver name,
        as saved by the RomRBFTrainer) and the rows of the basis of its unknowns."""
        interface = self.__class__.__new__(self.__class__)
        interface._SetDataFromModelFile(self.phi[rows, :], self.model_path / f"{sub_solver_name}_model_data.npz")
        return interface

    @classmethod
    def FromNumpyFiles(cls, rom_data_folder, rbf_settings, file_prefix="", rows=None):
        """Creates the interface without a RomManager database, e.g. for an RBF model trained outside the RomManager.
        The folder must contain:
        - 'RightBasisMatrix.npy': the solution basis, with at least n_sup columns and the rows in the order used for training.
        - '<file_prefix>SingularValues.npy': the singular values used to scale the basis during training (the RomManager uses sigma/sqrt(number_of_training_cases)).
        - '<file_prefix>rbf_centers.npy': the centers of the RBF, i.e. the scaled inferior reduced coordinates of the training snapshots (number of centers x n_inf).
        - '<file_prefix>rbf_weights.npy': the weights of the RBF (number of centers x (n_sup - n_inf)).
        'rbf_settings' are the Parameters with the "kernel" ("gaussian" or "imq") and its "epsilon".
        For a coupled sub-solver, 'file_prefix' is '<sub_solver_name>_' and 'rows' are the rows of the basis of its unknowns.
        """
        rom_data_folder = pathlib.Path(rom_data_folder)
        phi = np.load(rom_data_folder / 'RightBasisMatrix.npy')
        interface = cls.__new__(cls)
        interface._SetData(
            phi if rows is None else phi[rows, :],
            np.load(rom_data_folder / f'{file_prefix}SingularValues.npy'),
            np.load(rom_data_folder / f'{file_prefix}rbf_weights.npy'),
            np.load(rom_data_folder / f'{file_prefix}rbf_centers.npy'),
            rbf_settings["kernel"].GetString(),
            rbf_settings["epsilon"].GetDouble())
        return interface

    def SaveToNumpyFiles(self, rom_data_folder, file_prefix=""):
        """Saves the files of the RBF model read by FromNumpyFiles (all but the basis) and returns its 'rbf_settings'."""
        rom_data_folder = pathlib.Path(rom_data_folder)
        np.save(rom_data_folder / f'{file_prefix}SingularValues.npy', self.sigma[:self.n_sup])
        np.save(rom_data_folder / f'{file_prefix}rbf_weights.npy', self.W_mat)
        np.save(rom_data_folder / f'{file_prefix}rbf_centers.npy', self.centers_mat)
        return {"kernel" : self.kernel_name, "epsilon" : self.kernel_eps}

    def _SetDataFromModelFile(self, phi, model_file_path):
        model_data = np.load(model_file_path)
        self._SetData(phi, model_data["singular_values"], model_data["W"], model_data["centers_mat"], model_data["kernel_name"][0], model_data["kernel_eps"][0])

    def _SetData(self, phi, sigma, W_mat, centers_mat, kernel_name, kernel_eps):
        # The centers are the inferior reduced coordinates of the training snapshots, and the RBF returns the superior ones
        self.n_inf = centers_mat.shape[1]
        self.n_sup = self.n_inf + W_mat.shape[1]
        self.W_mat = W_mat
        self.centers_mat = centers_mat
        self.kernel_name = str(kernel_name)
        self.kernel_eps = float(kernel_eps)
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
        sigma_inf_inv = np.linalg.inv(np.diag(sigma_inf))
        phi_sup_weighted = phi_sup @ np.diag(sigma_sup)

        return KratosMultiphysics.Matrix(phi_inf), KratosMultiphysics.Matrix(phi_sup_weighted), KratosMultiphysics.Matrix(sigma_inf_inv)

    def get_RBF_data_for_kratos(self):
        kernel_types = {"gaussian": 0, "imq": 1}
        if self.kernel_name not in kernel_types:
            raise Exception(f"The kernel of the RBF model ('{self.kernel_name}') is not one of {list(kernel_types.keys())}.")
        return {
            "W_mat": KratosMultiphysics.Matrix(self.W_mat),
            "centers_mat": KratosMultiphysics.Matrix(self.centers_mat),
            "kernel_type": kernel_types[self.kernel_name],
            "kernel_eps": self.kernel_eps
        }

    def get_ref_snapshot(self):
        return self.ref_snapshot
