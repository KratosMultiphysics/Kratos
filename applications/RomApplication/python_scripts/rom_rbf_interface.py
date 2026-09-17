import numpy as np
import pathlib

import KratosMultiphysics


class RBF_ROM_Interface():
    def __init__(self, mu_train, data_base):

        model_name, _ = data_base.get_hashed_file_name_for_table("RBF", mu_train)
        model_path=pathlib.Path(data_base.database_root_directory / 'saved_rbf_models' / model_name)

        model_data = np.load(model_path / "model_data.npz", allow_pickle=True)

        self.n_inf = int(model_data["modes"][0])
        self.n_sup = int(model_data["modes"][1])

        self.W_mat = model_data["W"]
        self.kernel_name = model_data["kernel_name"]
        self.centers_mat = model_data["centers_mat"]
        self.kernel_eps = model_data["kernel_eps"]

        _, hash_basis = data_base.check_if_in_database("RightBasis", mu_train)
        self.phi = data_base.get_single_numpy_from_database(hash_basis)
        _, hash_sigma = data_base.check_if_in_database("SingularValues_Solution", mu_train)
        self.sigma =  data_base.get_single_numpy_from_database(hash_sigma)/np.sqrt(len(mu_train))
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
    
    def get_RBF_data_for_kratos(self):
        return {
            "W_mat": KratosMultiphysics.Matrix(self.W_mat),
            "centers_mat": KratosMultiphysics.Matrix(self.centers_mat),
            "kernel_name": self.kernel_name,
            "kernel_eps": self.kernel_eps
        }
    
    def get_ref_snapshot(self):
        return self.ref_snapshot