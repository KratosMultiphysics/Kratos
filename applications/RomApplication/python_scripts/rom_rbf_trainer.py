import numpy as np
import pathlib
import json

class RomRBFTrainer(object):

    def __init__(self, general_rom_manager_parameters, mu_train, mu_validation, data_base):

        self.general_rom_manager_parameters = general_rom_manager_parameters
        self.rbf_parameters = self.general_rom_manager_parameters["ROM"]["rbf_enhanced_settings"]
        self.mu_train = mu_train
        self.mu_validation = mu_validation
        self.data_base = data_base

    def _CheckNumberOfModes(self,n_inf,n_sup,n_max):
        if n_inf >= n_max:
            err_msg = f'Specified number of inferior modes ({n_inf}) is higher than or equal to the available ones from the Phi matrix ({n_max}).'
            raise Exception(err_msg)
        elif n_sup > n_max:
            err_msg = f'Specified number of superior modes ({n_sup}) is higher than the available ones from the Phi matrix ({n_max}).'
            raise Exception(err_msg)

    def _GetTrainingData(self, n_inf, n_sup):

        S_train = self.data_base.get_snapshots_matrix_from_database(self.mu_train, table_name=f'FOM')
        S_val = self.data_base.get_snapshots_matrix_from_database(self.mu_validation, table_name=f'FOM')

        _, hash_basis = self.data_base.check_if_in_database("RightBasis", self.mu_train)
        phi = self.data_base.get_single_numpy_from_database(hash_basis)
        _, hash_sigma = self.data_base.check_if_in_database("SingularValues_Solution", self.mu_train)
        sigma_vec =  self.data_base.get_single_numpy_from_database(hash_sigma)/np.sqrt(len(self.mu_train))

        self._CheckNumberOfModes(n_inf,n_sup,sigma_vec.shape[0])

        phisig_inv_inf = np.linalg.inv(np.diag(sigma_vec[:n_inf]))@phi[:,:n_inf].T
        phisig_inv_sup = np.linalg.inv(np.diag(sigma_vec[n_inf:n_sup]))@phi[:,n_inf:n_sup].T
        phisig_inf = phi[:,:n_inf]@np.diag(sigma_vec[:n_inf])
        phisig_sup = phi[:,n_inf:n_sup]@np.diag(sigma_vec[n_inf:n_sup])

        Q_inf_train = (phisig_inv_inf@S_train).T
        Q_inf_val = (phisig_inv_inf@S_val).T
        Q_sup_train = (phisig_inv_sup@S_train).T
        Q_sup_val = (phisig_inv_sup@S_val).T
        Q_inf_train_original = Q_inf_train.copy()

        if self.general_rom_manager_parameters["ROM"]["use_non_converged_sols"].GetBool():
            #fetching nonconverged sols for enlarging training samples in ann enhanced prom
            data = self.data_base.get_snapshots_matrix_from_database(self.mu_train, table_name='NonconvergedFOM') #TODO this might be too large. Add partitioned approached or a limit size
            Q_inf_train = np.r_[Q_inf_train, (phisig_inv_inf@data).T]
            Q_sup_train = np.r_[Q_sup_train, (phisig_inv_sup@data).T]

        phisig_norm_matrix = phisig_sup.T @ phisig_sup

        rescaling_factor = np.mean(np.square((phisig_inf@Q_inf_train_original.T)-S_train))
        rescaling_factor *= S_train.shape[0]/Q_inf_train_original.shape[1]

        return Q_inf_train, Q_inf_val, Q_sup_train, Q_sup_val, phisig_norm_matrix, rescaling_factor

    def _GetEvaluationData(self, n_inf, n_sup):

        S_val = self.data_base.get_snapshots_matrix_from_database(self.mu_validation, table_name=f'FOM')

        _, hash_basis = self.data_base.check_if_in_database("RightBasis", self.mu_train)
        phi = self.data_base.get_single_numpy_from_database(hash_basis)
        _, hash_sigma = self.data_base.check_if_in_database("SingularValues_Solution", self.mu_train)
        sigma_vec =  self.data_base.get_single_numpy_from_database(hash_sigma)/np.sqrt(len(self.mu_train))

        phisig_inv_inf = np.linalg.inv(np.diag(sigma_vec[:n_inf]))@phi[:,:n_inf].T
        phisig_inv_sup = np.linalg.inv(np.diag(sigma_vec[n_inf:n_sup]))@phi[:,n_inf:n_sup].T
        phisig_inf = phi[:,:n_inf]@np.diag(sigma_vec[:n_inf])
        phisig_sup = phi[:,n_inf:n_sup]@np.diag(sigma_vec[n_inf:n_sup])

        Q_inf_val = (phisig_inv_inf@S_val).T
        Q_sup_val = (phisig_inv_sup@S_val).T

        return S_val, Q_inf_val, Q_sup_val, phisig_inf, phisig_sup

    def _GetRBFKernelFunctions(self):
        # Kernel definitions
        def gaussian_rbf(r, epsilon):
            return np.exp(-(epsilon * r) ** 2)

        def inverse_multiquadric_rbf(r, epsilon):
            return 1.0 / np.sqrt(1.0 + (epsilon * r) ** 2)

        return {"gaussian": gaussian_rbf, "imq": inverse_multiquadric_rbf}

    def _GetPairwiseDistanceFunction(self):
        def pairwise_dist(XA, XB):
            """This computes the distances from each point in XA to each point in XB and returns them in format:
            [[XA1-XB1,XA1-XB2,XA1-XB3,...],[XA2-XB1,XA2-XB2,XA2-XB3,...],...]
            """
            return np.linalg.norm(XA[:, None, :] - XB[None, :, :], axis=2)
        return pairwise_dist

    def _GetRBFFunction(self, W_mat, centers_mat, kernel_name, kernel_epsilon):
        pairwise_distance = self._GetPairwiseDistanceFunction()
        kernel_fun = self._GetRBFKernelFunctions()[kernel_name]
        def rbf(x_mat):
            d = pairwise_distance(x_mat, centers_mat)
            return kernel_fun(d,kernel_epsilon) @ W_mat
        return rbf

    def TrainRBF(self):

        rbf_kernels = self._GetRBFKernelFunctions()
        pairwise_dist = self._GetPairwiseDistanceFunction()

        svd_rcond = 1e-12 # SVD cutoff (only used if rbf_solver is "svd")
        lambda_values = np.logspace(-12, -8, 3) # Ridge grid (only used if rbf_solver is "ridge")

        def rel_fro_pct(A, B):
            return 100.0 * (np.linalg.norm(A - B, ord="fro") / (np.linalg.norm(A, ord="fro") + 1e-14))

        # Possible epsilon values
        epsilon_values = np.logspace(np.log10(0.2), np.log10(5.0), 15)

        # Possible solvers:
        def svd_solve_full_rank(Phi, Y, rcond=1e-12):
            U, s, Vt = np.linalg.svd(Phi, full_matrices=False)
            if s.size == 0 or s[0] <= 0.0:
                return np.zeros((Phi.shape[1], Y.shape[1]), dtype=Y.dtype)
            tol = rcond * s[0]
            s_inv = np.where(s > tol, 1.0 / s, 0.0)
            return Vt.T @ (s_inv[:, None] * (U.T @ Y))

        def ridge_solve(Phi, Y, lam):
            n = Phi.shape[0]
            A = Phi + lam * np.eye(n)
            return np.linalg.solve(A, Y)

        # Get training data
        rbf_training_parameters = self.rbf_parameters
        
        n_inf = int(rbf_training_parameters['modes'].GetVector()[0])
        n_sup = int(rbf_training_parameters['modes'].GetVector()[1])
        Q_inf_train, Q_inf_val, Q_sup_train, Q_sup_val, phisig_norm_matrix, rescaling_factor = self._GetTrainingData(n_inf, n_sup) # Check that Q_inf_tran, etc are indeed (num. snapshots) x (num. modes), otherwise they should be transposed

        rbf_solver = rbf_training_parameters['rbf_solver'].GetString()

        # Get pairwise distances
        d_tr = pairwise_dist(Q_inf_train, Q_inf_train)
        d_va = pairwise_dist(Q_inf_val, Q_inf_train)

        # Grid search
        best = {"score": np.inf, "kernel": None, "eps": None, "lam": None}

        for kernel_name, kernel_fun in rbf_kernels.items():
            for eps in epsilon_values:
                kernel_result_train = kernel_fun(d_tr, eps)
                kernel_result_val = kernel_fun(d_va, eps)

                if rbf_solver == "svd":
                    W = svd_solve_full_rank(kernel_result_train, Q_sup_train, rcond=svd_rcond)
                    Q_sup_pred = kernel_result_val @ W
                    score = rel_fro_pct(Q_sup_val, Q_sup_pred)
                    print(f"[STAGE3] kernel={kernel_name:8s} eps={eps:.4e} -> Val RPE={score:.2f}%")
                    if score < best["score"]:
                        best.update(score=score, kernel=kernel_name, eps=float(eps), lam=None)

                elif rbf_solver == "ridge":
                    for lam in lambda_values:
                        W = ridge_solve(kernel_result_train, Q_sup_train, lam)
                        Q_sup_pred = kernel_result_val @ W
                        score = rel_fro_pct(Q_sup_val, Q_sup_pred)
                        print(f"[STAGE3] kernel={kernel_name:8s} eps={eps:.4e} lam={lam:.1e} -> Val RPE={score:.2f}%")
                        if score < best["score"]:
                            best.update(score=score, kernel=kernel_name, eps=float(eps), lam=float(lam))
                else:
                    raise ValueError("'rbf_solver' in configuration must be 'svd' or 'ridge'")

        print("\n[STAGE3] BEST hyperparameters:")
        if rbf_solver == "ridge":
            print(f"        kernel={best['kernel']}, eps={best['eps']:.6e}, lam={best['lam']:.3e}, Val RPE={best['score']:.2f}%")
        else:
            print(f"        kernel={best['kernel']}, eps={best['eps']:.6e}, Val RPE={best['score']:.2f}%")

        # -----------------------------
        # Save model (NOW includes q_p_train_scaled)
        # -----------------------------

        # Things to save for inference:
        # - Weights vector (W)
        # - Centers (Q_inf_train)
        # - Kernel type (kernel_name)
        # - Kernel epsilon
        # - n_inf
        # - n_sup
        # - solver

        model_name, _ = self.data_base.get_hashed_file_name_for_table("RBF", self.mu_train)
        model_path=pathlib.Path(self.data_base.database_root_directory / 'saved_rbf_models' / model_name)
        model_path.mkdir(parents=True, exist_ok=True)

        np.savez(
            model_path / "model_data.npz",
            kernel_name=np.array([best["kernel"]], dtype=object),
            kernel_eps=np.array([best["eps"]], dtype=float),
            W=W,
            centers_mat=Q_inf_train,
            modes = np.array([n_inf,n_sup], dtype=int),
            solver=np.array([rbf_solver], dtype=object)
        )


    def EvaluateRBF(self):

        model_name, _ = self.data_base.get_hashed_file_name_for_table("RBF", self.mu_train)
        model_path=pathlib.Path(self.data_base.database_root_directory / 'saved_rbf_models' / model_name)

        model_data = np.load(model_path / "model_data.npz", allow_pickle=True)

        W_mat = model_data["W"]
        centers_mat = model_data["centers_mat"]
        kernel_name = model_data["kernel_name"][0]
        kernel_eps = model_data["kernel_eps"][0]
        n_inf = model_data['modes'][0]
        n_sup = model_data['modes'][1]

        rbf = self._GetRBFFunction(W_mat, centers_mat, kernel_name, kernel_eps)

        S_val, Q_inf_val, Q_sup_val, phisig_inf, phisig_sup = self._GetEvaluationData(n_inf, n_sup)

        print('RECONSTRUCTION RESULTS, VALIDATION DATASET:')
        print(' - Relative Frobenius error:')

        S_recons_val = phisig_sup@rbf(Q_inf_val).T+phisig_inf@Q_inf_val.T
        err_rel_recons = np.linalg.norm(S_recons_val-S_val)/np.linalg.norm(S_val)
        print('     RBF-PROM: ', err_rel_recons)

        S_pod_sup_recons_val = phisig_sup@Q_sup_val.T + phisig_inf@Q_inf_val.T
        print('     POD Sup: ', np.linalg.norm(S_pod_sup_recons_val-S_val)/np.linalg.norm(S_val))

        S_pod_inf_recons_val = phisig_inf@Q_inf_val.T
        print('     POD Inf: ', np.linalg.norm(S_pod_inf_recons_val-S_val)/np.linalg.norm(S_val))

        print(' - Relative geometric mean of the L2 error of each snapshot:')
        sample_l2_err_list=[]
        for i in range(S_recons_val.shape[1]):
            sample_l2_err_list.append(np.linalg.norm(S_recons_val[:,i]-S_val[:,i])/np.linalg.norm(S_val[:,i]))
        print('     RBF-PROM: ', np.linalg.norm(np.exp(np.mean(np.log(sample_l2_err_list)))))

        sample_l2_err_list_pod_sup=[]
        for i in range(S_pod_sup_recons_val.shape[1]):
            sample_l2_err_list_pod_sup.append(np.linalg.norm(S_pod_sup_recons_val[:,i]-S_val[:,i])/np.linalg.norm(S_val[:,i]))
        print('     POD Sup: ', np.linalg.norm(np.exp(np.mean(np.log(sample_l2_err_list_pod_sup)))))

        sample_l2_err_list_pod_inf=[]
        for i in range(S_pod_inf_recons_val.shape[1]):
            sample_l2_err_list_pod_inf.append(np.linalg.norm(S_pod_inf_recons_val[:,i]-S_val[:,i])/np.linalg.norm(S_val[:,i]))
        print('     POD Inf: ', np.linalg.norm(np.exp(np.mean(np.log(sample_l2_err_list_pod_inf)))))

        return err_rel_recons