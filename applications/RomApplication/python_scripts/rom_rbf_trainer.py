import numpy as np
import pickle
from sklearn.metrics import mean_squared_error
from itertools import product
import json
import pathlib

class RomRBFTrainer:
    def __init__(self, general_rom_manager_parameters, mu_train, mu_validation, data_base):
        self.general_rom_manager_parameters = general_rom_manager_parameters
        self.mu_train = mu_train
        self.mu_validation = mu_validation
        self.data_base = data_base

        # Extract RBF hyperparameters
        self.rbf_parameters = self.general_rom_manager_parameters["ROM"]["rbf_enhanced_settings"]

        # Generate epsilon samples (logarithmic scale) from epsilon_range
        self.epsilon_samples = np.logspace(
            np.log10(self.rbf_parameters["epsilon_range"][0].GetDouble()), 
            np.log10(self.rbf_parameters["epsilon_range"][1].GetDouble()), 
            int(self.rbf_parameters["num_epsilon_samples"].GetInt())
        )

        # Extract the list of kernels to test
        self.kernel_list = self.rbf_parameters["kernel_list"].GetStringArray()

        # Print for debugging purposes
        print(f"Epsilon Samples: {self.epsilon_samples}")
        print(f"Kernel List: {self.kernel_list}")

    def _get_kernel_function(self, kernel_name):
        """Retrieve the RBF kernel function based on the name."""
        rbf_kernels = {
            'gaussian': lambda r, epsilon: np.exp(-(epsilon * r)**2),
            'imq': lambda r, epsilon: 1.0 / np.sqrt(1 + (epsilon * r)**2),
            # Add other kernels as needed...
        }
        if kernel_name not in rbf_kernels:
            raise ValueError(f"Unsupported kernel: {kernel_name}")
        return rbf_kernels[kernel_name]

    def _get_training_data(self):
        """Load and preprocess training data."""
        # Get the snapshot matrices for training and validation
        S_train = self.data_base.get_snapshots_matrix_from_database(self.mu_train, table_name="FOM")
        S_val = self.data_base.get_snapshots_matrix_from_database(self.mu_validation, table_name="FOM")

        # Load POD basis and singular values
        _, hash_basis = self.data_base.check_if_in_database("RightBasis", self.mu_train)
        phi = self.data_base.get_single_numpy_from_database(hash_basis)

        # Define primary and secondary modes
        self.n_inf_modes = int(self.rbf_parameters["modes"][0].GetInt())
        total_modes = int(self.rbf_parameters["modes"][1].GetInt())
        self.n_sup_modes = int(total_modes - self.n_inf_modes)

        phi_inf = phi[:, :self.n_inf_modes]  # Inferior (primary) modes
        phi_sup = phi[:, self.n_inf_modes:total_modes]  # Superior (secondary) modes

        # Compute reduced coordinates
        q_inf_train = (phi_inf.T @ S_train).T  # Primary reduced coordinates
        q_inf_val = (phi_inf.T @ S_val).T
        q_sup_train = (phi_sup.T @ S_train).T  # Secondary reduced coordinates
        q_sup_val = (phi_sup.T @ S_val).T

        return q_inf_train, q_inf_val, q_sup_train, q_sup_val

    def TrainRBF(self):
        """Train the RBF model by performing a grid search for optimal weights W."""
        # Load training data
        self.q_inf_train, q_inf_val, q_sup_train, q_sup_val = self._get_training_data()

        # Initialize variables for grid search
        best_epsilon = None
        best_kernel_name = None
        lowest_error = np.inf
        best_W = None
        self.regularization = 1e-8

        print("Starting grid search for RBF training...")
        for epsilon, kernel_name in product(self.epsilon_samples, self.kernel_list):
            kernel_func = self._get_kernel_function(kernel_name)

            # Compute distance matrix between training points
            dists_train = np.linalg.norm(self.q_inf_train[:, np.newaxis, :] - self.q_inf_train[np.newaxis, :, :], axis=2)
            Phi_train = kernel_func(dists_train, epsilon) + np.eye(dists_train.shape[0]) * self.regularization

            # Solve for W
            try:
                W = np.linalg.solve(Phi_train, q_sup_train)
            except np.linalg.LinAlgError:
                print(f"LinAlgError at epsilon={epsilon:.5f}, kernel={kernel_name}. Skipping.")
                continue

            # Validation step
            dists_val = np.linalg.norm(q_inf_val[:, np.newaxis, :] - self.q_inf_train[np.newaxis, :, :], axis=2)
            Phi_val = kernel_func(dists_val, epsilon)
            q_sup_pred = Phi_val @ W
            error = mean_squared_error(q_sup_val, q_sup_pred)

            print(f"Epsilon: {epsilon:.5f}, Kernel: {kernel_name}, Validation MSE: {error:.5e}")

            # Update best parameters
            if error < lowest_error:
                lowest_error = error
                self.best_epsilon = epsilon
                self.best_kernel_name = kernel_name
                self.best_W = W.copy()

        # Raise an error if no parameters found
        if self.best_epsilon is None or self.best_kernel_name is None:
            raise ValueError("No suitable epsilon and kernel combination found during grid search.")

        print(f"Best epsilon: {self.best_epsilon:.5f}")
        print(f"Best kernel: {self.best_kernel_name} with Validation MSE: {lowest_error:.5e}")

        # Save the model configuration and weights
        self._save_model()


    def _save_model(self):
        """Save the trained RBF model and relevant parameters."""
        # Generate hashed model name and path
        model_name, _ = self.data_base.get_hashed_file_name_for_table("RBF_Model", self.mu_train)
        model_path = pathlib.Path(self.data_base.database_root_directory / "saved_models" / model_name)
        model_path.mkdir(parents=True, exist_ok=True)

        # Save model weights, kernel name, epsilon, and training coordinates
        model_data = {
            "W": self.best_W,  # Best RBF weights
            "q_inf_train": self.q_inf_train,  # Save training primary coordinates
            "kernel_name": self.best_kernel_name,  # Best kernel
            "epsilon": self.best_epsilon  # Best epsilon
        }
        with open(model_path / "rbf_weights.pkl", "wb") as f:
            pickle.dump(model_data, f)
        print(f"RBF model data saved to: {model_path / 'rbf_weights.pkl'}")

        # Save training configuration (train_config.json) for reference
        training_parameters_dict = {
            "modes": [self.n_inf_modes, self.n_sup_modes],
            "kernel": self.best_kernel_name,
            "epsilon": self.best_epsilon,
            "regularization": self.regularization
        }
        with open(model_path / "train_config.json", "w") as config_file:
            json.dump(training_parameters_dict, config_file, indent=4)
        print(f"RBF training configuration saved to: {model_path / 'train_config.json'}")

    def _GetEvaluationData(self):
        """
        Prepares validation data for evaluating the RBF-enhanced PROM.

        Returns:
            Tuple: Validation snapshots, primary coordinates (q_inf), 
                secondary coordinates (q_sup), phi_inf, phi_sup.
        """
        # Retrieve the validation snapshots from the database
        S_val = self.data_base.get_snapshots_matrix_from_database(self.mu_validation, table_name="FOM")

        # Load the basis (phi) from the database
        _, hash_basis = self.data_base.check_if_in_database("RightBasis", self.mu_train)
        phi = self.data_base.get_single_numpy_from_database(hash_basis)

        # Split the basis into primary and secondary components
        phi_inf = phi[:, :self.n_inf_modes]  # Primary modes
        phi_sup = phi[:, self.n_inf_modes:self.n_inf_modes + self.n_sup_modes]  # Secondary modes

        # Project the snapshots onto the primary and secondary bases
        Q_inf_val = (phi_inf.T @ S_val).T  # Primary reduced coordinates
        Q_sup_val = (phi_sup.T @ S_val).T  # Secondary reduced coordinates

        return S_val, Q_inf_val, Q_sup_val, phi_inf, phi_sup


    def EvaluateRBF(self):
        """
        Evaluates the RBF-enhanced PROM on the validation dataset.
        Computes relative Frobenius error and geometric mean of L2 errors.
        """
        # Load the validation data
        S_val, Q_inf_val, Q_sup_val, phi_inf, phi_sup = self._GetEvaluationData()

        # Load RBF model (weights, kernel, epsilon)
        rbf_model_name, _ = self.data_base.get_hashed_file_name_for_table("RBF_Model", self.mu_train)
        rbf_model_path = pathlib.Path(self.data_base.database_root_directory / "saved_models" / rbf_model_name)
        rbf_weights_path = rbf_model_path / "rbf_weights.pkl"

        with open(rbf_weights_path, 'rb') as f:
            rbf_data = pickle.load(f)

        W = rbf_data['W']  # Global RBF weights
        Q_inf_train = rbf_data['q_inf_train']  # Training primary coordinates
        kernel_name = rbf_data['kernel_name']  # Saved kernel name
        epsilon = rbf_data['epsilon']  # Saved epsilon value
        kernel_func = self._get_kernel_function(kernel_name)

        # Compute RBF kernel matrix between validation and training points
        dists_val = np.linalg.norm(Q_inf_val[:, np.newaxis, :] - Q_inf_train[np.newaxis, :, :], axis=2)
        Phi_val = kernel_func(dists_val, epsilon)

        # Predict secondary coordinates using RBF
        Q_sup_val_pred = Phi_val @ W

        # Reconstruct snapshots using RBF-enhanced PROM
        S_recons_val = phi_sup @ Q_sup_val_pred.T + phi_inf @ Q_inf_val.T

        # Compute POD Sup and POD Inf reconstructions
        S_pod_sup_recons_val = phi_sup @ Q_sup_val.T + phi_inf @ Q_inf_val.T
        S_pod_inf_recons_val = phi_inf @ Q_inf_val.T

        # Compute Relative Frobenius Errors
        err_rel_frob_rbf = np.linalg.norm(S_recons_val - S_val) / np.linalg.norm(S_val)
        err_rel_frob_pod_sup = np.linalg.norm(S_pod_sup_recons_val - S_val) / np.linalg.norm(S_val)
        err_rel_frob_pod_inf = np.linalg.norm(S_pod_inf_recons_val - S_val) / np.linalg.norm(S_val)

        print('RECONSTRUCTION RESULTS, VALIDATION DATASET:')
        print(' - Relative Frobenius error:')
        print('     RBF-PROM: ', err_rel_frob_rbf)
        print('     POD Sup: ', err_rel_frob_pod_sup)
        print('     POD Inf: ', err_rel_frob_pod_inf)

        # Compute Geometric Mean of Relative L2 Errors for each snapshot
        sample_l2_err_list_rbf = [
            np.linalg.norm(S_recons_val[:, i] - S_val[:, i]) / np.linalg.norm(S_val[:, i])
            for i in range(S_val.shape[1])
        ]
        sample_l2_err_list_pod_sup = [
            np.linalg.norm(S_pod_sup_recons_val[:, i] - S_val[:, i]) / np.linalg.norm(S_val[:, i])
            for i in range(S_val.shape[1])
        ]
        sample_l2_err_list_pod_inf = [
            np.linalg.norm(S_pod_inf_recons_val[:, i] - S_val[:, i]) / np.linalg.norm(S_val[:, i])
            for i in range(S_val.shape[1])
        ]

        geom_mean_rbf = np.exp(np.mean(np.log(sample_l2_err_list_rbf)))
        geom_mean_pod_sup = np.exp(np.mean(np.log(sample_l2_err_list_pod_sup)))
        geom_mean_pod_inf = np.exp(np.mean(np.log(sample_l2_err_list_pod_inf)))

        print(' - Relative geometric mean of the L2 error of each snapshot:')
        print('     RBF-PROM: ', geom_mean_rbf)
        print('     POD Sup: ', geom_mean_pod_sup)
        print('     POD Inf: ', geom_mean_pod_inf)

        return err_rel_frob_rbf


