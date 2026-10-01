import numpy as np
import pathlib
import json
import sqlite3

import tensorflow as tf
from keras.models import Model
from keras import layers
from keras.optimizers import AdamW
from keras.callbacks import LearningRateScheduler
from keras.utils import set_random_seed as keras_set_random_seed
from keras import metrics

tf.keras.backend.set_floatx('float64')

class ANNPROM_Keras_Model(Model):

    def __init__(self, phisig_norm_matrix, rescaling_factor, w_gradNN, *args, **kwargs):
        super(ANNPROM_Keras_Model,self).__init__(*args, **kwargs)

        self.run_eagerly = False

        self.phisig_norm_matrix = phisig_norm_matrix
        self.w_gradNN = w_gradNN
        self.rescaling_factor_x = rescaling_factor

        self.loss_tracker = metrics.Mean(name="loss")
        self.loss_x_tracker = metrics.Mean(name="loss_x")
        self.loss_gradNN_tracker = metrics.Mean(name="loss_grad")


    def train_step(self,data):
        input_batch, x_true_batch = data # target_aux is the reference force or residual, depending on the settings

        with tf.GradientTape() as tape_d:

            with tf.GradientTape() as tape_e:
                tape_e.watch(input_batch)
                x_pred_batch = self(input_batch, training=True)

            mse_gradNN = tf.math.reduce_mean(tf.math.square(tape_e.gradient(x_pred_batch, input_batch)))

            x_diff_batch=x_true_batch-x_pred_batch
            mse_x = tf.math.reduce_mean(tf.math.multiply(x_diff_batch,tf.transpose(tf.matmul(self.phisig_norm_matrix,tf.transpose(x_diff_batch)))))
            mse_x /= self.rescaling_factor_x

            mse = mse_x + self.w_gradNN * mse_gradNN

        gradients = tape_d.gradient(mse, self.trainable_variables)

        self.optimizer.apply_gradients(zip(gradients, self.trainable_variables))

        # Compute our own metrics
        self.loss_tracker.update_state(mse)
        self.loss_x_tracker.update_state(mse_x)
        self.loss_gradNN_tracker.update_state(mse_gradNN)
        return {"loss": self.loss_tracker.result(), "loss_x": self.loss_x_tracker.result(), "loss_grad": self.loss_gradNN_tracker.result()}

    def test_step(self, data):
        input_batch, x_true_batch = data

        with tf.GradientTape() as tape_e:
            tape_e.watch(input_batch)
            x_pred_batch = self(input_batch, training=True)

        mse_gradNN = tf.math.reduce_mean(tf.math.square(tape_e.gradient(x_pred_batch, input_batch)))

        x_diff_batch=x_true_batch-x_pred_batch
        mse_x = tf.math.reduce_mean(tf.math.multiply(x_diff_batch,tf.transpose(tf.matmul(self.phisig_norm_matrix,tf.transpose(x_diff_batch)))))
        mse_x /= self.rescaling_factor_x

        mse = mse_x + self.w_gradNN * mse_gradNN

        # Compute our own metrics
        self.loss_tracker.update_state(mse)
        self.loss_x_tracker.update_state(mse_x)
        self.loss_gradNN_tracker.update_state(mse_gradNN)
        return {"loss": self.loss_tracker.result(), "loss_x": self.loss_x_tracker.result(), "loss_grad": self.loss_gradNN_tracker.result()}

    @property
    def metrics(self):
        # We list our `Metric` objects here so that `reset_states()` can be
        # called automatically at the start of each epoch
        # or at the start of `evaluate()`.
        # If you don't implement this property, you have to call
        # `reset_states()` yourself at the time of your choosing.

        return [self.loss_tracker, self.loss_x_tracker, self.loss_gradNN_tracker]


class RomNeuralNetworkTrainer(object):
    def __init__(self, general_rom_manager_parameters, mu_train, mu_validation, data_base, solver_info=None):
        self.general_rom_manager_parameters = general_rom_manager_parameters
        self.nn_parameters = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]
        self.mu_train = mu_train
        self.mu_validation = mu_validation
        self.data_base = data_base
        self.solver_info = solver_info

        # Create a prefix if this is a coupled solver (e.g., "fluid_solver_")
        self.file_prefix = ""
        if self.solver_info is not None:
            self.file_prefix = self.solver_info["sub_solver_name"].GetString() + "_"

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

        if getattr(self, "solver_info", None) is not None:
            # 1. Match Kratos C++ Alphabetical Sorting
            global_unknowns = self.general_rom_manager_parameters["ROM"]["nodal_unknowns"].GetStringArray()
            global_unknowns.sort()
            num_vars = len(global_unknowns)
            num_nodes = S_train.shape[0] // num_vars

            physics_unknowns = self.solver_info["nodal_unknowns"].GetStringArray()
            physics_unknowns.sort()

            local_indices = [global_unknowns.index(var) for var in physics_unknowns]
            row_indices = [node_idx * num_vars + loc_idx
                           for node_idx in range(num_nodes)
                           for loc_idx in local_indices]

            # 2. Slice the matrices
            S_train = S_train[row_indices, :]
            S_val = S_val[row_indices, :]
            phi = phi[np.ix_(row_indices, range(phi.shape[1]))]

            # 3. Strip zero-padding to find TRUE available modes
            valid_modes = 0
            for col in range(phi.shape[1]):
                if np.any(np.abs(phi[:, col]) > 1e-12):
                    valid_modes += 1
                else:
                    break

            if n_sup > valid_modes:
                print(f"\n[WARNING] Capping n_sup from {n_sup} to {valid_modes} to avoid zero-padding.\n")
                n_sup = valid_modes

            # 4. Recompute the EXACT True Sigma
            phi = phi[:, :n_sup]
            sigma_vec = np.linalg.norm(phi.T @ S_train, axis=1) / np.sqrt(len(self.mu_train))
        else:
            _, hash_sigma = self.data_base.check_if_in_database("SingularValues_Solution", self.mu_train)
            sigma_vec = self.data_base.get_single_numpy_from_database(hash_sigma) / np.sqrt(len(self.mu_train))

        self._CheckNumberOfModes(n_inf, n_sup, sigma_vec.shape[0])

        # Kept to be saved with the network: the online decoder must use the same scaling
        self.sigma_vec = sigma_vec[:n_sup]

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
            data = self.data_base.get_snapshots_matrix_from_database(self.mu_train, table_name='NonconvergedFOM')
            if getattr(self, "solver_info", None) is not None:
                data = data[row_indices, :]
            Q_inf_train = np.r_[Q_inf_train, (phisig_inv_inf@data).T]
            Q_sup_train = np.r_[Q_sup_train, (phisig_inv_sup@data).T]

        phisig_norm_matrix = phisig_sup.T @ phisig_sup
        rescaling_factor = np.mean(np.square((phisig_inf@Q_inf_train_original.T)-S_train))
        rescaling_factor *= S_train.shape[0]/Q_inf_train_original.shape[1]

        # IMPORTANT: Returning capped n_sup!
        return Q_inf_train, Q_inf_val, Q_sup_train, Q_sup_val, phisig_norm_matrix, rescaling_factor, n_sup

    def _GetEvaluationData(self, model_properties):
        S_val = self.data_base.get_snapshots_matrix_from_database(self.mu_validation, table_name=f'FOM')

        n_inf = model_properties['modes'][0]
        n_sup = model_properties['modes'][1]

        _, hash_basis = self.data_base.check_if_in_database("RightBasis", self.mu_train)
        phi = self.data_base.get_single_numpy_from_database(hash_basis)

        if getattr(self, "solver_info", None) is not None:
            # 1. Match Kratos C++ Alphabetical Sorting
            global_unknowns = self.general_rom_manager_parameters["ROM"]["nodal_unknowns"].GetStringArray()
            global_unknowns.sort()
            num_vars = len(global_unknowns)
            num_nodes = S_val.shape[0] // num_vars

            physics_unknowns = self.solver_info["nodal_unknowns"].GetStringArray()
            physics_unknowns.sort()

            local_indices = [global_unknowns.index(var) for var in physics_unknowns]
            row_indices = [node_idx * num_vars + loc_idx
                           for node_idx in range(num_nodes)
                           for loc_idx in local_indices]

            # 2. Slice Matrices
            S_val = S_val[row_indices, :]
            phi = phi[np.ix_(row_indices, range(n_sup))]

            # 3. Recompute Exact Sigma
            S_train = self.data_base.get_snapshots_matrix_from_database(self.mu_train, table_name=f'FOM')
            S_train = S_train[row_indices, :]
            sigma_vec = np.linalg.norm(phi.T @ S_train, axis=1) / np.sqrt(len(self.mu_train))
        else:
            _, hash_sigma = self.data_base.check_if_in_database("SingularValues_Solution", self.mu_train)
            sigma_vec =  self.data_base.get_single_numpy_from_database(hash_sigma)/np.sqrt(len(self.mu_train))

        phisig_inv_inf = np.linalg.inv(np.diag(sigma_vec[:n_inf]))@phi[:,:n_inf].T
        phisig_inv_sup = np.linalg.inv(np.diag(sigma_vec[n_inf:n_sup]))@phi[:,n_inf:n_sup].T
        phisig_inf = phi[:,:n_inf]@np.diag(sigma_vec[:n_inf])
        phisig_sup = phi[:,n_inf:n_sup]@np.diag(sigma_vec[n_inf:n_sup])

        Q_inf_val = (phisig_inv_inf@S_val).T
        Q_sup_val = (phisig_inv_sup@S_val).T

        return S_val, Q_inf_val, Q_sup_val, phisig_inf, phisig_sup

    def _SelectScheduler(self, strategy_name, base_lr, additional_params):

        def lr_const_scheduler(epoch, lr):
            new_lr= base_lr
            return new_lr

        def lr_steps_scheduler(epoch, lr):
            if epoch==0:
                lr=base_lr
            elif epoch%additional_params[2]==0:
                lr/=additional_params[1]
            if lr<=additional_params[0]:
                lr = additional_params[0]
            return lr

        def lr_sgdr_scheduler(epoch, lr):
            if epoch==0:
                new_lr=base_lr
            else:
                max_lr=base_lr
                min_lr=additional_params[0]
                cycle_length=additional_params[2]
                scale_factor=additional_params[1]
                cycle=np.floor(epoch/cycle_length)
                x=np.abs(epoch/cycle_length-cycle)
                new_lr = min_lr+(max_lr-min_lr)*0.5*(1+np.cos(x*np.pi))/scale_factor**cycle
            return new_lr

        schedulers_dict={"const": lr_const_scheduler, "steps": lr_steps_scheduler, "sgdr": lr_sgdr_scheduler}

        return schedulers_dict[strategy_name]

    def _DefineNetwork(self, n_inf, n_sup, layers_size, phisig_norm_matrix=None, rescaling_factor=None, w_gradNN=None):
        input_layer=layers.Input((n_inf,), dtype=tf.float64)
        layer_out=input_layer
        for layer_size in layers_size:
            layer_out=layers.Dense(int(layer_size), 'elu', use_bias=False, kernel_initializer="he_normal", dtype=tf.float64)(layer_out)
        output_layer=layers.Dense(n_sup-n_inf, 'linear', use_bias=False, kernel_initializer="he_normal", dtype=tf.float64)(layer_out)

        network=ANNPROM_Keras_Model(phisig_norm_matrix, rescaling_factor, w_gradNN, input_layer, output_layer)
        return network

    def _SaveWeightsKratosFormat(self, network, weights_path):
        layers=[]
        for layer in network.trainable_variables:
            layers.append(layer.numpy())

        np.save(weights_path, np.array(layers, dtype=object), allow_pickle=True)



    def TrainNetwork(self, seed=None):

        if seed is not None:
            keras_set_random_seed(seed)


        nn_training_parameters = self.nn_parameters

        n_inf = int(nn_training_parameters['modes'].GetVector()[0])
        n_sup = int(nn_training_parameters['modes'].GetVector()[1])
        layers_size = nn_training_parameters['layers_size'].GetVector()
        batch_size = nn_training_parameters['batch_size'].GetInt()
        w_gradNN = nn_training_parameters['NN_gradient_regularisation_weight'].GetDouble()
        lr_scheme = nn_training_parameters['lr_strategy']['scheduler'].GetString()
        base_lr = nn_training_parameters['lr_strategy']['base_lr'].GetDouble()
        lr_additional_params = nn_training_parameters['lr_strategy']['additional_params'].GetVector()
        epochs = nn_training_parameters['epochs'].GetInt()

        # Always use "Neural_Network" for the DB table
        model_name, _ = self.data_base.get_hashed_file_name_for_table("Neural_Network", self.mu_train)
        model_path = pathlib.Path(self.data_base.database_root_directory / 'saved_nn_models' / model_name)

        # exist_ok=True is critical here so the second loop doesn't crash
        model_path.mkdir(parents=True, exist_ok=True)


        Q_inf_train, Q_inf_val, Q_sup_train, Q_sup_val, phisig_norm_matrix, rescaling_factor, n_sup = self._GetTrainingData(n_inf, n_sup)

        network = self._DefineNetwork(n_inf, n_sup, layers_size, phisig_norm_matrix, rescaling_factor, w_gradNN)

        # def scaled_phinorm_mse_loss(y_true, y_pred):
        #     y_diff=y_true-y_pred
        #     mse = tf.math.reduce_mean(tf.math.multiply(y_diff,tf.transpose(tf.matmul(phisig_norm_matrix,tf.transpose(y_diff)))))
        #     mse /= rescaling_factor

        #     return mse

        network.compile(AdamW(epsilon=1e-17), run_eagerly=False)
        network.summary()

        callbacks = [LearningRateScheduler(self._SelectScheduler(lr_scheme, base_lr, lr_additional_params), verbose=0)]

        # history = network.fit(Q_inf_train, Q_sup_train, batch_size=batch_size, epochs=epochs, validation_data=(Q_inf_val,Q_sup_val), shuffle=False, validation_batch_size=1, callbacks=callbacks)
        history = network.fit(Q_inf_train, Q_sup_train, batch_size=batch_size, epochs=epochs, validation_data=(Q_inf_val,Q_sup_val), shuffle=True, callbacks=callbacks)


        training_parameters_dict = {
            "modes":[n_inf,n_sup],
            "layers_size":list(layers_size),
            "batch_size":batch_size,
            "epochs":epochs,
            "NN_gradient_regularisation_weight": w_gradNN,
            "lr_strategy":{
                "scheduler": lr_scheme,
                "base_lr": base_lr,
                "additional_params": list(lr_additional_params)
            }
        }
        # Apply the prefix to the saved files
        with open(str(model_path) + f"/{self.file_prefix}train_config.json", "w") as ae_config_json_file:
            json.dump(training_parameters_dict, ae_config_json_file)

        network.save_weights(str(model_path) + f"/{self.file_prefix}model.weights.h5")
        with open(str(model_path) + f"/{self.file_prefix}history.json", "w") as history_file:
            json.dump(str(history.history), history_file)

        self._SaveWeightsKratosFormat(network, str(model_path) + f"/{self.file_prefix}model_weights.npy")
        np.save(str(model_path) + f"/{self.file_prefix}SingularValues.npy", self.sigma_vec)


    def EvaluateNetwork(self):

        model_name, _ = self.data_base.get_hashed_file_name_for_table("Neural_Network", self.mu_train)
        model_path = pathlib.Path(self.data_base.database_root_directory / 'saved_nn_models' / model_name)

        # Load using the prefix
        with open(str(model_path) + f'/{self.file_prefix}train_config.json', "r") as config_file:
            model_properties = json.load(config_file)

        n_inf = model_properties['modes'][0]
        n_sup = model_properties['modes'][1]
        layers_size = model_properties['layers_size']

        network = self._DefineNetwork(n_inf, n_sup, layers_size)
        network.summary()

        # Load weights using the prefix
        network.load_weights(str(model_path) + f'/{self.file_prefix}model.weights.h5')

        S_val, Q_inf_val, Q_sup_val, phisig_inf, phisig_sup = self._GetEvaluationData(model_properties)

        print('RECONSTRUCTION RESULTS, VALIDATION DATASET:')
        print(' - Relative Frobenius error:')

        S_recons_val = phisig_sup@network(Q_inf_val).numpy().T+phisig_inf@Q_inf_val.T
        err_rel_recons = np.linalg.norm(S_recons_val-S_val)/np.linalg.norm(S_val)
        print('     ANN-PROM: ', err_rel_recons)

        S_pod_sup_recons_val = phisig_sup@Q_sup_val.T+phisig_inf@Q_inf_val.T
        print('     POD Sup: ', np.linalg.norm(S_pod_sup_recons_val-S_val)/np.linalg.norm(S_val))

        S_pod_inf_recons_val = phisig_inf@Q_inf_val.T
        print('     POD Inf: ', np.linalg.norm(S_pod_inf_recons_val-S_val)/np.linalg.norm(S_val))


        print(' - Relative geometric mean of the L2 error of each snapshot:')
        sample_l2_err_list=[]
        for i in range(S_recons_val.shape[1]):
            sample_l2_err_list.append(np.linalg.norm(S_recons_val[:,i]-S_val[:,i])/np.linalg.norm(S_val[:,i]))
        print('     ANN-PROM: ', np.linalg.norm(np.exp(np.mean(np.log(sample_l2_err_list)))))

        sample_l2_err_list_pod_sup=[]
        for i in range(S_pod_sup_recons_val.shape[1]):
            sample_l2_err_list_pod_sup.append(np.linalg.norm(S_pod_sup_recons_val[:,i]-S_val[:,i])/np.linalg.norm(S_val[:,i]))
        print('     POD Sup: ', np.linalg.norm(np.exp(np.mean(np.log(sample_l2_err_list_pod_sup)))))

        sample_l2_err_list_pod_inf=[]
        for i in range(S_pod_inf_recons_val.shape[1]):
            sample_l2_err_list_pod_inf.append(np.linalg.norm(S_pod_inf_recons_val[:,i]-S_val[:,i])/np.linalg.norm(S_val[:,i]))
        print('     POD Inf: ', np.linalg.norm(np.exp(np.mean(np.log(sample_l2_err_list_pod_inf)))))

        return err_rel_recons