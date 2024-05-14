import numpy as np
import pathlib
import json

import tensorflow as tf
from keras.models import Model
from keras import layers
from keras.optimizers import AdamW
from keras.callbacks import LearningRateScheduler
from keras.utils import set_random_seed as keras_set_random_seed

import KratosMultiphysics


class RomNeuralNetworkTrainer(object):

    def __init__(self, general_rom_manager_parameters, mu_train, mu_validation, data_base):

        self.general_rom_manager_parameters = general_rom_manager_parameters
        self.nn_parameters = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]
        self.nn_parameters.RecursivelyValidateAndAssignDefaults(self._GetDefaultNeuralNetworkParameters())
        self.mu_train = mu_train
        self.mu_validation = mu_validation
        self.data_base = data_base

    @classmethod
    def _GetDefaultNeuralNetworkParameters(self):
        nn_training_parameters = KratosMultiphysics.Parameters("""{
            "modes":[5,50],
            "layers_size":[200,200],
            "batch_size":2,
            "epochs":800,
            "lr_strategy":{
                "scheduler": "sgdr",
                "base_lr": 0.001,
                "additional_params": [1e-4, 10, 400]
            },
            "training":{
                "retrain_if_exists" : false,  // If false only one model will be trained for each the mu_train and NN hyperparameters combination
                "model_number" : 0     // this part of the parameters will be updated with the number of trained models that exist for the same mu_train and NN hyperparameters combination
            },
            "online":{
                "model_number": 0   // out of the models existing for the same parameters, this is the model that will be lauched
            }
         }""")
        return nn_training_parameters

    def _CheckNumberOfModes(self,n_inf,n_sup,n_max):
        if n_inf >= n_max:
            err_msg = f'Specified number of inferior modes ({n_inf}) is higher than or equal to the available ones from the Phi matrix ({n_max}).'
            raise Exception(err_msg)
        elif n_sup > n_max:
            err_msg = f'Specified number of superior modes ({n_sup}) is higher than the available ones from the Phi matrix ({n_max}).'
            raise Exception(err_msg)

    def _GetTrainingData(self, n_inf, n_sup):

        S_train = self.data_base.get_snapshots_matrix_from_database(self.mu_train, table_name=f'FOM_Fit')
        S_val = self.data_base.get_snapshots_matrix_from_database(self.mu_validation, table_name=f'FOM_Fit')

        _, hash_basis = self.data_base.check_if_in_database("RightBasis", self.mu_train)
        phi = self.data_base.get_single_numpy_from_database(hash_basis)
        _, hash_sigma = self.data_base.check_if_in_database("SingularValues_Solution", self.mu_train)
        sigma_vec =  self.data_base.get_single_numpy_from_database(hash_sigma)

        self._CheckNumberOfModes(n_inf,n_sup,sigma_vec.shape[0])

        phisig_inv_inf = np.linalg.inv(np.diag(sigma_vec[:n_inf]))@phi[:,:n_inf].T
        phisig_inv_sup = np.linalg.inv(np.diag(sigma_vec[n_inf:n_sup]))@phi[:,n_inf:n_sup].T
        phisig_inf = phi[:,:n_inf]@np.diag(sigma_vec[:n_inf])
        phisig_sup = phi[:,n_inf:n_sup]@np.diag(sigma_vec[n_inf:n_sup])

        Q_inf_train = (phisig_inv_inf@S_train).T
        Q_inf_val = (phisig_inv_inf@S_val).T
        Q_sup_train = (phisig_inv_sup@S_train).T
        Q_sup_val = (phisig_inv_sup@S_val).T

        phisig_norm_matrix = phisig_sup.T @ phisig_sup

        rescaling_factor = np.mean(np.square((phisig_inf@Q_inf_train.T)-S_train))
        rescaling_factor *= S_train.shape[0]/Q_sup_train.shape[1]

        return Q_inf_train, Q_inf_val, Q_sup_train, Q_sup_val, phisig_norm_matrix, rescaling_factor

    def _GetEvaluationData(self, model_properties):

        S_val = self.data_base.get_snapshots_matrix_from_database(self.mu_validation, table_name=f'FOM_Fit')

        n_inf = model_properties['modes'][0]
        n_sup = model_properties['modes'][1]

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

    def _DefineNetwork(self, n_inf, n_sup, layers_size):
        input_layer=layers.Input((n_inf,), dtype=tf.float64)
        layer_out=input_layer
        for layer_size in layers_size:
            layer_out=layers.Dense(layer_size, 'elu', use_bias=False, kernel_initializer="he_normal", dtype=tf.float64)(layer_out)
        output_layer=layers.Dense(n_sup-n_inf, 'linear', use_bias=False, kernel_initializer="he_normal", dtype=tf.float64)(layer_out)

        network=Model(input_layer, output_layer)
        return network

    def TrainNetwork(self, seed=None):

        if seed is not None:
            keras_set_random_seed(seed)


        nn_training_parameters = self.nn_parameters

        n_inf = int(nn_training_parameters['modes'].GetVector()[0])
        n_sup = int(nn_training_parameters['modes'].GetVector()[1])
        layers_size = nn_training_parameters['layers_size'].GetVector()
        batch_size = nn_training_parameters['batch_size'].GetInt()
        lr_scheme = nn_training_parameters['lr_strategy']['scheduler'].GetString()
        base_lr = nn_training_parameters['lr_strategy']['base_lr'].GetDouble()
        lr_additional_params = nn_training_parameters['lr_strategy']['additional_params'].GetVector()
        epochs = nn_training_parameters['epochs'].GetInt()

        model_name, _ = self.data_base.get_hashed_mu_for_table("Neural_Network", self.mu_train)
        model_path=pathlib.Path(self.data_base.database_root_directory / 'saved_nn_models' / model_name)
        model_path.mkdir(parents=True, exist_ok=False)


        Q_inf_train, Q_inf_val, Q_sup_train, Q_sup_val, phisig_norm_matrix, rescaling_factor = self._GetTrainingData(n_inf, n_sup)

        network = self._DefineNetwork(n_inf, n_sup, layers_size)

        def scaled_phinorm_mse_loss(y_true, y_pred):
            y_diff=y_true-y_pred
            mse = tf.math.reduce_mean(tf.math.multiply(y_diff,tf.transpose(tf.matmul(phisig_norm_matrix,tf.transpose(y_diff)))))
            mse /= rescaling_factor
            return mse

        network.compile(AdamW(epsilon=1e-17), loss=scaled_phinorm_mse_loss, run_eagerly=False)
        network.summary()

        callbacks = [LearningRateScheduler(self._SelectScheduler(lr_scheme, base_lr, lr_additional_params), verbose=0)]

        # history = network.fit(Q_inf_train, Q_sup_train, batch_size=batch_size, epochs=epochs, validation_data=(Q_inf_val,Q_sup_val), shuffle=False, validation_batch_size=1, callbacks=callbacks)
        history = network.fit(Q_inf_train, Q_sup_train, batch_size=batch_size, epochs=epochs, validation_data=(Q_inf_val,Q_sup_val), shuffle=True, callbacks=callbacks)


        training_parameters_dict = {
            "modes":[n_inf,n_sup],
            "layers_size":list(layers_size),
            "batch_size":batch_size,
            "epochs":epochs,
            "lr_strategy":{
                "scheduler": lr_scheme,
                "base_lr": base_lr,
                "additional_params": list(lr_additional_params)
            }
        }

        with open(str(model_path)+"/train_config.json", "w") as ae_config_json_file:
            json.dump(training_parameters_dict, ae_config_json_file)

        network.save_weights(str(model_path)+"/model_weights.h5")
        with open(str(model_path)+"/history.json", "w") as history_file:
            json.dump(str(history.history), history_file)

        return model_name

    def EvaluateNetwork(self, model_name):

        model_name, _ = self.data_base.get_hashed_mu_for_table("Neural_Network", self.mu_train)
        model_path=pathlib.Path(self.data_base.database_root_directory / 'saved_nn_models' / model_name)

        with open(str(model_path)+'/train_config.json', "r") as config_file:
            model_properties = json.load(config_file)

        n_inf = model_properties['modes'][0]
        n_sup = model_properties['modes'][1]
        layers_size = model_properties['layers_size']

        network = self._DefineNetwork(n_inf, n_sup, layers_size)
        network.summary()

        network.load_weights(str(model_path)+'/model_weights.h5')

        S_val, Q_inf_val, Q_sup_val, phisig_inf, phisig_sup = self._GetEvaluationData(model_properties)

        S_recons_val = phisig_sup@network(Q_inf_val).numpy().T+phisig_inf@Q_inf_val.T
        err_rel_recons = np.linalg.norm(S_recons_val-S_val)/np.linalg.norm(S_val)
        print('ANN-PROM Validation Reconstruction error (Rel. Frob): ', err_rel_recons)

        sample_l2_err_list=[]
        for i in range(S_recons_val.shape[1]):
            sample_l2_err_list.append(np.linalg.norm(S_recons_val[:,i]-S_val[:,i])/np.linalg.norm(S_val[:,i]))
        print('ANN-PROM Validation Reconstruction error (Geometric Rel. L2): ', np.linalg.norm(np.exp(np.mean(np.log(sample_l2_err_list)))))

        S_pod_sup_recons_val = phisig_sup@Q_sup_val.T+phisig_inf@Q_inf_val.T
        print('POD Sup Validation Reconstruction error (Rel. Frob): ', np.linalg.norm(S_pod_sup_recons_val-S_val)/np.linalg.norm(S_val))

        sample_l2_err_list=[]
        for i in range(S_pod_sup_recons_val.shape[1]):
            sample_l2_err_list.append(np.linalg.norm(S_pod_sup_recons_val[:,i]-S_val[:,i])/np.linalg.norm(S_val[:,i]))
        print('POD Sup Validation Reconstruction error (Geometric Rel. L2): ', np.linalg.norm(np.exp(np.mean(np.log(sample_l2_err_list)))))

        S_pod_inf_recons_val = phisig_inf@Q_inf_val.T
        print('POD Inf Validation Reconstruction error (Rel. Frob): ', np.linalg.norm(S_pod_inf_recons_val-S_val)/np.linalg.norm(S_val))

        sample_l2_err_list=[]
        for i in range(S_pod_inf_recons_val.shape[1]):
            sample_l2_err_list.append(np.linalg.norm(S_pod_inf_recons_val[:,i]-S_val[:,i])/np.linalg.norm(S_val[:,i]))
        print('POD Inf Validation Reconstruction error (Geometric Rel. L2): ', np.linalg.norm(np.exp(np.mean(np.log(sample_l2_err_list)))))

        return err_rel_recons