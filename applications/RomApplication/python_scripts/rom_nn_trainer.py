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

    def __init__(self, general_rom_manager_parameters):

        self.general_rom_manager_parameters = general_rom_manager_parameters
        self.nn_parameters = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]
        self.nn_parameters.RecursivelyValidateAndAssignDefaults(self._GetDefaultNeuralNetworkParameters())

    @classmethod
    def _GetDefaultNeuralNetworkParameters(self):
        nn_training_parameters = KratosMultiphysics.Parameters("""{
            "saved_models_root_path": "rom_data/saved_nn_models/",
            "training":{
                "modes":[5,50],
                "layers_size":[200,200],
                "batch_size":2,
                "epochs":800,
                "lr_strategy":{
                    "scheduler": "sgdr",
                    "base_lr": 0.001,
                    "additional_params": [1e-4, 10, 400]
                },
                "database":{
                    "training_set": "rom_data/SnapshotsMatrices/fom_snapshots.npy",
                    "validation_set": "rom_data/SnapshotsMatrices/fom_snapshots_val.npy",
                    "phi_matrix": "rom_data/RightBasisMatrix.npy",
                    "sigma_vector": "rom_data/SingularValuesVector.npy"
                },
                "use_automatic_name": false,
                "custom_name": "test_neural_network"                         
            },
            "online":{
                "model_name": "test_neural_network"
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

    def _GetTrainingData(self, n_inf, n_sup, database_settings):
        
        S_train = np.load(database_settings['training_set'].GetString())
        S_val = np.load(database_settings['validation_set'].GetString())

        phi = np.load(database_settings['phi_matrix'].GetString())
        sigma_vec = np.load(database_settings['sigma_vector'].GetString())/np.sqrt(S_train.shape[1])

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

        S_train = np.load(model_properties['database']['training_set'])
        S_val = np.load(model_properties['database']['validation_set'])

        n_inf = model_properties['modes'][0]
        n_sup = model_properties['modes'][1]
        phi = np.load(model_properties['database']['phi_matrix'])
        sigma_vec = np.load(model_properties['database']['sigma_vector'])/np.sqrt(S_train.shape[1])

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
            layer_out=layers.Dense(int(layer_size), 'elu', use_bias=False, kernel_initializer="he_normal", dtype=tf.float64)(layer_out)
        output_layer=layers.Dense(n_sup-n_inf, 'linear', use_bias=False, kernel_initializer="he_normal", dtype=tf.float64)(layer_out)

        network=Model(input_layer, output_layer)
        return network

    def TrainNetwork(self, seed=None):
        
        if seed is not None:
            keras_set_random_seed(seed)


        nn_training_parameters = self.nn_parameters['training']
        
        n_inf = int(nn_training_parameters['modes'].GetVector()[0])
        n_sup = int(nn_training_parameters['modes'].GetVector()[1])
        layers_size = nn_training_parameters['layers_size'].GetVector()
        batch_size = nn_training_parameters['batch_size'].GetInt()
        lr_scheme = nn_training_parameters['lr_strategy']['scheduler'].GetString()
        base_lr = nn_training_parameters['lr_strategy']['base_lr'].GetDouble()
        lr_additional_params = nn_training_parameters['lr_strategy']['additional_params'].GetVector()
        epochs = nn_training_parameters['epochs'].GetInt()
        
        if nn_training_parameters['use_automatic_name'].GetBool():
            model_name='NN_model_'+str(n_inf)+'.'+str(n_sup)+'_'+str(layers_size)+'_lr'+lr_scheme+'.'+str(base_lr)+'_batchsize'+str(batch_size)
        else:
            model_name=nn_training_parameters['custom_name'].GetString()

        model_path=pathlib.Path(self.nn_parameters['saved_models_root_path'].GetString()+model_name)
        model_path.mkdir(parents=True, exist_ok=False)

        database_settings = nn_training_parameters['database']
        Q_inf_train, Q_inf_val, Q_sup_train, Q_sup_val, phisig_norm_matrix, rescaling_factor = self._GetTrainingData(n_inf, n_sup, database_settings)

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
            },
            "database":{
                "training_set": nn_training_parameters["database"]["training_set"].GetString(),
                "validation_set": nn_training_parameters["database"]["validation_set"].GetString(),
                "phi_matrix": nn_training_parameters["database"]["phi_matrix"].GetString(),
                "sigma_vector": nn_training_parameters["database"]["sigma_vector"].GetString()
            }
        }

        with open(str(model_path)+"/train_config.json", "w") as ae_config_json_file:
            json.dump(training_parameters_dict, ae_config_json_file)

        network.save_weights(str(model_path)+"/model.weights.h5")
        with open(str(model_path)+"/history.json", "w") as history_file:
            json.dump(str(history.history), history_file)

        return model_name

    def EvaluateNetwork(self, model_name):

        model_path=pathlib.Path(self.nn_parameters['saved_models_root_path'].GetString()+model_name)

        with open(str(model_path)+'/train_config.json', "r") as config_file:
            model_properties = json.load(config_file)
        
        n_inf = model_properties['modes'][0]
        n_sup = model_properties['modes'][1]
        layers_size = model_properties['layers_size']
                                                  
        network = self._DefineNetwork(n_inf, n_sup, layers_size)
        network.summary()

        network.load_weights(str(model_path)+'/model.weights.h5')

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