import numpy as np
import os
import json

import tensorflow as tf
from keras.models import Model
from keras import layers
from keras.optimizers import AdamW, SGD
from keras.callbacks import LearningRateScheduler

from matplotlib import pyplot as plt


class Rom_NN_trainer(object):

    def __init__(self, general_rom_manager_parameters):

        self.general_rom_manager_parameters = general_rom_manager_parameters
        self.nn_training_parameters = self._SetNNTrainingParameters()

    def _SetNNTrainingParameters(self):
        defaults = self._GetDefaultNNTrainingParameters()

        nn_params = self.general_rom_manager_parameters["NN"]

        keys_to_copy = [
             'modes',
             'layers_size',
             'batch_size',
             'epochs',
             'lr_strategy',
             'saved_models_root_path',
             'database'
        ]

        for key in keys_to_copy:
            if key in nn_params.keys():
                defaults["Parameters"][key] = nn_params[key]

        return defaults

    def _GetTrainingData(self):
        
        database_settings = self.nn_training_parameters['database']

        S_train = np.load(database_settings['training_set'])
        S_val = np.load(database_settings['validation_set'])

        n_inf = self.nn_training_parameters['modes'][0]
        n_sup = self.nn_training_parameters['modes'][1]
        phi = np.load(database_settings['phi_matrix'])
        sigma_vec = np.load(database_settings['sigma_vector'])/np.sqrt(S_train.shape[1])

        phisig_inv_inf = np.linalg.inv(np.diag(sigma_vec[:n_inf]))@phi[:,:n_inf].T
        phisig_inv_sup = np.linalg.inv(np.diag(sigma_vec[:n_sup]))@phi[:,:n_sup].T
        phisig_inf = phi[:,:n_inf]@np.diag(sigma_vec[:n_inf])
        phisig_sup = phi[:,:n_sup]@np.diag(sigma_vec[:n_sup])

        Q_inf_train = (phisig_inv_inf@S_train).T
        Q_inf_val = (phisig_inv_inf@S_val).T
        Q_sup_train = (phisig_inv_sup@S_train).T
        Q_sup_val = (phisig_inv_sup@S_val).T

        phisig_norm_matrix = phisig_sup.T @ phisig_sup

        rescaling_factor = np.mean(np.square((phisig_inf@Q_inf_train.T)-S_train)) 
        rescaling_factor *= S_train.shape[1]/Q_sup_train.shape[1]

        Q_inf_train, Q_inf_val, Q_sup_train, Q_sup_val, phisig_norm_matrix, rescaling_factor

    def _GetEvaluationData(model_properties):

        S_train = np.load(model_properties['training_set'])
        S_val = np.load(model_properties['validation_set'])

        n_inf = model_properties['modes'][0]
        n_sup = model_properties['modes'][1]
        phi = np.load(model_properties['phi_matrix'])
        sigma_vec = np.load(model_properties['sigma_vector'])/np.sqrt(S_train.shape[1])

        phisig_inv_inf = np.linalg.inv(np.diag(sigma_vec[:n_inf]))@phi[:,:n_inf].T
        phisig_inv_sup = np.linalg.inv(np.diag(sigma_vec[:n_sup]))@phi[:,:n_sup].T
        phisig_inf = phi[:,:n_inf]@np.diag(sigma_vec[:n_inf])
        phisig_sup = phi[:,:n_sup]@np.diag(sigma_vec[:n_sup])

        Q_inf_val = (phisig_inv_inf@S_val).T
        Q_sup_val = (phisig_inv_sup@S_val).T

        return S_val, Q_inf_val, Q_sup_val, phisig_inf, phisig_sup
    
    def select_scheduler(self, strategy_name):

        def lr_const_scheduler(epoch, lr):
            new_lr= 0.001
            return new_lr
        
        def lr_steps_scheduler(epoch, lr):
            if epoch<=400:
                lr=0.001
            else:
                lr=0.0001
            return lr
        
        def lr_sgdr_scheduler(epoch, lr):
            """
            Blog post: jeremyjordan.me/nn-learning-rate
            Original paper: http://arxiv.org/abs/1608.03983
            """
            if epoch==0:
                new_lr=0.001
            else:
                max_lr=0.001
                min_lr=1e-4
                cycle_length=400
                scale_factor=10
                cycle=np.floor(epoch/cycle_length)
                x=np.abs(epoch/cycle_length-cycle)
                new_lr = min_lr+(max_lr-min_lr)*0.5*(1+np.cos(x*np.pi))/scale_factor**cycle
            return new_lr
        
        schedulers_dict={"const": lr_const_scheduler, "steps": lr_steps_scheduler, "sgdr": lr_sgdr_scheduler}
        
        return schedulers_dict[strategy_name]
    
    def define_network(self, n_inf, n_sup, layers_size):
        input_layer=layers.Input((n_inf,), dtype=tf.float64)
        layer_out=input_layer
        for layer_size in layers_size:
            layer_out=layers.Dense(layer_size, 'elu', use_bias=False, kernel_initializer="he_normal", dtype=tf.float64)(layer_out)
        output_layer=layers.Dense(n_sup-n_inf, 'linear', use_bias=False, kernel_initializer="he_normal", dtype=tf.float64)(layer_out)

        network=Model(input_layer, output_layer)
        return network

    def train_network(self):
        
        n_inf = self.nn_training_parameters['modes'][0]
        n_sup = self.nn_training_parameters['modes'][1]
        layers_size = self.nn_training_parameters['layers_size']
        batch_size = self.nn_training_parameters['batch_size']
        lr_strategy = self.nn_training_parameters['lr_strategy']
        epochs = self.nn_training_parameters['epochs']
        
        model_name='NN_model_'+str(n_inf)+'.'+str(n_sup)+'_'+str(layers_size)+'_lr'+lr_strategy+'_batchsize'+str(batch_size)
        model_path=self.nn_training_parameters['saved_models_root_path']+model_name+'/'

        os.makedirs(model_path, exist_ok=False)

        Q_inf_train, Q_inf_val, Q_sup_train, Q_sup_val, phisig_norm_matrix, rescaling_factor = self._GetTrainingData()

        network = self.define_network(n_inf, n_sup, layers_size)

        def scaled_phinorm_mse_loss(y_true, y_pred):
            y_diff=y_true-y_pred
            mse = tf.math.reduce_mean(tf.math.multiply(y_diff,tf.transpose(tf.matmul(phisig_norm_matrix,tf.transpose(y_diff)))))
            mse /= rescaling_factor
            return mse

        network.compile(AdamW(epsilon=1e-17), loss=scaled_phinorm_mse_loss, run_eagerly=False)
        network.summary()

        callbacks = [LearningRateScheduler(self.select_scheduler(lr_strategy), verbose=0)]

        history = network.fit(Q_inf_train, Q_sup_train, batch_size=batch_size, epochs=epochs, validation_data=(Q_inf_val,Q_sup_val), shuffle=True, validation_batch_size=1, callbacks=callbacks)
        
        with open(model_path+"train_config.json", "w") as ae_config_json_file:
            json.dump(self.nn_training_parameters, ae_config_json_file)

        network.save_weights(model_path+"model_weights.h5")
        with open(model_path+"history.json", "w") as history_file:
            json.dump(str(history.history), history_file)

        return model_name

    def evaluate_network(self, model_name):

        model_path=self.nn_training_parameters['saved_models_root_path']+model_name+'/'

        with open(model_path+'train_config.json', "r") as config_file:
            model_properties = json.load(config_file)
        
        n_inf = model_properties['modes'][0]
        n_sup = model_properties['modes'][1]
        layers_size = model_properties['layers_size']
                                                  
        network = self.define_network(n_inf, n_sup, layers_size)
        network.summary()

        network.load_weights(model_path+'model_weights.h5')

        S_val, Q_inf_val, Q_sup_val, phisig_inf, phisig_sup = self._GetEvaluationData(model_properties)

        S_recons_val = (phisig_sup@network(Q_inf_val).numpy().T).T+(phisig_inf@Q_inf_val.T).T
        print('ANN-PROM Validation Reconstruction error (Rel. Frob): ', np.linalg.norm(S_recons_val-S_val)/np.linalg.norm(S_val))

        S_pod_sup_recons_val = (phisig_sup@Q_sup_val.T).T+(phisig_inf@Q_inf_val.T).T
        print('POD Sup Validation Reconstruction error (Rel. Frob): ', np.linalg.norm(S_pod_sup_recons_val-S_val)/np.linalg.norm(S_val))
        S_pod_inf_recons_val = (phisig_inf@Q_inf_val.T).T
        print('POD Inf Validation Reconstruction error (Rel. Frob): ', np.linalg.norm(S_pod_inf_recons_val-S_val)/np.linalg.norm(S_val))

# if __name__=="__main__":

#     nn_training_config = {
#         'modes':[6,60],
#         'layers_size':[200,200],
#         'batch_size':16,
#         'epochs':200,
#         'lr_strategy':'sgdr',
#         'saved_models_root_path': 'rom_data/saved_nn_models/',
#         'database':{
#             'training_set': 'rom_data/SnapshotsMatrices/fom_snapshots.npy',
#             'validation_set': 'rom_data/SnapshotsMatrices/fom_snapshots_val.npy',
#             'phi_matrix': 'rom_data/RightBasisMatrix.npy',
#             'sigma_vector': 'rom_data/SingularValuesVector.npy'
#         }
#     }