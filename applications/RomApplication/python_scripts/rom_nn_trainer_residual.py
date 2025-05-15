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

from KratosMultiphysics.RomApplication.rom_nn_interface_structural import StructuralMechanics_NN_Interface
from KratosMultiphysics.RomApplication.rom_nn_prepost_processor import PODANN_prepost_processor

tf.keras.backend.set_floatx('float64')


class R_Only_Strategy_KerasModel(Model):

    def __init__(self, prepost_processor, structural_nn_interface, *args, **kwargs):
        super(R_Only_Strategy_KerasModel,self).__init__(*args,**kwargs)

        self.prepost_processor = prepost_processor
        self.structural_nn_interface = structural_nn_interface

        self.get_v_loss_r = self.structural_nn_interface.get_v_loss_rdiff_batch
        self.get_err_r = self.structural_nn_interface.get_err_rdiff_batch

        self.run_eagerly = False
        
        # self.rescaling_factor_r=1.0
        # self.rescaling_factor_x=1.0

        self.loss_x_tracker = metrics.Mean(name="loss_x")
        self.loss_r_tracker = metrics.Mean(name="loss_r")
    
    @tf.function
    def r_loss_scale(self, err_r, v_loss_r=None):
        loss_r = tf.math.reduce_mean(tf.math.square(err_r), axis=1)
        return loss_r, v_loss_r
    
    @tf.function
    def get_gradients(self, trainable_vars, input_batch, v_loss_r_batch):

        v_loss = 2*v_loss_r_batch/self.rescaling_factor_r

        with tf.GradientTape(persistent=False) as tape_d:
            tape_d.watch(trainable_vars)
            x_pred = self(input_batch, training=True)
            x_pred_denorm = self.prepost_processor.postprocess_output_data_tf(x_pred, (input_batch,None))
            v_u_dotprod = tf.math.reduce_mean(tf.math.multiply(v_loss, x_pred_denorm), axis=1)
            v_u_dotprod_mean = tf.math.reduce_mean(v_u_dotprod)
        grad_loss=tape_d.gradient(v_u_dotprod_mean, trainable_vars)

        return grad_loss

    def update_rescaling_factors(self, S_true, R_true):

        S_recons_aux1 = self.prepost_processor.preprocess_nn_output_data(S_true)
        S_recons_aux2, _ =self.prepost_processor.preprocess_input_data(S_true)
        S_recons = self.prepost_processor.postprocess_output_data(np.zeros(S_recons_aux1.shape), (S_recons_aux2, None))
        
        err_r_batch, _ = self.get_v_loss_r(tf.constant(S_recons), tf.constant(R_true))
        rescaling_factor_r = np.mean(np.square(err_r_batch.numpy()))
        rescaling_factor_x = np.mean(np.square(S_recons-S_true))

        self.rescaling_factor_r = rescaling_factor_r
        self.rescaling_factor_x = rescaling_factor_x

        print('Updated gradient rescaling factors. r: ', self.rescaling_factor_r)
        print('Updated rescaling factor X: ', self.rescaling_factor_x)

    def train_step(self,data):

        input_batch, (target_snapshot_batch,target_aux_batch) = data # target_aux is the residual
        trainable_vars = self.trainable_variables

        x_pred_batch = self(input_batch, training=True)
        x_pred_denorm_batch = self.prepost_processor.postprocess_output_data_tf(x_pred_batch,(input_batch,None))
        err_x_batch = target_snapshot_batch - x_pred_denorm_batch
        loss_x_batch = tf.math.reduce_mean(tf.math.square(err_x_batch), axis=1)

        err_r_batch, v_loss_r_batch = self.get_v_loss_r(x_pred_denorm_batch,target_aux_batch)
        loss_r_batch, v_loss_r_batch = self.r_loss_scale(err_r_batch, v_loss_r_batch)

        total_loss_x=tf.math.reduce_mean(loss_x_batch)/self.rescaling_factor_x
        total_loss_r=tf.math.reduce_mean(loss_r_batch)/self.rescaling_factor_r

        grad_loss = self.get_gradients(trainable_vars, input_batch, v_loss_r_batch)

        self.optimizer.apply_gradients(zip(grad_loss, trainable_vars))

        # Compute our own metrics
        self.loss_x_tracker.update_state(total_loss_x)
        self.loss_r_tracker.update_state(total_loss_r)
        
        return {"loss_x": self.loss_x_tracker.result(), "loss_r": self.loss_r_tracker.result()}

    def test_step(self, data):
        input_batch, (target_snapshot_batch,target_aux_batch) = data

        x_pred_batch = self(input_batch, training=False)
        x_pred_denorm_batch = self.prepost_processor.postprocess_output_data_tf(x_pred_batch,(input_batch,None))
        err_x_batch = target_snapshot_batch - x_pred_denorm_batch
        loss_x_batch = tf.math.reduce_mean(tf.math.square(err_x_batch), axis=1)

        err_r_batch = self.get_err_r(x_pred_denorm_batch, target_aux_batch)
        loss_r_batch, _ = self.r_loss_scale(err_r_batch)
        
        total_loss_r=tf.math.reduce_mean(loss_r_batch)/self.rescaling_factor_r
        total_loss_x=tf.math.reduce_mean(loss_x_batch)/self.rescaling_factor_x

        # Compute our own metrics
        self.loss_x_tracker.update_state(total_loss_x)
        self.loss_r_tracker.update_state(total_loss_r)
        return {"loss_x": self.loss_x_tracker.result(), "loss_r": self.loss_r_tracker.result()}

    @property
    def metrics(self):
        # We list our `Metric` objects here so that `reset_states()` can be
        # called automatically at the start of each epoch
        # or at the start of `evaluate()`.
        # If you don't implement this property, you have to call
        # `reset_states()` yourself at the time of your choosing.
        return [self.loss_x_tracker, self.loss_r_tracker]
    

class RomNeuralNetworkTrainerResidual(object):

    def __init__(self, general_rom_manager_parameters, mu_train, mu_validation, data_base, kratos_simulation):

        self.general_rom_manager_parameters = general_rom_manager_parameters
        self.nn_parameters = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]
        self.mu_train = mu_train
        self.mu_validation = mu_validation
        self.data_base = data_base
        self.structural_nn_interface = StructuralMechanics_NN_Interface(kratos_simulation)

    def _CheckNumberOfModes(self,n_inf,n_sup,n_max):
        if n_inf >= n_max:
            err_msg = f'Specified number of inferior modes ({n_inf}) is higher than or equal to the available ones from the Phi matrix ({n_max}).'
            raise Exception(err_msg)
        elif n_sup > n_max:
            err_msg = f'Specified number of superior modes ({n_sup}) is higher than the available ones from the Phi matrix ({n_max}).'
            raise Exception(err_msg)

    def _GetTrainingData(self, n_inf, n_sup, kratos_simulation):

        S_train = self.data_base.get_snapshots_matrix_from_database(self.mu_train, table_name=f'FOM').T
        S_val = self.data_base.get_snapshots_matrix_from_database(self.mu_validation, table_name=f'FOM').T
        print('Shape S_train: ', S_train.shape)
        print('Shape S_val: ', S_val.shape)

        R_train =  kratos_simulation.get_r_batch_(S_train)
        R_val =  kratos_simulation.get_r_batch_(S_val)
        print('Shape R_train: ', R_train.shape)
        print('Shape R_val: ', R_val.shape)

        _, hash_basis = self.data_base.check_if_in_database("RightBasis", self.mu_train)
        phi = self.data_base.get_single_numpy_from_database(hash_basis)
        _, hash_sigma = self.data_base.check_if_in_database("SingularValues_Solution", self.mu_train)
        sigma_vec =  self.data_base.get_single_numpy_from_database(hash_sigma)
        print('Shape phi: ', phi.shape)
        print('Shape sigma_vec: ', sigma_vec.shape)

        self._CheckNumberOfModes(n_inf,n_sup,sigma_vec.shape[0])

        return S_train, S_val, R_train, R_val, phi, sigma_vec

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

    def _DefineNetwork(self, n_inf, n_sup, layers_size, prepost_processor, structural_nn_interface):
        input_layer=layers.Input((n_inf,), dtype=tf.float64)
        layer_out=input_layer
        for layer_size in layers_size:
            layer_out=layers.Dense(int(layer_size), 'elu', use_bias=False, kernel_initializer="he_normal", dtype=tf.float64)(layer_out)
        output_layer=layers.Dense(n_sup-n_inf, 'linear', use_bias=False, kernel_initializer="he_normal", dtype=tf.float64)(layer_out)

        network=R_Only_Strategy_KerasModel(prepost_processor, structural_nn_interface, input_layer, output_layer)
        return network

    def _SaveWeightsKratosFormat(self, network, weights_path):
        layers=[]
        for layer in network.trainable_variables:
            layers.append(layer.numpy())

        np.save(weights_path, np.array(layers, dtype=object), allow_pickle=True)

    def TrainNetwork(self):

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

        model_name, _ = self.data_base.get_hashed_file_name_for_table("Neural_Network_Residual", self.mu_train)
        model_path=pathlib.Path(self.data_base.database_root_directory / 'saved_nn_models_residual' / model_name)
        model_path.mkdir(parents=True, exist_ok=False)
        print('New_model_path: ', model_path)

        S_train, S_val, R_train, R_val, phi, sigma_vec = self._GetTrainingData(n_inf, n_sup, self.structural_nn_interface)

        prepost_processor=PODANN_prepost_processor(phi, sigma_vec, S_train, n_inf, n_sup)

        network = self._DefineNetwork(n_inf, n_sup, layers_size, prepost_processor, self.structural_nn_interface)

        network.compile(AdamW(epsilon=1e-17), run_eagerly=False)
        network.summary()

        callbacks = [LearningRateScheduler(self._SelectScheduler(lr_scheme, base_lr, lr_additional_params), verbose=0),
                     tf.keras.callbacks.ModelCheckpoint(str(model_path)+"/model.weights.h5",save_weights_only=True,save_best_only=True,monitor="val_loss_r",mode="min")]
        
        pretrained_model_path = pathlib.Path(nn_training_parameters['online']['custom_model_path'].GetString())
        # pretrained_model_path = pathlib.Path(self.data_base.database_root_directory / 'saved_nn_models' / pretrained_model_name)
        network.load_weights(pretrained_model_path / 'model.weights.h5')
        print('Using pre-trained weights from: ', pretrained_model_path)
        
        input_data, target_data, val_input, val_target = prepost_processor.get_training_data(S_train, S_val, R_train, R_val)
        
        network.update_rescaling_factors(target_data[0],target_data[1])
        history = network.fit(input_data, target_data, batch_size=batch_size, epochs=epochs, validation_data=(val_input,val_target), shuffle=True, callbacks=callbacks)

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

        with open(str(model_path)+"/train_config.json", "w") as ae_config_json_file:
            json.dump(training_parameters_dict, ae_config_json_file)

        network.load_weights(str(model_path)+"/model.weights.h5")
        self._SaveWeightsKratosFormat(network, str(model_path)+"/model_weights.npy")