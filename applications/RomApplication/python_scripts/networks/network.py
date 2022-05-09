import os
import abc

import numpy as np

import keras
import tensorflow as tf

from keras import layers

class Network(abc.ABC):

    def __init__(self):
        """
        Initis the network class

        model_name:     name of the directory in which the model will be saved.
        valid:          percentage of the input which is valid. (First steps are typically)
        """

        self.model_name = "./saved_models/emtpy_network"
        self.valid = 1.0

    @abc.abstractmethod
    def define_network(self, input_data, custom_loss):
        """ Define the network configuration """

    @abc.abstractmethod
    def calculate_gradients(self, input_data, encoder, decoder, lossf, trainable_variables):
        """ Compute the gradients of the input network list for the given input variables """

    @abc.abstractmethod
    def compute_full_gradient(self, network, all_gradients):
        """ Operates all gradients fir the given input variables """

    def check_gradient(self, SEncoded, SDecoded):
        """ Checks the correctness of a given gradient """
        print("Method not implemented.")

    def calculate_data_limits(self, input_data):
        self.data_min = np.min(input_data)
        self.data_max = np.max(input_data)

        self.data_min = 0
        # self.data_max = 1

    def normalize_data(self, input_data):
        return (input_data - self.data_min) / (self.data_max - self.data_min)

    def denormalize_data(self, input_data):
        return input_data * (self.data_max - self.data_min) + self.data_min

    # def normalize_data(self, input_data):
    #     return (input_data - self.data_min) / (self.data_max - self.data_min)

    # def denormalize_data(self, input_data):
    #     return input_data * (self.data_max - self.data_min) + self.data_min

    def normalize_data(self, input_data):
        return input_data - self.data_min

    def denormalize_data(self, input_data):
        return input_data + self.data_min

    def expand_dataset(self, data, times):
        noise_data = np.tile(data, (times, 1))
        noise_data = noise_data + np.random.normal(0, 0.1, noise_data.shape)
        return np.concatenate((data, noise_data), axis=0)

    def prepare_data(self, input_data, num_files):
        data = np.transpose(input_data)

        print("Data shape:", data.shape)
        print("NOR - Min:", np.min(data), "Max:", np.max(data))

        # Select some of the snapshots to train and some others to validate
        train_cut = len(data) / num_files

        train_pre = [data[i] for i in range(0, data.shape[0]) if (i % train_cut) <  (train_cut * self.valid)]
        valid_pre = [data[i] for i in range(0, data.shape[0]) if (i % train_cut) >= (train_cut * self.valid)]

        train_samples = np.array(train_pre)
        valid_samples = np.array(valid_pre)

        train_dataset = np.asarray(train_samples)
        valid_dataset = np.asarray(valid_samples)

        return train_samples, valid_samples

    def train_network(self, model, input_data, grad_data, num_files):
        # train_dataset, valid_dataset = self.prepare_data(input_data, num_files)

        # Shuffle the snapshots to prevent batches from the same clusters
        # np.random.shuffle(train_dataset)
        # np.random.shuffle(valid_dataset)

        # Train the model
        model.grads = grad_data
        model.fit(
            input_data.T, grad_data.T,
            epochs=10,
            batch_size=1,
            # shuffle=True,
            # validation_data=(valid_dataset, valid_dataset),
        )

    def predict_vector(self, network, input_vector):

        tmp = input_vector.reshape(1,len(input_vector))
        predicted_vector = network.predict(tmp)

        return np.asarray(predicted_vector.T)

    def get_gradients(self, model, models_variables, input_data, output_data):
        ''' Calculates the derivative of a neural network model respect of the input '''
        in_tf_var = tf.Variable(input_data)
        with tf.GradientTape(persistent=True) as tape:
            tape.watch(in_tf_var)
            output = model(in_tf_var, training=False)
        auto_grad = tape.batch_jacobian(output, in_tf_var, unconnected_gradients=tf.UnconnectedGradients.ZERO, experimental_use_pfor=False) 

        # Compute gradients
        return auto_grad[0].numpy()