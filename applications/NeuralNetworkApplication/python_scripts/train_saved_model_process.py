import tensorflow.keras as keras
import numpy as np
import KratosMultiphysics
import KratosMultiphysics.NeuralNetworkApplication as NeuralNetworkApplication
from KratosMultiphysics.NeuralNetworkApplication.neural_network_training_process import NeuralNetworkTrainingProcess

def Factory(settings):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return TrainerFromSavedModel(settings["parameters"])

class TrainerFromSavedModel(NeuralNetworkTrainingProcess):
    """ 
    This class creates a process that trains a neural network located at a saved model.
    """
    def __init__(self, parameters):
        super().__init__(parameters)
        self.model_name = parameters["model"].GetString()
        self.epochs = parameters["epochs"].GetInt()
        self.training_log = parameters["training_log"].GetBool()
        self.shuffle = parameters["shuffle"].GetBool()
        if self.training_log:
            self.training_log_file = parameters["training_log_file"].GetString()

    def ExecuteTraining(self, loss_function, optimizer, callbacks_list, metrics_list):
        """ Processes to act directly during the training step. """
        self.reconstructed_model = keras.models.load_model(self.model_name)
        self.reconstructed_model.compile(loss=loss_function, optimizer=optimizer, metrics = metrics_list)
        if self.validation:
            history = self.reconstructed_model.fit(self.test_input.ExportAsArray(), self.test_output.ExportAsArray(), epochs = self.epochs, 
            validation_data = (self.val_input, self.val_output), shuffle=self.shuffle, callbacks = callbacks_list)
        else:
            history = self.reconstructed_model.fit(self.test_input.ExportAsArray(), self.test_output.ExportAsArray(), epochs = self.epochs, shuffle=self.shuffle,
            validation_split = self.validation_split, callbacks = callbacks_list) 
        self.reconstructed_model.save(self.model_name)
        if self.training_log:
            with open(self.training_log_file,'w') as f:
                print(history.history, file = f)

