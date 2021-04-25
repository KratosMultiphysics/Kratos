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

    def ExecuteTraining(self):
        self.reconstructed_model = keras.models.load_model(self.model_name)
        self.reconstructed_model.compile(loss="mean_squared_error", optimizer=keras.optimizers.Adam(learning_rate=0.001, decay = 1e-3 / 200))
        self.reconstructed_model.fit(self.test_input, self.test_output, epochs = self.epochs, shuffle=False)
        self.reconstructed_model.save(self.model_name)

