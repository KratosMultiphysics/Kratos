import KratosMultiphysics as KM
import tensorflow.keras as keras
from KratosMultiphysics.NeuralNetworkApplication.neural_network_training_process import NeuralNetworkTrainingProcess


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return GenerateSimpleModelProcess(settings["parameters"])

class GenerateSimpleModelProcess(NeuralNetworkTrainingProcess):
    """
    This class creates a process that generates a model of a simple fully connected layer neural network 
    based on the data provided for the trianing.
    """
    def __init__(self, parameters):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        settings -- Kratos parameters containing layer settings.
        """
        super().__init__(parameters)
        self.model_name = parameters["model_name"].GetString() 

    def ExecuteInitialize(self):
        self.model=self.__getmodel(self.test_input[0,:].shape,self.test_output[0,:].size)
        self.model.save(self.model_name)

    def __getmodel(self, input_shape, output_shape):
        """ Model definition """
        inputs = keras.Input(shape=input_shape)
        outputs = keras.layers.Dense(output_shape)(inputs)
        model = keras.Model(inputs, outputs)
        model.compile(optimizer="adam", loss="mean_squared_error")
        return model
