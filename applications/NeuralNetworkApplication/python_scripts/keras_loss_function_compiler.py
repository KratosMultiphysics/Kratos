import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
import numpy as np


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return KerasLossFunctionCompilerProcess(settings["parameters"])


class KerasLossFunctionCompilerProcess(NeuralNetworkProcess):
    """
    This class compiles the loss function among those available in the Keras API.

    This process is not tunable. Use custom_loss_function_compiler instead for tuning.

    """
    def __init__(self, parameters):
        super().__init__()
        if parameters.Has("loss_function"): 
            self.loss_function = parameters["loss_function"].GetString()
        else:
            raise Exception("No loss function specified")
            
    def Compile(self, loss_function, optimizer):
        """ Process for compiling a network. """
        loss_function = self.loss_function
        return [loss_function, optimizer]