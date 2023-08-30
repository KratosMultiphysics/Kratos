import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
import numpy as np


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return KerasOptimizerCompilerProcess(settings["parameters"])


class KerasOptimizerCompilerProcess(NeuralNetworkProcess):
    """
    This class compiles the optimizer among those available in the Keras API.

    This process is not tunable. Use custom_optimizer_compiler instead for tuning.

    """
    def __init__(self, parameters):
        super().__init__()
        if parameters.Has("optimizer"): 
            self.optimizer = parameters["optimizer"].GetString()
        else:
            raise Exception("No optimizer specified")

    def Compile(self, loss_function, optimizer):
        """ Process for compiling a network. """
        optimizer = self.optimizer
        return [loss_function, optimizer]