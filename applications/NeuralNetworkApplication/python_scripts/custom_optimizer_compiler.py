import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
import numpy as np
import tensorflow.keras.optimizers
from importlib import import_module


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return CustomOptimizerCompilerProcess(settings["parameters"])


class CustomOptimizerCompilerProcess(NeuralNetworkProcess):
    """
    This class compiles the loss function among those available in the Keras API.

    """
    def __init__(self, parameters):
        super().__init__()

        default_settings = KM.Parameters("""{
            "help"                     : "This process compiles a custom optimizer.",
            "optimizer"                : "",
            "module"                   : "keras",
            "learning_rate"            : 0.001
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        if parameters.Has("optimizer"): 
            self.optimizer_name = parameters["optimizer"].GetString()
        else:
            raise Exception("No optimizer specified")
        self.module = parameters["module"].GetString()
        self.learning_rate = parameters["learning_rate"].GetDouble()
        optimizer_module = import_module(self.module)
        if self.module == "keras":
            self.optimizer = getattr(tensorflow.keras.optimizers,self.optimizer_name)(learning_rate = self.learning_rate)
        else:
            self.optimizer = getattr(optimizer_module,optimizer)
            
    def Compile(self, loss_function, optimizer):
        """ Process for compiling a network. """
        optimizer = self.optimizer
        return [loss_function, optimizer]