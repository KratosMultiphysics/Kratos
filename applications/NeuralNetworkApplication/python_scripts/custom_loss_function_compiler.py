import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
import numpy as np
import tensorflow.keras.losses
from importlib import import_module


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return CustomLossFunctionCompilerProcess(settings["parameters"])


class CustomLossFunctionCompilerProcess(NeuralNetworkProcess):
    """
    This class compiles ta custom loss function.

    """
    def __init__(self, parameters):
        super().__init__()

        default_settings = KM.Parameters("""{
            "help"                     : "This process compiles a custom loss function.",
            "loss_function"            : "",
            "module"                   : "keras",
            "reduction"                : "sum_over_batch_size"
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        if parameters.Has("loss_function"): 
            self.loss_function_name = parameters["loss_function"].GetString()
        else:
            raise Exception("No loss function specified")
        self.module = parameters["module"].GetString()
        self.reduction = parameters["reduction"].GetString()
        loss_module = import_module(self.module)
        if self.module == "keras":
            self.loss_function = getattr(tensorflow.keras.losses,self.loss_function_name)(reduction = self.reduction)
        else:
            self.loss_function = getattr(loss_module,loss_function)
            
    def Compile(self, loss_function, optimizer):
        """ Process for compiling a network. """
        loss_function = self.loss_function
        return [loss_function, optimizer]