import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
from KratosMultiphysics.NeuralNetworkApplication.load_tunable_parameters_utility import LoadTunableParametersUtility
import numpy as np
import tensorflow.keras.losses
from importlib import import_module


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return CustomLossFunctionCompilerProcess(settings["parameters"])


class CustomLossFunctionCompilerProcess(NeuralNetworkProcess):
    """
    This class compiles a custom loss function.

    """
    def __init__(self, parameters):
        super().__init__()

        default_settings = KM.Parameters("""{
            "help"                     : "This process compiles a custom loss function.",
            "loss_function"            : "",
            "module"                   : "keras",
            "reduction"                : "sum_over_batch_size",
            "tunable"                  : false,
            "tunable_variable"         : "",
            "tuning_parameters"        : {}
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        self.tunable = parameters["tunable"].GetBool()
        if parameters.Has("loss_function"): 
            self.loss_function_name = parameters["loss_function"].GetString()
        elif not parameters["tunable"]:
            raise Exception("No loss function specified")
        self.module = parameters["module"].GetString()
        self.reduction = parameters["reduction"].GetString()
        loss_module = import_module(self.module)
        if self.module == "keras":
            self.loss_function = getattr(tensorflow.keras.losses,self.loss_function_name)(reduction = self.reduction)
        else:
            self.loss_function = getattr(loss_module,self.loss_function_name)

        if self.tunable:
            self.tunable_variable = parameters["tunable_variable"].GetString()
            self.tuning_parameters = parameters["tuning_parameters"]
            
    def Compile(self, loss_function, optimizer, hp = None):
        """ Process for compiling a network. """
        if self.tunable:
            setattr(self,self.tunable_variable,LoadTunableParametersUtility(self.tuning_parameters).Load(hp))
        loss_function = self.loss_function
        return [loss_function, optimizer]