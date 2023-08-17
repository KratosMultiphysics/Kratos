import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
from KratosMultiphysics.NeuralNetworkApplication.load_tunable_parameters_utility import LoadTunableParametersUtility

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
            "learning_rate"            : 0.001,
            "decay"                    : 0,
            "tunable"                  : false,
            "tunable_variable"         : "",
            "tuning_parameters"        : {}
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        self.tunable = parameters["tunable"].GetBool()

        if parameters.Has("optimizer"): 
            self.optimizer_name = parameters["optimizer"].GetString()
        elif not self.tunable:
            raise Exception("No optimizer specified")
        self.module = parameters["module"].GetString()
        self.learning_rate = parameters["learning_rate"].GetDouble()
        self.decay = parameters["decay"].GetDouble()
        optimizer_module = import_module(self.module)
        if self.module == "keras":
            self.optimizer = getattr(tensorflow.keras.optimizers.legacy,self.optimizer_name)(learning_rate = self.learning_rate, decay = self.decay)
        else:
            self.optimizer = getattr(optimizer_module,self.optimizer_name)(learning_rate = self.learning_rate)

        # The tuner can only tune the LR or the optimizer model at a time
        if self.tunable:
            self.tunable_variable = parameters["tunable_variable"].GetString()
            self.tuning_parameters = parameters["tuning_parameters"] 

    def Compile(self, loss_function, optimizer, hp = None):
        """ Process for compiling a network. """
        if self.tunable:
            setattr(self,self.tunable_variable,LoadTunableParametersUtility(self.tuning_parameters).Load(hp))
        optimizer = self.optimizer
        return [loss_function, optimizer]