import KratosMultiphysics as KM
import KratosMultiphysics.NeuralNetworkApplication as NeuralNetwork
from KratosMultiphysics.NeuralNetworkApplication.neural_network_training_process import NeuralNetworkTrainingProcess

from importlib import import_module
from tensorflow.keras import layers

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AddLayerProcess(settings["parameters"])

class AddLayerProcess(KM.Process):

    def __init__(self, settings):
        super().__init__()
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """

        default_settings = KM.Parameters("""{
            "help"                     : "This process adds a layer to a network.",
            "layer_type"               : "",
            "layer_parameters"         : {}
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.layer_type = settings["layer_type"].GetString()
        self.layer_parameters = settings["layer_parameters"]
        # Import layer classes
        module_name = "KratosMultiphysics.NeuralNetworkApplication." + self.layer_type
        class_module = import_module(module_name)
        layer_class = class_module.Factory(self.layer_parameters)

        # Build layer
        self.layer = layer_class.Build()

    def Add(self,model):
        return self.layer(model)

    def Save(self,model):
        """ Process for saving a network. """
        pass
    def ExecuteInitialize(self):
        """ Processes to act on the initialization. """
        pass

    def ExecuteFinalize(self):
        """ Processes to act on the finalization. """
        pass

    def ExecuteBeforeTraining(self):
        """ Processes to act just before the training. """
        pass
    
    def ExecuteTraining(self):
        """ Processes to act directly during the training step. """
        pass
