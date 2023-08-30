import KratosMultiphysics as KM
import KratosMultiphysics.NeuralNetworkApplication as NeuralNetwork
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess

from importlib import import_module
from tensorflow.keras import layers

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SaveModelProcess(settings["parameters"])

class SaveModelProcess(NeuralNetworkProcess):

    def __init__(self, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        super().__init__()
        default_settings = KM.Parameters("""{
            "help"                     : "This process saves a network.",
            "model_name"               : "",
            "format"                   : ""
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_name = settings["model_name"].GetString()
        self.format = settings["format"].GetString()


    def Save(self,model):
        """ Processes for saving a network. """
        if self.format == "SavedModel":
            model.save(self.model_name)
        elif self.format == "h5":
            model.save(self.model_name + ".h5")
        else:
            raise Exception("Saving format not recognized. Available formats are SavedModel and h5.")
