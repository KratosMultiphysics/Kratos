import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_layer import NeuralNetworkLayerClass
from KratosMultiphysics.NeuralNetworkApplication.load_tunable_parameters_utility import LoadTunableParametersUtility
import KratosMultiphysics.NeuralNetworkApplication.add_layer_process as AddLayer

from tensorflow.keras import layers
import tensorflow.keras as keras
import kerastuner as kt

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return BlockDenseLayer(settings)

class BlockDenseLayer(NeuralNetworkLayerClass):
    """This class generates a block of Dense layers

    Public member variables:
    settings -- Kratos parameters containing process settings.
    """
    def __init__(self, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        settings -- Kratos parameters containing layer settings.
        """
        default_settings = KM.Parameters("""{
            "number_layers"            : 1,
            "neurons_per_layer"        : [],
            "batch_normalization"      : true,
            "activation"               : "",
            "layer_parameters"         : {},
            "tunable"                  : false,
            "tunable_variable"         : "",
            "tuning_parameters"        : {}
        }""")
   
        settings.ValidateAndAssignDefaults(default_settings)

        self.tunable = settings["tunable"].GetBool()
        self.number_layers = settings["number_layers"].GetInt()
        self.neurons_per_layer = settings["neurons_per_layer"].GetVector()
        self.batch_normalization = settings["batch_normalization"].GetBool()
        if not settings["activation"].GetString() == "":
            self.activation = settings["activation"].GetString()
        else:
            self.activation = None
        self.layer_parameters = settings["layer_parameters"]
        if self.tunable:
            self.tunable_variable = settings["tunable_variable"].GetString()
            self.tuning_parameters = settings["tuning_parameters"]

    def Build(self, hp = None):
        """ This method builds the layer when called by a process."""
        if self.tunable:
            setattr(self,self.tunable_variable,LoadTunableParametersUtility(self.tuning_parameters).Load(hp))
        self.layer_block = keras.Sequential()
        for i in range(self.number_layers):
            dense_layer_parameters = KM.Parameters()
            dense_layer_parameters.AddEmptyValue("parameters")
            dense_layer_parameters["parameters"].AddEmptyValue("layer_type")
            dense_layer_parameters["parameters"]["layer_type"].SetString("dense_layer")
            dense_layer_parameters["parameters"].AddEmptyArray("layer_parameters")
            dense_layer_parameters["parameters"]["layer_parameters"] = self.layer_parameters
            if not dense_layer_parameters["parameters"]["layer_parameters"].Has("tunable"):
                dense_layer_parameters["parameters"]["layer_parameters"].AddEmptyValue("tunable")
                dense_layer_parameters["parameters"]["layer_parameters"]["tunable"].SetBool(False)
            if dense_layer_parameters["parameters"]["layer_parameters"].Has("units") and not dense_layer_parameters["parameters"]["layer_parameters"]["tunable"].GetBool():
                raise Exception("The units per layer must be defined on the block level.")
            if not dense_layer_parameters["parameters"]["layer_parameters"]["tunable"].GetBool():
                dense_layer_parameters["parameters"]["layer_parameters"].AddEmptyValue("units")
                dense_layer_parameters["parameters"]["layer_parameters"]["units"].SetInt(int(self.neurons_per_layer[i]))
            else:
                dense_layer_parameters["parameters"]["layer_parameters"]["tuning_parameters"].AddEmptyValue("name")
                dense_layer_parameters["parameters"]["layer_parameters"]["tuning_parameters"]["name"].SetString("units_"+str(i))
            dense_layer = AddLayer.Factory(dense_layer_parameters)
        
            self.layer_block = dense_layer.Add(self.layer_block, hp = hp)

            if self.batch_normalization:
                normalization_layer_parameters = KM.Parameters()
                normalization_layer_parameters.AddEmptyValue("parameters")
                normalization_layer_parameters["parameters"].AddEmptyValue("layer_type")
                normalization_layer_parameters["parameters"]["layer_type"].SetString("batch_normalization_layer")
                normalization_layer_parameters["parameters"].AddEmptyArray("layer_parameters")
                normalization_layer_parameters["parameters"]["layer_parameters"]=KM.Parameters()
                normalization_layer = AddLayer.Factory(normalization_layer_parameters)
                self.layer_block = normalization_layer.Add(self.layer_block, hp = hp)

            if not self.activation == None:
                activation_layer_parameters = KM.Parameters()
                activation_layer_parameters.AddEmptyValue("parameters")
                activation_layer_parameters["parameters"].AddEmptyValue("layer_type")
                activation_layer_parameters["parameters"]["layer_type"].SetString("activation_layer")
                activation_layer_parameters["parameters"].AddEmptyArray("layer_parameters")
                activation_layer_parameters["parameters"]["layer_parameters"]=KM.Parameters()
                activation_layer_parameters["parameters"]["layer_parameters"].AddEmptyValue("activation")
                activation_layer_parameters["parameters"]["layer_parameters"]["activation"].SetString(self.activation)
                activation_layer = AddLayer.Factory(activation_layer_parameters)
                self.layer_block = activation_layer.Add(self.layer_block, hp = hp)
        
        return self.layer_block



