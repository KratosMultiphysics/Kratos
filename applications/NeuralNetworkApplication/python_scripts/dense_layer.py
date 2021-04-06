import KratosMultiphysics as KM
import KratosMultiphysics.NeuralNetworkApplication as NeuralNetwork

from tensorflow.keras import layers

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DenseLayer(settings["Parameters"])

class DenseLayer(NeuralNetwork.NeuralNetworkLayer):
    """This class generates a base class for a Dense layer

    Public member variables:
    settings -- Kratos parameters containing process settings.
    """
    def __init__(self, settings)
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        settings -- Kratos parameters containing layer settings.
        """

        default_settings = KM.Parameters("""{
            "layer_name"               : "",
            "trainable"                : True,
            "dtype"                    : "",
            "dynamic"                  : False
            "units"                    : 1,
            "activation"               : "",
            "use_bias"                 : True,
            "kernel_initializer"       : "glorot-uniform",
            "bias_initializer"         : "zeros",
            "kernel_regularizer"       : "",
            "bias_regularizer"         : "",
            "activity_regularizer"     : "",
            "kernel_constraint"        : "",
            "bias_constraint"          : ""
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.layer_name = settings["layer_name"].GetString()
        self.trainable = settings["trainable"].GetBool()
        self.dtype = settings["dtype"].GetString()
        self.dynamic = settings["dynamic"].GetBool()
        self.units = settings["units"].GetInt()
        self.activation = settings["activation"].GetString()
        self.use_bias = settings["use_bias"].GetBool()
        self.kernel_initializer = settings["kernel_initializer"].GetString()
        self.kernel_regularizer = settings["kernel_regularizer"].GetString()
        self.kernel_constraint = settings["kernel_constraint"].GetString()
        self.bias_initializer = settings["bias_initializer"].GetString()
        self.bias_regularizer = settings["bias_regularizer"].GetString()
        self.bias_constraint = settings["bias_constraint"].GetString()
        self.activity_regularizer = settings["activity_regularizer"].GetString()

        # When called by the add_layer_process, input the parameters in the keras function




