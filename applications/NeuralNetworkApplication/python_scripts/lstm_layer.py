import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_layer import NeuralNetworkLayerClass
from KratosMultiphysics.NeuralNetworkApplication.load_tunable_parameters_utility import LoadTunableParametersUtility

from tensorflow.keras import layers

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return LSTMLayer(settings)

class LSTMLayer(NeuralNetworkLayerClass):
    """This class generates a base class for a Dense layer

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
            "layer_name"               : "",
            "trainable"                : true,
            "dtype"                    : "",
            "dynamic"                  : false,
            "units"                    : 1,
            "activation"               : "tanh",
            "recurrent_activation"     : "sigmoid",
            "use_bias"                 : true,
            "kernel_initializer"       : "glorot_uniform",
            "recurrent_initializer"    : "orthogonal",
            "bias_initializer"         : "zeros",
            "unit_forget_bias"         : true,
            "kernel_regularizer"       : 0.0,
            "recurrent_regularizer"    : 0.0,
            "bias_regularizer"         : 0.0,
            "activity_regularizer"     : 0.0,
            "dropout"                  : 0.0,
            "recurrent_dropout"         : 0.0,
            "return_sequences"         : false,
            "return_state"             : false,
            "go_backwards"             : false,
            "stateful"                 : false,
            "time_mayor"               : false,
            "unroll"                   : false,
            "tunable"                  : false,
            "tunable_variable"         : "",
            "tuning_parameters"        : {}
        }""")
    # Constraints are currently not supported. Every supported constrain in keras should be implemented individually
        settings.ValidateAndAssignDefaults(default_settings)

        self.layer_name = settings["layer_name"].GetString()
        self.trainable = settings["trainable"].GetBool()
        self.dtype = settings["dtype"].GetString()
        self.dynamic = settings["dynamic"].GetBool()
        self.units = settings["units"].GetInt()
        if not settings["activation"].GetString() == "":
            self.activation = settings["activation"].GetString()
        else:
            self.activation = None
        if not settings["recurrent_activation"].GetString() == "":
            self.recurrent_activation = settings["recurrent_activation"].GetString()
        else:
            self.recurrent_activation = None
        self.use_bias = settings["use_bias"].GetBool()
        self.unit_forget_bias = settings["unit_forget_bias"].GetBool()
        self.recurrent_initializer = settings["recurrent_initializer"].GetString()
        self.recurrent_regularizer = settings["recurrent_regularizer"].GetDouble()
        if abs(self.recurrent_regularizer) < 1e-8:
            self.recurrent_regularizer = None
        self.kernel_initializer = settings["kernel_initializer"].GetString()
        self.kernel_regularizer = settings["kernel_regularizer"].GetDouble()
        if abs(self.kernel_regularizer) < 1e-8:
            self.kernel_regularizer = None
        self.bias_initializer = settings["bias_initializer"].GetString()
        self.bias_regularizer = settings["bias_regularizer"].GetDouble()
        if abs(self.bias_regularizer) < 1e-8:
            self.bias_regularizer = None
        self.activity_regularizer = settings["activity_regularizer"].GetDouble()
        if abs(self.activity_regularizer) < 1e-8:
            self.activity_regularizer = None
        self.dropout = settings["dropout"].GetDouble()
        self.recurrent_dropout = settings["recurrent_dropout"].GetDouble()
        self.return_sequences = settings["return_sequences"].GetBool()
        self.return_state = settings["return_state"].GetBool()
        self.go_backwards = settings["go_backwards"].GetBool()
        self.stateful = settings["stateful"].GetBool()
        self.time_major = settings["time_mayor"].GetBool()
        self.unroll = settings["unroll"].GetBool()

        self.tunable = settings["tunable"].GetBool()
        
        if self.tunable:
            self.tunable_variable = settings["tunable_variable"].GetString()
            self.tuning_parameters = settings["tuning_parameters"]

    def Build(self, hp = None):
        """ This method builds the layer when called by a process."""
        if self.tunable:
            setattr(self,self.tunable_variable,LoadTunableParametersUtility(self.tuning_parameters).Load(hp))
        
        self.layer = layers.LSTM(self.units, activation=self.activation,recurrent_activation = self.recurrent_activation, 
        use_bias=self.use_bias, kernel_initializer=self.kernel_initializer, recurrent_initializer = self.recurrent_initializer,
        bias_initializer=self.bias_initializer,unit_forget_bias=self.unit_forget_bias,kernel_regularizer=self.kernel_regularizer,
        recurrent_regularizer=self.recurrent_regularizer,bias_regularizer=self.bias_regularizer,
        activity_regularizer=self.activity_regularizer,dropout=self.dropout,recurrent_dropout=self.recurrent_dropout,
        return_sequences=self.return_sequences,return_state=self.return_state,go_backwards=self.go_backwards,
        stateful=self.stateful,time_major=self.time_major,unroll=self.unroll,
        trainable=self.trainable, dtype=self.dtype, name=self.layer_name, dynamic=self.dynamic)
        return self.layer



