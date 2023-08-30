import KratosMultiphysics as KM
# import KratosMultiphysics.NeuralNetworkApplication as NeuralNetwork

from tensorflow.keras import layers

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return NeuralNetworkLayerClass(settings)

class NeuralNetworkLayerClass:
    """This class generates a base class for a Neural Network layer

    Public member variables:
    settings -- Kratos parameters containing process settings.
    """

    def __init__(self,settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        settings -- Kratos parameters containing layer settings.
        """

        default_settings = KM.Parameters("""{
            "layer_name"               : "",
            "trainable"                : true,
            "dtype"                    : "",
            "dynamic"                  : false
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.layer_name = settings["layer_name"].GetString()
        if self.layer_name == "":
            self.layer_name = None
        self.trainable = settings["trainable"].GetBool()
        self.dtype = settings["dtype"].GetString()
        self.dynamic = settings["dynamic"].GetBool()

    def Build(self, hp = None):
        """ This method builds the layer

            Keyword arguments:
            self -- It signifies an instance of a class.
        """
        self.layer = layers.Layer(trainable = self.trainable, name = self.layer_name, dtype = self.dtype, dynamic = self.dynamic)
        return self.layer

    def GetWeights(self):
        """ This method returns the weights of the layer

            Keyword arguments:
            self -- It signifies an instance of a class.
        """
        if hasattr(self, self.layer):
            return self.layer.get_weights()
        else:
            raise Exception('The layer is not built yet!')

    def SetWeights(self,weights):
        """ This method sets the weights of the layer

            Keyword arguments:
            self -- It signifies an instance of a class.
            weights -- It signifies the weights to set.
        """
        if hasattr(self, self.layer):
            return self.layer.set_weights(weights)
        else:
            raise Exception('The layer is not built yet!')

