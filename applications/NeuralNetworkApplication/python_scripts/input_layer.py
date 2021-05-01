import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_layer import NeuralNetworkLayerClass

from tensorflow.keras import layers
import numpy as np
import h5py

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return InputLayer(settings)

class InputLayer(NeuralNetworkLayerClass):
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
            "data_input"               : "",
            "batch_size"               : "",
            "sparse"                   : false,
            "ragged"                   : false
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        # Get the input shape from a data input
        self.input_file = settings["data_input"].GetString()
        test_input_raw =[]

        if self.input_file.endswith(".h5"):
            with h5py.File(self.input_file,'r') as f:
                time = "1"
                time_array = []
                for variables in f[time]['InputData']['NodalSolutionStepData'].keys():
                    value = f[time]['InputData']['NodalSolutionStepData'][variables][:]
                    time_array.append(value)
                test_input_raw.append(time_array)
                self.test_input = test_input_raw
                self.shape = self.test_input[0][0].shape
        elif self.input_file.endswith(".dat"):
            with open(self.input_file,'r') as test_input_file:
                line = test_input_file.readline().strip()
                test_input_raw.append(list(map(str, line.split())))

            test_input_raw_splits =[]
            for line in test_input_raw:
                test_input_raw_splits_line = []
                for data in line:
                    if len(data.split(' '))>1:
                        array = ''.join(data[4:-1])
                        array_floats = list(map(float,array.split(',',)))
                        test_input_raw_splits_line.append(array_floats)
                    else:
                        array_floats = float(data)
                        test_input_raw_splits_line.append(array_floats)
                test_input_raw_splits.append(test_input_raw_splits_line)
            test_input_raw_splits = np.array(test_input_raw_splits)
            self.test_input = test_input_raw_splits
            self.shape = self.test_input[0,:].shape
        else:
            raise Exception("Data input format not supported. Supported formats are .dat and .h5")

        self.layer_name = settings["layer_name"].GetString()
        if not settings["batch_size"].GetString() == '':
            self.batch_size = settings["batch_size"].GetInt()
        else:
            self.batch_size = None
        self.sparse = settings["sparse"].GetBool()
        # self.tensor = settings["tensor"].GetString()
        self.ragged = settings["ragged"].GetBool()

        # When called by the add_layer_process, input the parameters in the keras function

    def Build(self, hp = None):
        self.layer = layers.Input(shape = self.shape, batch_size = self.batch_size, sparse = self.sparse,ragged=self.ragged,name=self.layer_name)
        return self.layer

