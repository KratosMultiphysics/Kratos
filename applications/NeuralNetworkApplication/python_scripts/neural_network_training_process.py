from importlib import import_module
import KratosMultiphysics as KM
import numpy as np

def Factory(settings):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return NeuralNetworkTrainingProcess(settings["parameters"])


class NeuralNetworkTrainingProcess(KM.Process):
    """
    This class is the parent for the classes that deal with the training of the neural network.

    It is similar to the normal Kratos Process but it does not receive a model, only settings.
    """
    def __init__(self, parameters):
        super().__init__()
        self.input_file = parameters["data_input"].GetString()
        test_input_raw =[]
        with open(self.input_file,'r') as test_input_file:
            for line in test_input_file:
                line = line.strip()
                if not line:
                    continue
                test_input_raw.append(list(map(str, line.split())))

        test_input_raw_splits =[]
        for line in test_input_raw:
            test_input_raw_splits_line = []
            for data in line:
                #dimension = int(data[1])
                array = ''.join(data[4:-1])
                array_floats = list(map(float,array.split(',',)))
                test_input_raw_splits_line.append(array_floats)
            test_input_raw_splits.append(test_input_raw_splits_line)
        test_input_raw_splits = np.array(test_input_raw_splits)

        self.test_input = test_input_raw_splits[:,:,1]

        # Output
        self.output_file = parameters["data_output"].GetString()
        test_output_raw =[]
        with open(self.output_file,'r') as test_output_file:
            for line in test_output_file:
                line = line.strip()
                if not line:
                    continue
                test_output_raw.append(list(map(str, line.split())))

        test_output_raw_splits =[]
        for line in test_output_raw:
            test_output_raw_splits_line = []
            for data in line:
                #dimension = int(data[1])
                array = ''.join(data[4:-1])
                array_floats = list(map(float,array.split(',',)))
                test_output_raw_splits_line.append(array_floats)
            test_output_raw_splits.append(test_output_raw_splits_line)
        test_output_raw_splits = np.array(test_output_raw_splits)

        self.test_output = test_output_raw_splits[:,:,2]

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