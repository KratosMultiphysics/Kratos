from importlib import import_module
import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities as DataLoadingUtilities
import numpy as np
import h5py

def Factory(settings):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return NeuralNetworkTrainingProcess(settings["parameters"])


class NeuralNetworkTrainingProcess(NeuralNetworkProcess):
    """
    This class is the parent for the classes that deal with the training of the neural network.

    It is similar to the normal Kratos Process but it does not receive a model, only settings.
    """
    def __init__(self, parameters):
        super().__init__()
        self.input_file = parameters["data_input"].GetString()
        # Data loading for h5
        if self.input_file.endswith('.h5'):
            self.test_input = DataLoadingUtilities.ImportH5(self.input_file, "InputData")
        #   with h5py.File(self.input_file,'r') as f:
                # for time in f.keys():
                #     time_array = []
                #     for variables in f[time]['InputData']['NodalSolutionStepData'].keys():
                #         value = f[time]['InputData']['NodalSolutionStepData'][variables][:]
                #         # time_array.append(value[:])
                #         time_array = value[:]
                #     test_input_raw.append(time_array)
                # self.test_input = np.array(test_input_raw)
        # Data loading for dat
        elif self.input_file.endswith('.dat'):
            self.test_input = DataLoadingUtilities.ImportAscii(self.input_file)
            # with open(self.input_file,'r') as test_input_file:
            #     for line in test_input_file:
            #         line = line.strip()
            #         if not line:
            #             continue
            #         test_input_raw.append(list(map(str, line.split())))
            # test_input_raw_splits =[]
            # for line in test_input_raw:
            #     test_input_raw_splits_line = []
            #     for data in line:
            #         if len(data.split(' '))>1:
            #             array = ''.join(data[4:-1])
            #             array_floats = list(map(float,array.split(',',)))
            #             test_input_raw_splits_line.append(array_floats)
            #         else:
            #             array_floats = float(data)
            #             test_input_raw_splits_line.append(array_floats)
            #     test_input_raw_splits.append(test_input_raw_splits_line)
            # test_input_raw_splits = np.array(test_input_raw_splits)
            # self.test_input = test_input_raw_splits
        # Exception for non-supported formats
        else:
            raise Exception("Input data format not supported. Supported formats are .dat and .h5")

        # Output
        self.output_file = parameters["data_output"].GetString()
        # Data loading for h5
        if self.output_file.endswith('.h5'):
            self.test_output = DataLoadingUtilities.ImportH5(self.output_file,"OutputData")
        #   with h5py.File(self.output_file,'r') as f:
        #         for time in f.keys():
        #             time_array = []
        #             for variables in f[time]['OutputData']['NodalSolutionStepData'].keys():
        #                 value = f[time]['OutputData']['NodalSolutionStepData'][variables][:]
        #                 # time_array.append(value[:])
        #                 time_array = value[:]
        #             test_output_raw.append(time_array)
        #         self.test_output = np.array(test_output_raw)
        # # Data loading for dat
        elif self.output_file.endswith('.dat'):
            self.test_output = DataLoadingUtilities.ImportAscii(self.output_file)
            # with open(self.output_file,'r') as test_output_file:
            #     for line in test_output_file:
            #         line = line.strip()
            #         if not line:
            #             continue
            #         test_output_raw.append(list(map(str, line.split())))
            # test_output_raw_splits =[]
            # for line in test_output_raw:
            #     test_output_raw_splits_line = []
            #     for data in line:
            #         if len(data.split(' '))>1:
            #             array = ''.join(data[4:-1])
            #             array_floats = list(map(float,array.split(',',)))
            #             test_output_raw_splits_line.append(array_floats)
            #         else:
            #             array_floats = float(data)
            #             test_output_raw_splits_line.append(array_floats)
            #     test_output_raw_splits.append(test_output_raw_splits_line)
            # test_output_raw_splits = np.array(test_output_raw_splits)
            # self.test_output = test_output_raw_splits
        else:
            raise Exception("Output data format not supported. Supported formats are .dat and .h5")

