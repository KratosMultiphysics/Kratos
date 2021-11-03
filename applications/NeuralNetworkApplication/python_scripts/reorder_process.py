import numpy as np
import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
from KratosMultiphysics.NeuralNetworkApplication.input_dataclasses import ListDataWithLookback, DataWithLookback


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ReorderProcess(settings["parameters"])

class ReorderProcess(PreprocessingProcess):

    def __init__(self, settings):
        super().__init__(settings)
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        if not self.load_from_log:
            try:
                self.number_of_partitions = settings["number_of_partitions"].GetInt()
            except RuntimeError:
                self.number_of_partitions = 2            

    def Preprocess(self, data_structure_in, data_structure_out):
        """ This method modifies the order that the data are exported from the class.
        By default, data + lookback are directly concatenated. This method allows to alter the order."""

        data_structure_in.SetReorder(self.number_of_partitions)

        return [data_structure_in, data_structure_out]