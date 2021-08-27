import numpy as np
import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.input_dataclasses import ListNeuralNetworkData
from KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDictionaryFromText, UpdateDictionaryJson, KratosVectorToList


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DivideDatasetProcess(settings["parameters"])

class DivideDatasetProcess(PreprocessingProcess):

    def __init__(self, settings):
        super().__init__(settings)
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        try:
            self.proportion = settings["proportion"].GetDouble()
        except AttributeError:
            print("No proportion indicated for the division of the dataset. The chosen value is 0.0.")
        
        if (self.proportion >= 1.0) or (self.proportion <= 0.0):
            print("The division proportion must be a value between 0.0 and 1.0.")

    def Preprocess(self, data_structure_in, data_structure_out):
        """ This method divides the dataset in two parts """

        data_structure_in.UpdateData(data_structure_in[:np.int64(len(data_structure_in)*self.proportion)])
        data_structure_out.UpdateData(data_structure_out[:np.int64(len(data_structure_out)*self.proportion)])

        return [data_structure_in, data_structure_out]
