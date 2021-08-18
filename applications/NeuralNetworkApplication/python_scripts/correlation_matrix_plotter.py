from tensorflow.python.ops.variables import Variable
import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDataFromFile
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sn



def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return CorrelationMatrixPlotterProcess(settings["parameters"])

class CorrelationMatrixPlotterProcess(NeuralNetworkProcess):

    def __init__(self,parameters):

        super().__init__()

        default_settings = KM.Parameters("""{
            "file"                : "",
            "origin"              : "output_generation",
            "bias"                : true,
            "annotations"         : false,
            "variables"           : [],
            "vmin"                : [],
            "vmax"                : [],
            "output_name"         : "",
            "output_format"       : "png"          
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.file = parameters["file"].GetString()
        self.origin = parameters["origin"].GetString() # origin = input_generation, output_generation or prediction
        self.bias = parameters["bias"].GetBool()
        self.vmin = parameters["vmin"].GetVector()
        self.vmax = parameters["vmax"].GetVector()
        self.annotations = parameters["annotations"].GetBool()
        self.variables = parameters["variables"].GetStringArray()
        self.output_format = parameters["output_format"].GetString()
        self.output_name = parameters["output_name"].GetString()

    def Plot(self):
        
        if self.origin == "input_generation" or self.origin == "output_generation":
            data = ImportDataFromFile(self.file, "OutputData")
        elif self.origin == "prediction":
            if self.file.endswith('.npy'):
                data = np.load(self.file)
            else:
                data = np.genfromtxt(self.file)
        else:
            raise Exception("Data origin must be input_generation, output_generation or prediction.")
       
        for variable in self.variables:
            try:
                value_vmin = self.vmin[self.variables.index(variable)] 
            except IndexError:
                value_vmin = None
            try:
                value_vmax = self.vmax[self.variables.index(variable)] 
            except IndexError:
                value_vmax = None

            figure, ax = plt.subplots()
            if isinstance(data[0], (list, tuple, np.ndarray)):
                try:
                    cov_matrix = np.cov(np.transpose(data[:,:,self.variables.index(variable)]), bias = self.bias)
                except IndexError:
                    cov_matrix = np.cov(np.transpose(data[:,self.variables.index(variable):][:,::len(self.variables)]), bias = self.bias)
                sn_plot = sn.heatmap(cov_matrix, annot = self.annotations, fmt='g', vmin = value_vmin, vmax = value_vmax)
                figure = sn_plot.get_figure()
                figure.savefig(self.output_name + "_" + variable + "." + self.output_format)

            


