from tensorflow.python.ops.variables import Variable
import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDataFromFile
import numpy as np
import matplotlib.pyplot as plt



def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return PredictionErrorPlotterProcess(settings["parameters"])

class PredictionErrorPlotterProcess(NeuralNetworkProcess):

    def __init__(self,parameters):

        super().__init__()


        default_settings = KM.Parameters("""{
            "predictions_file"     : "",
            "input_file"          : "",
            "target_file"         : "",
            "node_id"             : 0,
            "input_variable"      : "",
            "input_variable_id"   : 0,
            "timesteps"           : 0,
            "variables"           : [],
            "axis"                : "plot",
            "marker"              : "+",
            "figure_size_inches_1": 8,
            "figure_size_inches_2": 6,
            "dpi"                 : 100,
            "output_name"         : "",
            "output_format"       : "png"          
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.predictions_file = parameters["predictions_file"].GetString()
        self.input_file = parameters["input_file"].GetString()
        self.target_file = parameters["target_file"].GetString()
        self.node_id = parameters["node_id"].GetInt()
        self.input_variable = parameters["input_variable"].GetString()
        self.input_variable_id = parameters["input_variable_id"].GetInt()
        self.timesteps = parameters["timesteps"].GetInt()
        self.variables = parameters["variables"].GetStringArray()
        self.marker = parameters["marker"].GetString()
        self.figure_size_inches_1 = parameters["figure_size_inches_1"].GetInt()
        self.figure_size_inches_2 = parameters["figure_size_inches_2"].GetInt()
        self.dpi = parameters["dpi"].GetInt()
        self.output_format = parameters["output_format"].GetString()
        self.output_name = parameters["output_name"].GetString()
        self.axis = parameters["axis"].GetString()

    def Plot(self):

        input = ImportDataFromFile(self.input_file, "InputData").ExportAsArray()
        target = ImportDataFromFile(self.target_file, "OutputData").ExportAsArray()

        if self.predictions_file.endswith('.npy'):
            predictions = np.load(self.predictions_file)
        else:
            predictions = np.genfromtxt(self.predictions_file)
        if self.timesteps > 0:
            error = predictions[:self.timesteps] - target[:self.timesteps]
            input = input[:self.timesteps]
        else:
            error = predictions - target

        for variable in self.variables:
            figure, ax = plt.subplots()
            if isinstance(error[0],(list, tuple, np.ndarray)):
                try:
                    getattr(ax,self.axis)(input[:,self.node_id + self.input_variable_id],error[:,self.node_id + self.variables.index(variable)],self.marker)
                except IndexError:
                    getattr(ax,self.axis)(input,error[:,self.node_id + self.variables.index(variable)],self.marker)
            else:
                try:
                    getattr(ax,self.axis)(input[:,self.node_id + self.input_variable_id],error[:],self.marker)
                except IndexError:
                    getattr(ax,self.axis)(input,error[:],self.marker)

            ax.set_xlabel(self.input_variable)
            ax.set_ylabel(variable)
            figure = plt.gcf() 
            figure.set_size_inches(self.figure_size_inches_1, self.figure_size_inches_2)
            figure.show()
            figure.savefig(self.output_name + "_" + variable + "." + self.output_format, bbox_inches='tight', dpi = self.dpi)

            


