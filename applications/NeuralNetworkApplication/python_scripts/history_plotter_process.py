import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities as DataLoadingUtilities
import numpy as np
import matplotlib.pyplot as plt


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return HistoryPlotterProcess(settings["parameters"])

class HistoryPlotterProcess(NeuralNetworkProcess):

    def __init__(self,parameters):

        super().__init__()


        default_settings = KM.Parameters("""{
            "history_file"               : "",
            "variables"                  : [],
            "axis"                       : "plot",
            "output_name"                : "",
            "output_format"              : "png",
            "burn_up"                    : 0,
            "show"                       : true                     
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.history_file = parameters["history_file"].GetString()
        self.variables = parameters["variables"].GetStringArray()
        self.output_format = parameters["output_format"].GetString()
        self.output_name = parameters["output_name"].GetString()
        self.axis = parameters["axis"].GetString()
        self.burn_up = parameters["burn_up"].GetInt()
        self.show = parameters["show"].GetBool()

    def ExecuteFinalize(self):

        dictionary = DataLoadingUtilities.ImportDictionaryFromText(self.history_file)

        for variable in self.variables:
            figure, ax = plt.subplots() 
            getattr(ax,self.axis)(dictionary[variable][self.burn_up:])
            ax.set_xlabel("Epochs")
            ax.set_ylabel(variable)
            figure.savefig(self.output_name + "_" + variable + "." + self.output_format)
            if self.show:
                plt.show()
            
