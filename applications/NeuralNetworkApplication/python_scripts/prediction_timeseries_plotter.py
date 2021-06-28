from tensorflow.python.ops.variables import Variable
import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDataFromFile
import numpy as np
import matplotlib.pyplot as plt



def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return PredictionTimeseriesPlotterProcess(settings["parameters"])

class PredictionTimeseriesPlotterProcess(NeuralNetworkProcess):

    def __init__(self,parameters):

        super().__init__()


        default_settings = KM.Parameters("""{
            "predictions_file"     : "",
            "target_file"         : "",
            "input_variable"      : "",
            "variables"           : [],
            "training_timesteps"  : 100,
            "axis"                : "plot",
            "output_name"         : "",
            "output_format"       : "png"          
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.predictions_file = parameters["predictions_file"].GetString()
        self.target_file = parameters["target_file"].GetString()
        self.input_variable = parameters["input_variable"].GetString()
        self.variables = parameters["variables"].GetStringArray()
        self.output_format = parameters["output_format"].GetString()
        self.output_name = parameters["output_name"].GetString()
        self.axis = parameters["axis"].GetString()
        self.training_timesteps = parameters["training_timesteps"].GetInt()

    def ExecuteFinalize(self):

        target = ImportDataFromFile(self.target_file, "OutputData")
        if self.predictions_file.endswith('.npy'):
            predictions = np.load(self.predictions_file)
        else:
            predictions = np.genfromtxt(self.predictions_file)
        predictions_times = range(self.training_timesteps,len(predictions[:]))
        for variable in self.variables:
            figure, ax = plt.subplots()
            getattr(ax,self.axis)(target[:self.training_timesteps-1,self.variables.index(variable)],'b-',label='Ground Truth')
            if isinstance(predictions[0],(list, tuple, np.ndarray)): 
                getattr(ax,self.axis)(range(self.training_timesteps-1),predictions[:self.training_timesteps-1,self.variables.index(variable)],'r-',label='Training')
            else:
                getattr(ax,self.axis)(range(self.training_timesteps-1),predictions[:self.training_timesteps-1],'r-',label='Training')
            if isinstance(predictions[0],(list, tuple, np.ndarray)): 
                getattr(ax,self.axis)(predictions_times,predictions[self.training_timesteps:,self.variables.index(variable)],'g-',label='Prediction')
            else:
                getattr(ax,self.axis)(predictions_times,predictions[self.training_timesteps:],'g-',label='Prediction')
            ax.set_xlabel(self.input_variable)
            ax.set_ylabel(variable)
            ax.legend()
            manager = plt.get_current_fig_manager()
            manager.resize(*manager.window.maxsize())
            figure.show()
            figure.savefig(self.output_name + "_" + variable + "." + self.output_format, bbox_inches='tight')

            


