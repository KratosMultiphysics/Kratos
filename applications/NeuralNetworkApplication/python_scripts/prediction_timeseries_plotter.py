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
            "node_id"             : 0,
            "variables"           : [],
            "training_timesteps"  : 100,
            "axis"                : "plot",
            "figure_size_inches_1": 8,
            "figure_size_inches_2": 6,
            "dpi"                 : 100,
            "predict_marker"      : "g-",
            "target_marker"       : "b-",
            "training_marker"     : "r-",
            "output_name"         : "",
            "output_format"       : "png",
            "only_test"           : false          
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.predictions_file = parameters["predictions_file"].GetString()
        self.target_file = parameters["target_file"].GetString()
        self.input_variable = parameters["input_variable"].GetString()
        self.variables = parameters["variables"].GetStringArray()
        self.node_id = parameters["node_id"].GetInt()
        self.figure_size_inches_1 = parameters["figure_size_inches_1"].GetInt()
        self.figure_size_inches_2 = parameters["figure_size_inches_2"].GetInt()
        self.dpi = parameters["dpi"].GetInt()
        self.predict_marker = parameters["predict_marker"].GetString()
        self.target_marker = parameters["target_marker"].GetString()
        self.training_marker = parameters["training_marker"].GetString()
        self.output_format = parameters["output_format"].GetString()
        self.output_name = parameters["output_name"].GetString()
        self.axis = parameters["axis"].GetString()
        self.training_timesteps = parameters["training_timesteps"].GetInt()
        self.only_test = parameters["only_test"].GetBool()

    def Plot(self):

        target = ImportDataFromFile(self.target_file, "OutputData").ExportAsArray()
        if self.predictions_file.endswith('.npy'):
            predictions = np.load(self.predictions_file)
            predictions = np.squeeze(predictions)
        else:
            predictions = np.genfromtxt(self.predictions_file)
        predictions_times = range(self.training_timesteps,len(predictions[:]))
        for variable in self.variables:
            figure, ax = plt.subplots()
            if not self.only_test:
                if isinstance(target[0],(list, tuple, np.ndarray)):
                    try:
                        getattr(ax,self.axis)(target[:,self.node_id,self.variables.index(variable)],self.target_marker,label='Ground Truth')
                    except IndexError:
                        getattr(ax,self.axis)(target[:,len(self.variables)*self.node_id+self.variables.index(variable)],self.target_marker,label='Ground Truth')
                else:
                    getattr(ax,self.axis)(target[:,len(self.variables)*self.node_id+self.variables.index(variable)],self.target_marker,label='Ground Truth')
                if isinstance(predictions[0],(list, tuple, np.ndarray)):
                    try: 
                        getattr(ax,self.axis)(range(self.training_timesteps-1),predictions[:self.training_timesteps-1,self.node_id,self.variables.index(variable)],self.training_marker,label='Training')
                    except IndexError:
                        getattr(ax,self.axis)(range(self.training_timesteps-1),predictions[:self.training_timesteps-1,len(self.variables)*self.node_id+self.variables.index(variable)],self.training_marker,label='Training')
                else:
                    try:
                        getattr(ax,self.axis)(range(self.training_timesteps-1),predictions[:self.training_timesteps-1,len(self.variables)*self.node_id+self.variables.index(variable)],self.training_marker,label='Training')
                    except IndexError:
                        getattr(ax,self.axis)(range(self.training_timesteps-1),predictions[:self.training_timesteps-1],self.training_marker,label='Training')
            if isinstance(predictions[0],(list, tuple, np.ndarray)):
                try: 
                    getattr(ax,self.axis)(predictions_times,predictions[self.training_timesteps:,self.node_id,self.variables.index(variable)],self.predict_marker,label='Prediction')
                except IndexError:
                    getattr(ax,self.axis)(predictions_times,predictions[self.training_timesteps:,len(self.variables)*self.node_id+self.variables.index(variable)],self.predict_marker,label='Prediction')
            else:
                try:
                    getattr(ax,self.axis)(predictions_times,predictions[self.training_timesteps:,len(self.variables)*self.node_id+self.variables.index(variable)],self.predict_marker,label='Prediction')
                except IndexError:
                    getattr(ax,self.axis)(predictions_times,predictions[self.training_timesteps:],self.predict_marker,label='Prediction')
            ax.set_xlabel(self.input_variable)
            ax.set_ylabel(variable)
            ax.legend(loc = 'upper left')
            figure = plt.gcf()
            figure.set_size_inches(self.figure_size_inches_1, self.figure_size_inches_2)
            figure.show()
            figure.savefig(self.output_name + "_" + variable + "." + self.output_format, bbox_inches='tight', dpi = self.dpi)

            


