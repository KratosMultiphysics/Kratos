from tensorflow.python.ops.variables import Variable
import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDataFromFile
import numpy as np
import matplotlib.pyplot as plt



def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return FrequencyPlotterProcess(settings["parameters"])

class FrequencyPlotterProcess(NeuralNetworkProcess):

    def __init__(self,parameters):

        super().__init__()


        default_settings = KM.Parameters("""{
            "predictions_file"     : "",
            "target_file"         : "",
            "plot_zero_frequency" : true,
            "mirror_negative_frequencies": true,
            "node_id"             : 0,
            "variables"           : [],
            "training_timesteps"  : 100,
            "sampling_rate"       : 0.01,
            "axis"                : "plot",
            "output_name"         : "",
            "output_format"       : "png"          
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.predictions_file = parameters["predictions_file"].GetString()
        self.target_file = parameters["target_file"].GetString()
        self.plot_zero_frequency = parameters["plot_zero_frequency"].GetBool()
        self.mirror_negative_frequencies = parameters["mirror_negative_frequencies"].GetBool()
        self.variables = parameters["variables"].GetStringArray()
        self.node_id = parameters["node_id"].GetInt()
        self.output_format = parameters["output_format"].GetString()
        self.output_name = parameters["output_name"].GetString()
        self.axis = parameters["axis"].GetString()
        self.training_timesteps = parameters["training_timesteps"].GetInt()
        self.sampling_rate = parameters["sampling_rate"].GetDouble()

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
            if isinstance(target[0],(list, tuple, np.ndarray)):
                try:
                    data = target[self.training_timesteps:,self.node_id,self.variables.index(variable)]
                except IndexError:
                    data = target[:self.training_timesteps-1,len(self.variables)*self.node_id+self.variables.index(variable)]
            else:
                data = target[:self.training_timesteps-1,len(self.variables)*self.node_id+self.variables.index(variable)]
            if not self.plot_zero_frequency:
                data = data - np.mean(data)
            abs_fourier_transform = abs(np.fft.rfft(data))
            frequency = np.linspace(0, 1.0/self.sampling_rate/2, len(abs_fourier_transform))
            if self.mirror_negative_frequencies:
                neg_abs_fourier_transform = np.array([x for x in reversed(abs_fourier_transform)])
                abs_fourier_transform = neg_abs_fourier_transform + abs_fourier_transform
                neg_frequency = np.array([-x for x in reversed(frequency)])
                frequency = neg_frequency + frequency
            getattr(ax,self.axis)(frequency, abs_fourier_transform,'b-',label='Ground Truth')

            if isinstance(predictions[0],(list, tuple, np.ndarray)):
                try: 
                    data = predictions[self.training_timesteps:,self.node_id,self.variables.index(variable)]
                except IndexError:
                    data = predictions[self.training_timesteps:,len(self.variables)*self.node_id+self.variables.index(variable)]
            else:
                try:
                    data = predictions[self.training_timesteps:,len(self.variables)*self.node_id+self.variables.index(variable)]
                except IndexError:
                    data = predictions_times,predictions[self.training_timesteps:]
            if not self.plot_zero_frequency:
                data = data - np.mean(data)
            abs_fourier_transform = abs(np.fft.rfft(data))
            frequency = np.linspace(0, 1.0/self.sampling_rate/2, len(abs_fourier_transform))
            if self.mirror_negative_frequencies:
                neg_abs_fourier_transform = np.array([x for x in reversed(abs_fourier_transform)])
                abs_fourier_transform = neg_abs_fourier_transform + abs_fourier_transform
                neg_frequency = np.array([-x for x in reversed(frequency)])
                frequency = neg_frequency + frequency
            getattr(ax,self.axis)(frequency, abs_fourier_transform,'r-',label='Predictions')
            ax.set_xlabel('Frequency')
            ax.set_ylabel('FFT')
            ax.legend()
            # THIS NO LONGER WORKS WITH NEWER VERSIONS OF MATPLOTLIB, LOOKING FOR AN ALTERNATIVE
            # manager = plt.get_current_fig_manager()
            # manager.resize(*manager.window.maxsize())
            figure.show()
            figure.savefig(self.output_name + "_" + variable + "." + self.output_format, bbox_inches='tight')

            


