import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
import tensorflow.keras.metrics
from importlib import import_module


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return MetricCompilerProcess(settings["parameters"])


class MetricCompilerProcess(NeuralNetworkProcess):
    """
    This class compiles a custom loss function.

    """
    def __init__(self, parameters):
        super().__init__()

        default_settings = KM.Parameters("""{
            "help"                     : "This process compiles a metric.",
            "metric"                   : "",
            "module"                   : "keras"
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        if parameters.Has("metric"): 
            self.metric_name = parameters["metric"].GetString()
        else:
            raise Exception("No metric specified")
        self.module = parameters["module"].GetString()
        metric_module = import_module(self.module)
        if self.module == "keras":
            self.metric = getattr(tensorflow.keras.metrics,self.metric_name)()
        else:
            self.metric = getattr(metric_module,self.metric_name)()

    def CompileMetric(self):
        """ Process for compiling a metric. """
        
        return self.metric