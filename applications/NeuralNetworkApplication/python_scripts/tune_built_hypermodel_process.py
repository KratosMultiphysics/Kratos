import tensorflow.keras as keras
import numpy as np
import KratosMultiphysics
import KratosMultiphysics.NeuralNetworkApplication as NeuralNetworkApplication
from KratosMultiphysics.NeuralNetworkApplication.neural_network_training_process import NeuralNetworkTrainingProcess
from importlib import import_module
import kerastuner.tuners as tuners

def Factory(settings):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return TunerFromBuiltModel(settings["parameters"])

class TunerFromBuiltModel(NeuralNetworkTrainingProcess):
    """ 
    This class creates a process that trains a neural network located at a saved model.
    """
    def __init__(self, parameters):
        super().__init__(parameters)
        self.model_name = parameters["model_name"].GetString()
        self.tuner = parameters["tuner"].GetString()
        self.epochs = parameters["epochs"].GetInt()
        self.num_models = parameters["num_models"].GetInt()
        self.objective = parameters["objective"].GetString()
        self.max_trials = parameters["max_trials"].GetInt()
        self.executions_per_trial = parameters["executions_per_trial"].GetInt()
        self.directory = parameters["directory"].GetString()
        self.project_name = parameters["project_name"].GetString()

    def ExecuteTuning(self, hypermodel):
        """ Processes to act directly during the training step. """
        self.tuner = getattr(tuners, self.tuner)(hypermodel, objective = self.objective,
                                    max_trials=self.max_trials, executions_per_trial = self.executions_per_trial,
                                    directory = self.directory, project_name = self.project_name)
        self.tuner.search_space_summary()
        self.tuner.search(self.test_input.ExportAsArray(),self.test_output.ExportAsArray(),epochs=self.epochs)
        
    
    def ExecuteFinalize(self):
        
        models = self.tuner.get_best_models(num_models=self.num_models)
        self.tuner.results_summary()
        for i in range(self.num_models):
            models[i].save(self.model_name + "_best_"+str(i))
        