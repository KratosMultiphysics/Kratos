import tensorflow.keras as keras
import numpy as np
import KratosMultiphysics
import KratosMultiphysics.NeuralNetworkApplication as NeuralNetworkApplication
from KratosMultiphysics.NeuralNetworkApplication.neural_network_training_process import NeuralNetworkTrainingProcess
from importlib import import_module
from kerastuner.tuners import RandomSearch

def Factory(settings):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return TunerFromSavedModel(settings["parameters"])

class TunerFromSavedModel(NeuralNetworkTrainingProcess):
    """ 
    This class creates a process that trains a neural network located at a saved model.
    """
    def __init__(self, parameters):
        super().__init__(parameters)
        self.model_name = parameters["model"].GetString()
        self.epochs = parameters["epochs"].GetInt()
        self.num_models = parameters["num_models"].GetInt()
        self.objective = parameters["objective"].GetString()
        self.max_trials = parameters["max_trials"].GetInt()
        self.executions_per_trial = parameters["executions_per_trial"].GetInt()
        self.directory = parameters["directory"].GetString()
        self.project_name = parameters["project_name"].GetString()

    def ExecuteTraining(self, loss_function, optimizer, callbacks_list, metrics_list):
        """ Processes to act directly during the training step. """
        model_import = import_module(self.model_name)
        self.hypermodel = getattr(model_import, "build_model")
        self.tuner = RandomSearch(self.hypermodel, loss = loss_function,
                                    optimizer = optimizer, objective = self.objective,
                                    max_trials=self.max_trials, executions_per_trial = self.executions_per_trial,
                                    directory = self.directory, project_name = self.project_name)
        self.tuner.search_space_summary()
        self.tuner.search(self.test_input,self.test_output,epochs=self.epochs)
        models = self.tuner.get_best_models(num_models=self.num_models)
        self.tuner.results_summary()
        for i in range(self.num_models):
            models[i].save(self.model_name + "_best_"+str(i))
        