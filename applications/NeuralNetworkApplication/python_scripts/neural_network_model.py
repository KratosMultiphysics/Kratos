import tensorflow.keras as keras
import numpy as np
import torch
from KratosMultiphysics.NeuralNetworkApplication.neural_net import Fullyconnected, ConstructModel

class MachineLearningModel():
    """
    This class deals with the machine learning models constructed using different ML libraries like Keras, Pytorch, Julia etc
    """
    def __init__(self, param):  

        self.paramaters = param        
        self.library_name = self.paramaters["library_name"].GetString()    
        # self.net_arch_filename = self.paramaters["net_arch_filename"].GetString()    
        self.model = None 
        
    def load_model(self, file_name):

        if self.library_name=="keras":
            self.model = keras.models.load_model(file_name)
            self.model.compile()
            self.model.summary()

        elif self.library_name=="pytorch":
            device = torch.device('cpu')
            model_constructor = ConstructModel(self.net_arch_filename)
            self.model = model_constructor.get_model()
            self.model.load_state_dict(torch.load(file_name, map_location=device))
        else: 
            print("Please use either keras or pytorch model")
        
        return self.model

    def create_model(self, inputs, outputs):
        if self.library_name=="keras":
            self.model = keras.Model(inputs = inputs, outputs = outputs)
            self.model.summary()
        return self.model


    def summary(self):
        if self.library_name=="keras":
            print(self.model.summary())

    def get_model(self):
        return self.model
    
    def get_model_library(self):
        return self.library_name
