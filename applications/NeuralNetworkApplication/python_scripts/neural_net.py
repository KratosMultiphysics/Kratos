import torch
import torch.nn as nn
import torch.nn.functional as F
import sys 
device = 'cuda' if torch.cuda.is_available() else 'cpu'
"""

Contains neural network models for Pytorch, only by loading models from this pytorch learned paramaters can be used
"""
# from neural_net.layer_init import layer_initialisation

# torch.set_default_dtype(torch.DoubleTensor) 
# this is one way to define a network

class Fullyconnected(torch.nn.Module):
    """
    
    A class which defines model architecture with ´n_layer´ number of hidden 
    layers.

    Xavier uniform initialization of weights. ReLu activation function. 
    
    """
    
    def __init__(self, n_hidden, n_layer):
        super(Fullyconnected, self).__init__()
        self.hidden1 = torch.nn.Linear(3, n_hidden).double().to(device)   # hidden layer
        nn.init.xavier_uniform_(self.hidden1.weight,gain=1)
        nn.init.ones_(self.hidden1.bias)
        
        self.hidden2 = torch.nn.Linear(n_hidden, n_hidden).double().to(device)   # hidden layer
        nn.init.xavier_uniform_(self.hidden2.weight,gain=1)
        nn.init.ones_(self.hidden2.bias)
        
        self.predict = torch.nn.Linear(n_hidden, 1).double().to(device)   # output layer
        nn.init.xavier_uniform_(self.predict.weight,gain=1)
        nn.init.ones_(self.predict.bias)

        self.num_layers = n_layer

        self.tanh = torch.nn.Tanh()
    
    def forward(self, x):
        x = self.hidden1(x)
        x = self.tanh(x)   

        # for creating number of layers dynamically
        for i in range(self.num_layers):
            x = self.hidden2(x) 
            x = self.tanh(x)   
        x = self.predict(x)
        return x

