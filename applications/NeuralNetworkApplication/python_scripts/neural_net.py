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

class FullyConnected_2d_t(torch.nn.Module):
    def __init__(self, n_hidden, n_layer):
        super(FullyConnected_2d_t, self).__init__()
        self.hidden1 = torch.nn.Linear(3, n_hidden).double().to(device) #.to(torch.cdouble)
        torch.nn.init.xavier_uniform_(self.hidden1.weight, gain=1)
        torch.nn.init.ones_(self.hidden1.bias)
        
        self.hidden2 = torch.nn.Linear(n_hidden, n_hidden).double().to(device) #.to(torch.cdouble)   # hidden layer
        torch.nn.init.xavier_uniform_(self.hidden2.weight, gain=1)
        torch.nn.init.ones_(self.hidden2.bias)
        
        self.predict = torch.nn.Linear(n_hidden, 3).double().to(device) # .to(torch.cdouble)   # output layer
        torch.nn.init.xavier_uniform_(self.predict.weight, gain=1)
        torch.nn.init.ones_(self.predict.bias)

        self.bn1 = torch.nn.BatchNorm1d(num_features=n_hidden).double()

        self.num_layers = n_layer
        self.relu = torch.nn.ReLU()
        self.tanh = torch.nn.Tanh()
        
        self.dropout = torch.nn.Dropout(p=0.2)
    
    def forward(self, x, train=False):
        x = self.hidden1(x)
        # x = self.bn1(x)
        x = self.tanh(x)
        # x = self.dropout(x)

        # for creating number of layers dynamically
        for i in range(self.num_layers):
            x = self.hidden2(x)
            # x = self.bn1(x)
            x = self.tanh(x)
            # x = self.dropout(x)
        x = self.predict(x)
        return x