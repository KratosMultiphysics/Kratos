from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import math
import os

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

import swimming_DEM_procedures

class Gauge:
    def __init__(self, model_part):
        self.model_part = model_part
        self.measurements = []

    def ConstructArrayOfNodes(self, condition):
        self.nodes = [node for node in self.model_part.Nodes if condition(node)]
        self.coors = [[node.X, node.Y, node.Z] for node in self.nodes]
        self.number_of_nodes = len(self.nodes)
        print(self.nodes)
        print(self.coors)


    def SetListOfNodalScalarVariables(self, variables_to_measure):
        self.variables = variables_to_measure
        self.n_variables = len(self.variables)

    def MakeNodalMeasurement():
        measurement = [[node.GetSolutionStepValue(var) for var in self.variables] for node in self.nodes]
        self.measurements += measurement
        print(measurement)

#import matplotlib.pyplot as plt
#from numpy import linspace
#import numpy as np
#import numpy.fft as fft

#dt = T[1]-T[0]

#mean = np.average(Ry)
#val = np.array( [v - mean for v in Ry] )

#freq = fft.rfftfreq(len(val),dt)
#tr = fft.rfft(val)

#psd = [ (np.real(i)**2 + np.imag(i)**2)**0.5 for i in tr ]

#plt.figure("psd")
#plt.plot(freq,psd,colors[0])
#plt.xlim(0,1.0)
#plt.xlabel(r"$Frequency \, (Hz)$")
#plt.ylabel(r"$PSD$")
#plt.savefig("psd.eps")
