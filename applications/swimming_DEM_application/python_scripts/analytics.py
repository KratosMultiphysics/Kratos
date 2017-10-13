from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import math
import os
import numpy as np

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

import swimming_DEM_procedures
import matplotlib.pyplot as plt
import numpy.fft as fft

class Gauge:
    def __init__(self, model_part, Dt, final_time, variables_to_measure, steps_between_measurements, activity = True):
        self.model_part = model_part
        self.steps_between_measurements = steps_between_measurements
        self.Dt = Dt * steps_between_measurements
        self.final_time = final_time
        self.variables = variables_to_measure
        self.n_vars = len(variables_to_measure)
        self.n_instants = int(math.ceil(self.final_time / self.Dt + 1))
        self.current_instant = 0
        self.counter = swimming_DEM_procedures.Counter(steps_between_measurements, 1, activity)

    def ConstructArrayOfNodes(self, condition):
        self.nodes = [node for node in self.model_part.Nodes if condition(node)]
        self.n_nodes = len(self.nodes)
        zeros = [0.0 for i in range(self.n_instants)]
        zeros_for_a_node = [zeros for var in self.variables]
        self.times = np.array(zeros)
        self.measurements = [np.array(zeros_for_a_node) for node in self.nodes]

    def MakeNodalMeasurement(self):
        if self.counter.Tick():
            self.times[self.current_instant] = self.model_part.ProcessInfo[TIME]
            print(self.model_part.ProcessInfo[TIME] / self.final_time)
            print(self.current_instant / self.n_instants)
            i_node = 0
            for node in self.nodes:
                i_var = 0
                for var in self.variables:
                    self.measurements[i_node][i_var][self.current_instant] = node.GetSolutionStepValue(var)
                    i_var += 1
                i_node += 1
            self.current_instant += 1

    def PrintMeasurements(self, path):
        for i_node in range(self.n_nodes):
            with open("analytics_node_" + str(self.nodes[i_node].Id) + "_measurements.txt", 'w') as f:
                for i_instant in range(self.current_instant):
                    f.write(str(self.times[i_instant]))
                    for i_var in range(self.n_vars):
                        entry = self.measurements[i_node][i_var][i_instant]
                        f.write(" " + str(entry))
                    f.write("\n")

    def PlotPSD(self):
        colors = ['r']
        dt = self.Dt
        Ry = self.measurements[0][0]
        mean = np.average(Ry)
        val = np.array( [v - mean for v in Ry] )
        freq = fft.rfftfreq(len(val), dt)
        tr = fft.rfft(val)
        tr_vect = [np.real(i) for i in tr]

        psd = [(np.real(i)**2 + np.imag(i)**2)**0.5 for i in tr]

        plt.clf()
        plt.figure("psd")
        plt.plot(freq[1:],tr_vect[1:],colors[0])
        plt.xlim(0, 10.0)
        plt.xlabel(r"$Frequency \, (Hz)$")
        plt.ylabel(r"$PSD$")
        plt.savefig("psd.eps")
