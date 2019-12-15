from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.ConvectionDiffusionApplication


from convection_diffusion_analysis import ConvectionDiffusionAnalysis

import sys
import time

class ConvectionDiffusionAnalysisWithFlush(ConvectionDiffusionAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=10.0):
        super(ConvectionDiffusionAnalysisWithFlush,self).__init__(model,project_parameters)

        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def FinalizeSolutionStep(self):
        super(ConvectionDiffusionAnalysisWithFlush,self).FinalizeSolutionStep()

        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = ConvectionDiffusionAnalysisWithFlush(model,parameters)
    simulation.Run()
