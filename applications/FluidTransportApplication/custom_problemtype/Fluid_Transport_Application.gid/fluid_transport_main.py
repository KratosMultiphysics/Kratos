from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.ConvectionDiffusionApplication
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.FluidTransportApplication

from fluid_transport_analysis import FluidTransportAnalysis

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = FluidTransportAnalysis(model,parameters)
    simulation.Run()