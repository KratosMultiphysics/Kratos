import KratosMultiphysics
from KratosMultiphysics.RomApplication.fluid_dynamics_analysis_rom import FluidDynamicsAnalysisROM

import numpy as np

class TestFluidDynamicsROM(FluidDynamicsAnalysisROM):

    def __init__(self,model,project_parameters):
        super().__init__(model,project_parameters)

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()
    simulation = TestFluidDynamicsROM(model,parameters)
    simulation.Run()
