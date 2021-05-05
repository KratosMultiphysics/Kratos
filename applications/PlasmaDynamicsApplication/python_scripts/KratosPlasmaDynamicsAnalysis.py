import KratosMultiphysics
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
# from KratosMultiphysics.FluidTransportApplication import *
from KratosMultiphysics.PlasmaDynamicsApplication import *

from KratosMultiphysics.PlasmaDynamicsApplication.plasma_dynamics_analysis import PlasmaDynamicsAnalysis

class PlasmaDynamicsAnalysisWithFlush(PlasmaDynamicsAnalysis):
    
    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)


if __name__ == "__main__":
    with open("ProjectParameters.json",'r') as parameter_file:
        project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = PlasmaDynamicsAnalysisWithFlush(model, project_parameters)
    simulation.Run()
