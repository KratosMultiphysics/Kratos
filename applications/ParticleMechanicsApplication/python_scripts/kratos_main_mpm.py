
import KratosMultiphysics
from KratosMultiphysics.ParticleMechanicsApplication.particle_mechanics_analysis import ParticleMechanicsAnalysis

"""
For user-scripting it is intended that a new class is derived
from ParticleMechanicsAnalysis to do modifications
"""

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = ParticleMechanicsAnalysis(model,parameters)
    simulation.Run()
