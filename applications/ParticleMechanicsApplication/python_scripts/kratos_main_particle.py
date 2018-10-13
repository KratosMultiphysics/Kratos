from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.ParticleMechanicsApplication

from particle_mechanics_analysis import ParticleMechanicsAnalysis

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
