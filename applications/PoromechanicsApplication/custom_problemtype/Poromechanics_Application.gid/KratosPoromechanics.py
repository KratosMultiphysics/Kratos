from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
# import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.PoromechanicsApplication

from KratosMultiphysics.PoromechanicsApplication.poromechanics_analysis import PoromechanicsAnalysis

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = PoromechanicsAnalysis(model,parameters)
    simulation.Run()
