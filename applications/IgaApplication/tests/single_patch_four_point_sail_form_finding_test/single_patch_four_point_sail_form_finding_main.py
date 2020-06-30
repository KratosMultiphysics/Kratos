from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import KratosMultiphysics.IgaApplication

class ModifiedStructuralMechanicsAnalysis(StructuralMechanicsAnalysis):
    def Initialize(self):
        super(ModifiedStructuralMechanicsAnalysis, self).Initialize()

if __name__ == "__main__":
    with open("single_patch_four_point_sail_form_finding_project_parameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = ModifiedStructuralMechanicsAnalysis(model,parameters)
    simulation.Run()
    print("KRATOS TERMINATED")
