from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
# import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.GeoMechanicsApplication

from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import GeoMechanicsAnalysis

if __name__ == "__main__":

    model = KratosMultiphysics.Model()
    with open("ProjectParameters_stage1.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
        stage_1 = GeoMechanicsAnalysis(model, parameters)

    with open("ProjectParameters_stage2.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
        stage_2 = GeoMechanicsAnalysis(model, parameters)


    list_of_stages = [stage_1, stage_2]

    for stage in list_of_stages:
        stage.Run()
