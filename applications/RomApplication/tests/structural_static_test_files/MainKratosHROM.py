from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
from KratosMultiphysics.RomApplication.structural_mechanics_analysis_rom import StructuralMechanicsAnalysisROM
import json

class TestStructuralMechanicsStaticHROM(StructuralMechanicsAnalysisROM):

    def __init__(self,model,project_parameters):
        super(TestStructuralMechanicsStaticHROM,self).__init__(model,project_parameters)

    def ModifyInitialGeometry(self):
        super(TestStructuralMechanicsStaticHROM,self).ModifyInitialGeometry()
        computing_model_part = self._solver.GetComputingModelPart()

        ## Adding the weights to the corresponding elements
        with open('ElementsAndWeights.json') as f:
            HR_data = json.load(f)
            for key in HR_data["Elements"].keys():
                computing_model_part.GetElement(int(key)+1).SetValue(romapp.HROM_WEIGHT, HR_data["Elements"][key])
            for key in HR_data["Conditions"].keys():
                computing_model_part.GetCondition(int(key)+1).SetValue(romapp.HROM_WEIGHT, HR_data["Conditions"][key])



if __name__ == "__main__":
    model = KratosMultiphysics.Model()
    with open("ProjectParametersHROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    Simulation = TestStructuralMechanicsStaticHROM(model,parameters)
    Simulation.Run()