from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
from KratosMultiphysics.RomApplication.convection_diffusion_analysis_rom import ConvectionDiffusionAnalysisROM

class TestConvectionDiffusionStationaryROM(ConvectionDiffusionAnalysisROM):

    def __init__(self,model,project_parameters):
        super(TestConvectionDiffusionStationaryROM,self).__init__(model,project_parameters)

    def EvaluateQuantityOfInterest(self):
       ##############################################################################################
       # Functions evaluating the QoI of the problem: Array of temperature at every node on the mesh #
       #                    and nodal area                                                           #
       ##############################################################################################
        ArrayOfTemperatures = []
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            ArrayOfTemperatures.append(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        return ArrayOfTemperatures

    def EvaluateQuantityOfInterest2(self):
        computing_model_part = self.model["ThermalModelPart"]
        dimension = self._GetSolver().settings["domain_size"].GetInt()
        area_calculator = KratosMultiphysics.CalculateNodalAreaProcess(computing_model_part , dimension)
        area_calculator.Execute()        

        ArrayOfAreas = []
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            ArrayOfAreas.append(node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA))
        return ArrayOfAreas

############################################################################################################

if __name__ == "__main__":

    with open("ProjectParametersROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = TestConvectionDiffusionStationaryROM(model,parameters)
    simulation.Run()
