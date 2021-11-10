import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
from KratosMultiphysics.RomApplication.structural_mechanics_analysis_rom import StructuralMechanicsAnalysisROM
import json

class TestStructuralMechanicsStaticHROM(StructuralMechanicsAnalysisROM):

    def __init__(self,model,project_parameters):
        super().__init__(model,project_parameters)

    def ModifyAfterSolverInitialize(self):
        super().ModifyAfterSolverInitialize()
        computing_model_part = self._solver.GetComputingModelPart()

        ## Adding the weights to the corresponding elements
        with open('ElementsAndWeights.json') as f:
            HR_data = json.load(f)
            for key in HR_data["Elements"].keys():
                computing_model_part.GetElement(int(key)+1).SetValue(romapp.HROM_WEIGHT, HR_data["Elements"][key])
            for key in HR_data["Conditions"].keys():
                computing_model_part.GetCondition(int(key)+1).SetValue(romapp.HROM_WEIGHT, HR_data["Conditions"][key])

    def EvaluateQuantityOfInterest(self):
       ##############################################################################################
       # Functions evaluating the QoI of the problem: Array of displacement at every node on the mesh #
       #                    and nodal area                                                           #
       ##############################################################################################
        ArrayOfDisplacements = []
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            ArrayOfDisplacements.append(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0))
            ArrayOfDisplacements.append(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0))
        return ArrayOfDisplacements

    def EvaluateQuantityOfInterest2(self):
        computing_model_part = self._solver.GetComputingModelPart()
        dimension = self._GetSolver().settings["domain_size"].GetInt()
        area_calculator = KratosMultiphysics.CalculateNodalAreaProcess(computing_model_part , dimension)
        area_calculator.Execute()

        ArrayOfAreas = []
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            ArrayOfAreas.append(node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA))
            ArrayOfAreas.append(node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA))
        return ArrayOfAreas


if __name__ == "__main__":
    model = KratosMultiphysics.Model()
    with open("ProjectParametersHROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    Simulation = TestStructuralMechanicsStaticHROM(model,parameters)
    Simulation.Run()