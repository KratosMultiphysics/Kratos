import KratosMultiphysics as Kratos
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

class CustomStructuralMechanicsAnalysis(StructuralMechanicsAnalysis):
    def FinalizeSolutionStep(self):
        for element in self._GetSolver().GetComputingModelPart().Elements:
            element.SetValue(Kratos.YOUNG_MODULUS, element.Properties[Kratos.YOUNG_MODULUS])

if __name__ == "__main__":
    model = Kratos.Model()

    with open("primal_project_parameters.json", "r") as file_input:
        parameters = Kratos.Parameters(file_input.read())

    analysis = CustomStructuralMechanicsAnalysis(model, parameters)
    analysis.Run()