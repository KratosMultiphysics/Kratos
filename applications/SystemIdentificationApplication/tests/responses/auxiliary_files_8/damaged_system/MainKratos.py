import KratosMultiphysics as Kratos
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

class CustomStructuralMechanicsAnalysis(StructuralMechanicsAnalysis):
    def FinalizeSolutionStep(self):
        for element in self._GetSolver().GetComputingModelPart().Elements:
            element.SetValue(Kratos.YOUNG_MODULUS, element.Properties[Kratos.YOUNG_MODULUS])

if __name__ == "__main__":
    model1 = Kratos.Model()

    with open("primal_project_parameters_strain.json", "r") as file_input:
        parameters = Kratos.Parameters(file_input.read())

    analysis1 = CustomStructuralMechanicsAnalysis(model1, parameters)
    analysis1.Run()

    model2 = Kratos.Model()

    with open("primal_project_parameters_temp.json", "r") as file_input2:
        parameters2 = Kratos.Parameters(file_input2.read())

    analysis2 = CustomStructuralMechanicsAnalysis(model2, parameters2)
    analysis2.Run()