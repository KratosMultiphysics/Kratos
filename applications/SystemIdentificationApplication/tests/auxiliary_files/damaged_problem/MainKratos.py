import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as KratosSM
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

if __name__ == "__main__":
    model = KratosMultiphysics.Model()
    with open("auxiliary_files/damaged_problem/damaged_project_parameters.json", "r") as file_input:
        params = KratosMultiphysics.Parameters(file_input.read())

    analysis = StructuralMechanicsAnalysis(model, params)
    analysis.Run()