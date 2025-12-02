# --- Kratos Imports ---
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.KratosUnittest import WorkFolderScope

# --- STD Imports ---
import pathlib


class TestPretension(KratosMultiphysics.KratosUnittest.TestCase):

    def test_pretension_1d(self) -> None:
        with WorkFolderScope("constraints", pathlib.Path(__file__).absolute()):
            # Load config
            with open("pretension_1d.json", "r") as project_parameters_file:
                parameters = KratosMultiphysics.Parameters(project_parameters_file.read())
            model = KratosMultiphysics.Model()
            analysis = StructuralMechanicsAnalysis(model, parameters)
            analysis.Run()


if __name__ == "__main__":
    KratosMultiphysics.KratosUnittest.main()
