
import KratosMultiphysics as Kratos
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestVtuOutput(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model =  Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        with kratos_unittest.WorkFolderScope(".", __file__, True):
            ReadModelPart("auxiliar_files_for_python_unittest/mdpa_files/two_dim_symmetrical_square", cls.model_part)

        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.Id)
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([node.Id, node.Id + 1, node.Id + 2]))

    def test_WriteMeshAscii(self):
        vtu_output = Kratos.VtuOutput(self.model_part, True, Kratos.VtuOutput.ASCII, 9)
        vtu_output.AddHistoricalVariable(Kratos.PRESSURE)
        vtu_output.AddHistoricalVariable(Kratos.VELOCITY)

        a = Kratos.ContainerExpression.HistoricalExpression(self.model_part)
        a.Read(Kratos.PRESSURE)
        a *= 3

        vtu_output.AddContainerExpression("my_a", a)
        a.Read(Kratos.VELOCITY)
        a *= 100.0

        vtu_output.PrintOutput("test_1_ascii")

    def test_WriteMeshBinary(self):
        vtu_output = Kratos.VtuOutput(self.model_part, True, Kratos.VtuOutput.BINARY, 9)
        vtu_output.AddHistoricalVariable(Kratos.PRESSURE)
        vtu_output.AddHistoricalVariable(Kratos.VELOCITY)

        a = Kratos.ContainerExpression.HistoricalExpression(self.model_part)
        a.Read(Kratos.PRESSURE)
        a *= 3

        vtu_output.AddContainerExpression("my_a", a)
        a.Read(Kratos.VELOCITY)
        a *= 100.0

        vtu_output.PrintOutput("test_1_binary")

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()