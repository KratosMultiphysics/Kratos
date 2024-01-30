
import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess
from test_vtk_output_process import SetupModelPart2D, SetupModelPart3D

class TestVtuOutputBase:
    @classmethod
    def SetSolution(cls):
        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, 0,[node.X * 2, node.Y * 3, node.Z * 4])
            node.SetSolutionStepValue(Kratos.VELOCITY, 0,[node.X * 5, node.Y * 6, node.Z * 7])
            node.SetSolutionStepValue(Kratos.PRESSURE, 0, node.Id * 8.0)

        for i_elem, elem in enumerate(cls.model_part.Elements):
            elem.SetValue(Kratos.DETERMINANT, [i_elem*0.189, i_elem * 1.236, i_elem * 2.365])

        for i_cond, cond in enumerate(cls.model_part.Conditions):
            cond.SetValue(Kratos.DENSITY, i_cond * 4.362)
            cond.SetValue(Kratos.YOUNG_MODULUS, i_cond * 5.326)

    @classmethod
    def setUpClass(cls, output_prefix: str, setup_method) -> None:
        cls.model =  Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.output_prefix = output_prefix
        setup_method(cls.model_part)
        cls.SetSolution()

    def WriteVtu(self, output_format: Kratos.VtuOutput.WriterFormat):
        vtu_output = Kratos.VtuOutput(self.model_part, True, output_format, 9)
        vtu_output.AddHistoricalVariable(Kratos.PRESSURE)
        vtu_output.AddHistoricalVariable(Kratos.VELOCITY)
        vtu_output.AddHistoricalVariable(Kratos.DISPLACEMENT)
        vtu_output.AddNonHistoricalVariable(Kratos.DETERMINANT, vtu_output.ELEMENTS)

        a = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(a, Kratos.PRESSURE, True)
        a *= 3
        vtu_output.AddContainerExpression("hist_exp", a)

        a = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(a, Kratos.DETERMINANT)
        a *= 3
        vtu_output.AddContainerExpression("elem_exp", a)

        with kratos_unittest.WorkFolderScope("./auxiliar_files_for_python_unittest/vtk_output_process_ref_files", __file__, True):
            if output_format == Kratos.VtuOutput.ASCII:
                output_file_prefix = "ascii" + self.output_prefix + "/Main"
            else:
                output_file_prefix = "binary" + self.output_prefix + "/Main"
            vtu_output.PrintOutput(output_file_prefix + "_temp")
            self.Check(output_file_prefix + "_temp.vtu",  output_file_prefix + ".vtu")

    def test_WriteMeshAscii(self):
        self.WriteVtu(Kratos.VtuOutput.ASCII)

    def test_WriteMeshBinary(self):
        self.WriteVtu(Kratos.VtuOutput.BINARY)

    def Check(self, output_file, reference_file):
        ## Settings string in json format
        params = Kratos.Parameters("""{
            "reference_file_name" : "",
            "output_file_name"    : "",
            "comparison_type"     : "deterministic"
        }""")
        params["reference_file_name"].SetString(reference_file)
        params["output_file_name"].SetString(output_file)
        CompareTwoFilesCheckProcess(params).Execute()

class TestVtuOutput2D(TestVtuOutputBase, kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass("2D", SetupModelPart2D)

class TestVtuOutput3D(TestVtuOutputBase, kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass("3D", SetupModelPart3D)

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
    kratos_unittest.main()