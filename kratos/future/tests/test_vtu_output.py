from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

# import from the tests
from kratos.tests.test_vtk_output_process import SetupModelPart2D, SetupModelPart3D

class TestVtuOutputBase:
    @classmethod
    def SetSolution(cls):
        node: Kratos.Node
        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, 0,[node.X * 2, node.Y * 3, node.Z * 4])
            node.SetSolutionStepValue(Kratos.VELOCITY, 0,[node.X * 5, node.Y * 6, node.Z * 7])
            node.SetSolutionStepValue(Kratos.PRESSURE, 0, node.Id * 8.0)

        elem: Kratos.Element
        for i_elem, elem in enumerate(cls.model_part.Elements):
            elem.SetValue(Kratos.DETERMINANT, [i_elem*0.189, i_elem * 1.236, i_elem * 2.365])

        cond: Kratos.Condition
        for i_cond, cond in enumerate(cls.model_part.Conditions):
            cond.SetValue(Kratos.DENSITY, i_cond * 4.362)
            cond.SetValue(Kratos.YOUNG_MODULUS, i_cond * 5.326)

    @classmethod
    def setUpClass(cls, output_prefix: str, setup_method, output_sub_model_parts: bool) -> None:
        cls.model =  Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.STEP] = 0
        cls.model_part.ProcessInfo[Kratos.TIME] = 1.0
        cls.output_prefix = output_prefix
        cls.output_sub_model_parts = output_sub_model_parts
        setup_method(cls.model_part)
        cls.SetSolution()

    def WriteVtu(self, output_format: Kratos.Future.VtuOutput.WriterFormat):
        vtu_output = Kratos.Future.VtuOutput(self.model_part, Kratos.Configuration.Initial, output_format, 9, echo_level=0, output_sub_model_parts=self.output_sub_model_parts)
        vtu_output.AddVariable(Kratos.PRESSURE, Kratos.Globals.DataLocation.NodeHistorical)
        vtu_output.AddVariable(Kratos.VELOCITY, Kratos.Globals.DataLocation.NodeHistorical)
        vtu_output.AddVariable(Kratos.DISPLACEMENT, Kratos.Globals.DataLocation.NodeHistorical)
        vtu_output.AddVariable(Kratos.DETERMINANT, Kratos.Globals.DataLocation.Element)
        vtu_output.AddVariable(Kratos.DENSITY, Kratos.Globals.DataLocation.Condition)
        vtu_output.AddVariable(Kratos.YOUNG_MODULUS, Kratos.Globals.DataLocation.Condition)

        ta_1 = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
        ta_1.CollectData()
        ta_2 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.DETERMINANT)
        ta_2.CollectData()

        ta_1.data *= 3
        ta_2.data *= 3

        vtu_output.AddTensorAdaptor("hist_ta", ta_1)
        vtu_output.AddTensorAdaptor("elem_ta", ta_2)

        with kratos_unittest.WorkFolderScope("vtk_output_process_ref_files", __file__, True):
            output_file_prefix = output_format.name.lower() + self.output_prefix + "/Main"
            vtu_output.PrintOutput("temp/" + output_file_prefix)
            self.Check("temp/" + output_file_prefix,  output_file_prefix)

    def test_WriteMeshAscii(self):
        self.WriteVtu(Kratos.Future.VtuOutput.ASCII)

    def test_WriteMeshBinary(self):
        self.WriteVtu(Kratos.Future.VtuOutput.BINARY)

    def test_WriteMeshRaw(self):
        self.WriteVtu(Kratos.Future.VtuOutput.RAW)

    def test_WriteMeshCompressedRaw(self):
        self.WriteVtu(Kratos.Future.VtuOutput.COMPRESSED_RAW)

    def Check(self, output_prefix, reference_prefix):
        def check_file(output_file_name: str, reference_file_name: str):
            ## Settings string in json format
            params = Kratos.Parameters("""{
                "reference_file_name" : "",
                "output_file_name"    : "",
                "comparison_type"     : "deterministic"
            }""")
            params["reference_file_name"].SetString(reference_file_name)
            params["output_file_name"].SetString(output_file_name)
            CompareTwoFilesCheckProcess(params).Execute()

        for file_path in Path(reference_prefix).iterdir():
            self.assertTrue((Path(output_prefix) / file_path.name).is_file(), msg=f"\"{(Path(output_prefix) / file_path.name)}\" is not a file.")
            check_file(f"{output_prefix}/{file_path.name}", str(file_path))
        check_file(f"{output_prefix}.pvd", f"{reference_prefix}.pvd")

        kratos_utils.DeleteDirectoryIfExistingAndEmpty("temp")

class TestVtuOutput2D(TestVtuOutputBase, kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass("2D", SetupModelPart2D, output_sub_model_parts = True)

class TestVtuOutput3D(TestVtuOutputBase, kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        # the method SetupModelPart3D does not create sub model parts
        # with nodes which do not include nodes from its conditions or elements. It uses
        # some random nodes. Hence sub_model_part output is disabled.
        super().setUpClass("3D", SetupModelPart3D, output_sub_model_parts = False)

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.INFO)
    kratos_unittest.main()