# Import necessary modules from KratosMultiphysics
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.modelers.meshio_input_modeler import Factory as MeshioInputModelerFactory

# Import the pathlib and tempfile modules
from pathlib import Path
import tempfile


class TestMeshioInputModeler(KratosUnittest.TestCase):
    """Tests for the meshio++ input modeler."""

    def setUp(self):
        self.model = KratosMultiphysics.Model()

    def _WriteInputFile(self, file_name):
        write_model_part = self.model.CreateModelPart("to_write")
        props = write_model_part.CreateNewProperties(1)
        write_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        write_model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        write_model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        write_model_part.CreateNewNode(4, 0.0, 0.0, 1.0)
        write_model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], props)
        settings = KratosMultiphysics.Parameters("""{"time_series" : "single_file"}""")
        KratosMultiphysics.MeshioPlusPlusIO(file_name, settings).WriteModelPart(write_model_part)

    def testImportVtu(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            file_name = str(Path(temp_dir) / "input_mesh.vtu")
            self._WriteInputFile(file_name)

            settings = KratosMultiphysics.Parameters("""{
                "input_filename"  : "%s",
                "model_part_name" : "imported"
            }""" % file_name.replace("\\\\", "/"))
            modeler = MeshioInputModelerFactory(self.model, settings)

            # The destination model part must exist right after construction so
            # solvers can add nodal variables before the mesh is read
            self.assertTrue(self.model.HasModelPart("imported"))
            self.assertEqual(self.model["imported"].NumberOfNodes(), 0)

            # Run the modeler lifecycle (as AnalysisStage does)
            modeler.SetupGeometryModel()
            modeler.PrepareGeometryModel()
            modeler.SetupModelPart()

            imported_model_part = self.model["imported"]
            self.assertEqual(imported_model_part.NumberOfNodes(), 4)
            self.assertEqual(imported_model_part.NumberOfElements(), 1)

    def testMissingModelPartNameRaises(self):
        settings = KratosMultiphysics.Parameters("""{
            "input_filename" : "whatever.vtu"
        }""")
        with self.assertRaisesRegex(Exception, "model_part_name"):
            MeshioInputModelerFactory(self.model, settings)

    def testModelerIsRegistered(self):
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Modelers.KratosMultiphysics.MeshioInputModeler"))


if __name__ == '__main__':
    KratosUnittest.main()
