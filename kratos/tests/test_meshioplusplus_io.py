# Import necessary modules from KratosMultiphysics
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

# Import the pathlib and tempfile modules
from pathlib import Path
import tempfile

def _PopulateModelPart(model_part):
    """Creates 5 nodes, 2 tetrahedra elements and 1 triangle condition."""
    props = model_part.CreateNewProperties(1)
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
    model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
    model_part.CreateNewNode(4, 0.0, 0.0, 1.0)
    model_part.CreateNewNode(5, 1.0, 1.0, 1.0)
    model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], props)
    model_part.CreateNewElement("Element3D4N", 2, [2, 3, 4, 5], props)
    model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1, 2, 3], props)

class TestMeshioPlusPlusIO(KratosUnittest.TestCase):
    """Tests for the meshio++-based multi-format IO."""

    def setUp(self):
        self.model = KratosMultiphysics.Model()

    def testFormatIntrospection(self):
        supported = KratosMultiphysics.MeshioPlusPlusIO.GetSupportedFormats()
        for name in ("vtu", "vtk", "gmsh", "stl", "obj", "xdmf", "ensight", "vtp", "triangle"):
            self.assertIn(name, supported)

        read_formats = KratosMultiphysics.MeshioPlusPlusIO.GetSupportedReadFormats()
        write_formats = KratosMultiphysics.MeshioPlusPlusIO.GetSupportedWriteFormats()
        self.assertIn("openfoam", read_formats)   # read-only format
        self.assertNotIn("openfoam", write_formats)
        for name in ("svg", "tikz"):              # write-only formats
            self.assertIn(name, write_formats)
            self.assertNotIn(name, read_formats)

    def testFormatEnum(self):
        Format = KratosMultiphysics.MeshioPlusPlusIO.Format
        self.assertEqual(KratosMultiphysics.MeshioPlusPlusIO.FormatFromString("vtu"), Format.VTU)
        self.assertEqual(KratosMultiphysics.MeshioPlusPlusIO.FormatFromString("auto"), Format.AUTOMATIC)
        self.assertEqual(KratosMultiphysics.MeshioPlusPlusIO.FormatFromString("ensight"), Format.ENSIGHT)
        self.assertEqual(KratosMultiphysics.MeshioPlusPlusIO.FormatFromString("vtp"), Format.VTP)
        self.assertEqual(KratosMultiphysics.MeshioPlusPlusIO.FormatFromString("triangle"), Format.TRIANGLE)
        self.assertEqual(KratosMultiphysics.MeshioPlusPlusIO.FormatName(Format.GMSH), "gmsh")
        self.assertEqual(KratosMultiphysics.MeshioPlusPlusIO.FormatName(Format.SVG), "svg")
        self.assertEqual(KratosMultiphysics.MeshioPlusPlusIO.FormatName(Format.TIKZ), "tikz")
        self.assertEqual(KratosMultiphysics.MeshioPlusPlusIO.ResolveFormat("some_mesh.vtu"), Format.VTU)
        self.assertEqual(KratosMultiphysics.MeshioPlusPlusIO.ResolveFormat("some_mesh.case"), Format.ENSIGHT)
        self.assertEqual(KratosMultiphysics.MeshioPlusPlusIO.ResolveFormat("some_mesh.vtp"), Format.VTP)
        self.assertTrue(KratosMultiphysics.MeshioPlusPlusIO.IsFormatAvailable(Format.VTU))
        with self.assertRaisesRegex(RuntimeError, "Unknown format"):
            KratosMultiphysics.MeshioPlusPlusIO.FormatFromString("not_a_format")

    def testWriteReadRoundTripVtu(self):
        self._RunWriteReadRoundTrip(".vtu")

    def testWriteReadRoundTripGmsh(self):
        self._RunWriteReadRoundTrip(".msh")

    def testWriteReadRoundTripVtp(self):
        """vtp is a surface (PolyData) format: round trip a triangle mesh."""
        write_model_part = self.model.CreateModelPart("write_vtp")
        read_model_part = self.model.CreateModelPart("read_vtp")
        props = write_model_part.CreateNewProperties(1)
        write_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        write_model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        write_model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        write_model_part.CreateNewNode(4, 0.0, 1.0, 0.0)
        write_model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], props)
        write_model_part.CreateNewElement("Element2D3N", 2, [1, 3, 4], props)

        with tempfile.TemporaryDirectory() as temp_dir:
            file_name = str(Path(temp_dir) / "round_trip.vtp")
            write_settings = KratosMultiphysics.Parameters("""{"time_series" : "single_file"}""")
            KratosMultiphysics.MeshioPlusPlusIO(file_name, write_settings).WriteModelPart(write_model_part)
            KratosMultiphysics.MeshioPlusPlusIO(file_name).ReadModelPart(read_model_part)

        self.assertEqual(write_model_part.NumberOfNodes(), read_model_part.NumberOfNodes())
        self.assertEqual(write_model_part.NumberOfElements(), read_model_part.NumberOfElements())

    def testWriteReadRoundTripEnsight(self):
        """EnSight writes a .case file plus a sibling .geo file."""
        write_model_part = self.model.CreateModelPart("write_ensight")
        read_model_part = self.model.CreateModelPart("read_ensight")
        _PopulateModelPart(write_model_part)

        with tempfile.TemporaryDirectory() as temp_dir:
            file_name = str(Path(temp_dir) / "round_trip.case")
            write_settings = KratosMultiphysics.Parameters("""{"time_series" : "single_file"}""")
            KratosMultiphysics.MeshioPlusPlusIO(file_name, write_settings).WriteModelPart(write_model_part)
            self.assertTrue((Path(temp_dir) / "round_trip.geo").exists())
            KratosMultiphysics.MeshioPlusPlusIO(file_name).ReadModelPart(read_model_part)

        self.assertEqual(write_model_part.NumberOfNodes(), read_model_part.NumberOfNodes())
        self.assertEqual(write_model_part.NumberOfElements(), read_model_part.NumberOfElements())
        self.assertEqual(write_model_part.NumberOfConditions(), read_model_part.NumberOfConditions())

    def testStlSkinOfVolumeMesh(self):
        """STL writes the boundary skin of a volume mesh by default ("skin" opts out)."""
        model_part = self.model.CreateModelPart("skin_write")
        _PopulateModelPart(model_part)  # 2 tetrahedra sharing a face: 6 skin triangles

        with tempfile.TemporaryDirectory() as temp_dir:
            file_name = str(Path(temp_dir) / "skin.stl")

            settings = KratosMultiphysics.Parameters("""{"time_series" : "single_file"}""")
            KratosMultiphysics.MeshioPlusPlusIO(file_name, settings).WriteModelPart(model_part)
            read_skin = self.model.CreateModelPart("skin_read")
            KratosMultiphysics.MeshioPlusPlusIO(file_name).ReadModelPart(read_skin)
            self.assertEqual(read_skin.NumberOfElements(), 6)

            settings = KratosMultiphysics.Parameters("""{"time_series" : "single_file", "skin" : false}""")
            KratosMultiphysics.MeshioPlusPlusIO(file_name, settings).WriteModelPart(model_part)
            read_no_skin = self.model.CreateModelPart("no_skin_read")
            KratosMultiphysics.MeshioPlusPlusIO(file_name).ReadModelPart(read_no_skin)
            self.assertEqual(read_no_skin.NumberOfElements(), 1)  # only the existing triangle condition

    def testWriteReadRoundTripMed(self):
        Format = KratosMultiphysics.MeshioPlusPlusIO.Format
        if not KratosMultiphysics.MeshioPlusPlusIO.IsFormatAvailable(Format.MED):
            self.skipTest("med requires an HDF5-enabled build")
        self._RunWriteReadRoundTrip(".med")

    def _RunWriteReadRoundTrip(self, extension):
        write_model_part = self.model.CreateModelPart("write" + extension.replace(".", "_"))
        read_model_part = self.model.CreateModelPart("read" + extension.replace(".", "_"))
        _PopulateModelPart(write_model_part)

        with tempfile.TemporaryDirectory() as temp_dir:
            file_name = str(Path(temp_dir) / ("round_trip" + extension))
            write_settings = KratosMultiphysics.Parameters("""{"time_series" : "single_file"}""")
            KratosMultiphysics.MeshioPlusPlusIO(file_name, write_settings).WriteModelPart(write_model_part)
            KratosMultiphysics.MeshioPlusPlusIO(file_name).ReadModelPart(read_model_part)

        self.assertEqual(write_model_part.NumberOfNodes(), read_model_part.NumberOfNodes())
        self.assertEqual(write_model_part.NumberOfElements(), read_model_part.NumberOfElements())
        self.assertEqual(write_model_part.NumberOfConditions(), read_model_part.NumberOfConditions())
        for node_write, node_read in zip(write_model_part.Nodes, read_model_part.Nodes):
            self.assertAlmostEqual(node_write.X, node_read.X, 12)
            self.assertAlmostEqual(node_write.Y, node_read.Y, 12)
            self.assertAlmostEqual(node_write.Z, node_read.Z, 12)

    def testXdmfTimeSeries(self):
        model_part = self.model.CreateModelPart("transient")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        _PopulateModelPart(model_part)

        settings = KratosMultiphysics.Parameters("""{
            "output_control_type"                : "time",
            "xdmf_data_format"                   : "XML",
            "nodal_solution_step_data_variables" : ["TEMPERATURE"]
        }""")

        with tempfile.TemporaryDirectory() as temp_dir:
            file_name = str(Path(temp_dir) / "transient.xdmf")
            meshio_io = KratosMultiphysics.MeshioPlusPlusIO(file_name, settings.Clone())
            for step in range(1, 4):
                model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.1 * step
                for node in model_part.Nodes:
                    node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, step * node.Id)
                meshio_io.WriteModelPart(model_part)

            with open(file_name) as xdmf_file:
                content = xdmf_file.read()
            self.assertEqual(content.count("<Time Value="), 3)
            self.assertEqual(content.count('Name="mesh"'), 1)
            self.assertEqual(content.count('Name="TEMPERATURE"'), 3)

            # A new IO on the same file extends the existing time series
            appending_io = KratosMultiphysics.MeshioPlusPlusIO(file_name, settings.Clone())
            model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.4
            appending_io.WriteModelPart(model_part)

            with open(file_name) as xdmf_file:
                content = xdmf_file.read()
            self.assertEqual(content.count("<Time Value="), 4)
            self.assertEqual(content.count('Name="mesh"'), 1)

    def testReadThroughPythonSolver(self):
        """The simulation-loop entry point: any meshio++ format as "input_type"."""
        from KratosMultiphysics.python_solver import PythonSolver

        write_model_part = self.model.CreateModelPart("solver_write")
        _PopulateModelPart(write_model_part)

        with tempfile.TemporaryDirectory() as temp_dir:
            file_name = str(Path(temp_dir) / "solver_input.vtu")
            write_settings = KratosMultiphysics.Parameters("""{"time_series" : "single_file"}""")
            KratosMultiphysics.MeshioPlusPlusIO(file_name, write_settings).WriteModelPart(write_model_part)

            for input_type in ("vtu", "meshio", "auto"):
                read_model_part = self.model.CreateModelPart("solver_read_" + input_type)
                solver = PythonSolver(self.model, KratosMultiphysics.Parameters("""{"echo_level" : 0}"""))
                import_settings = KratosMultiphysics.Parameters("""{
                    "input_type"     : "%s",
                    "input_filename" : "%s"
                }""" % (input_type, file_name.replace("\\", "/")))
                solver._ImportModelPart(read_model_part, import_settings)
                self.assertEqual(read_model_part.NumberOfNodes(), write_model_part.NumberOfNodes())
                self.assertEqual(read_model_part.NumberOfElements(), write_model_part.NumberOfElements())

    def testEntityType(self):
        write_model_part = self.model.CreateModelPart("entity_type_write")
        _PopulateModelPart(write_model_part)

        with tempfile.TemporaryDirectory() as temp_dir:
            file_name = str(Path(temp_dir) / "entity_type.vtu")
            settings = KratosMultiphysics.Parameters("""{
                "time_series" : "single_file",
                "entity_type" : "element"
            }""")
            KratosMultiphysics.MeshioPlusPlusIO(file_name, settings).WriteModelPart(write_model_part)
            read_model_part = self.model.CreateModelPart("entity_type_read")
            KratosMultiphysics.MeshioPlusPlusIO(file_name).ReadModelPart(read_model_part)
            self.assertEqual(read_model_part.NumberOfElements(), 2)   # only the tetrahedra were written
            self.assertEqual(read_model_part.NumberOfConditions(), 0)

    def testWriteIdsAndCellData(self):
        model_part = self.model.CreateModelPart("ids_write")
        _PopulateModelPart(model_part)
        for element in model_part.Elements:
            element.SetValue(KratosMultiphysics.PRESSURE, 2.0 * element.Id)

        with tempfile.TemporaryDirectory() as temp_dir:
            file_name = str(Path(temp_dir) / "ids.vtu")
            settings = KratosMultiphysics.Parameters("""{
                "time_series"                  : "single_file",
                "file_format"                  : "ascii",
                "write_ids"                    : true,
                "nodal_flags"                  : ["TO_ERASE"],
                "element_data_value_variables" : ["PRESSURE"]
            }""")
            KratosMultiphysics.MeshioPlusPlusIO(file_name, settings).WriteModelPart(model_part)
            with open(file_name) as vtu_file:
                content = vtu_file.read()
            for token in ("KRATOS_NODE_ID", "KRATOS_ELEMENT_ID", "PROPERTIES_ID", "TO_ERASE", "PRESSURE"):
                self.assertIn(token, content)
            self.assertIn('format="ascii"', content)  # the "file_format" override

    def testOutputSubModelParts(self):
        model_part = self.model.CreateModelPart("smp_write")
        _PopulateModelPart(model_part)
        volume = model_part.CreateSubModelPart("Volume")
        volume.AddNodes([1, 2, 3, 4, 5])
        volume.AddElements([1, 2])
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 1

        with tempfile.TemporaryDirectory() as temp_dir:
            file_name = str(Path(temp_dir) / "smp.vtu")
            settings = KratosMultiphysics.Parameters("""{"output_sub_model_parts" : true}""")
            KratosMultiphysics.MeshioPlusPlusIO(file_name, settings).WriteModelPart(model_part)
            produced = sorted(f.name for f in Path(temp_dir).iterdir())
            self.assertEqual(produced, ["smp_1.vtu", "smp_smp_write_Volume_1.vtu"])

    def testUnsupportedInputTypeMessage(self):
        from KratosMultiphysics.python_solver import PythonSolver
        model_part = self.model.CreateModelPart("wrong_input_type")
        solver = PythonSolver(self.model, KratosMultiphysics.Parameters("""{"echo_level" : 0}"""))
        import_settings = KratosMultiphysics.Parameters("""{
            "input_type"     : "not_a_format",
            "input_filename" : "irrelevant"
        }""")
        with self.assertRaisesRegex(Exception, 'unsupported "input_type"'):
            solver._ImportModelPart(model_part, import_settings)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()