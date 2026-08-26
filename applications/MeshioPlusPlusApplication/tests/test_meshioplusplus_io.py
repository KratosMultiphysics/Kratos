# Import necessary modules from KratosMultiphysics
import KratosMultiphysics
import KratosMultiphysics.MeshioPlusPlusApplication as KratosMeshioPlusPlus
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
        supported = KratosMeshioPlusPlus.MeshioPlusPlusIO.GetSupportedFormats()
        for name in ("vtu", "vtk", "gmsh", "stl", "obj", "xdmf", "ensight", "vtp", "triangle"):
            self.assertIn(name, supported)

        read_formats = KratosMeshioPlusPlus.MeshioPlusPlusIO.GetSupportedReadFormats()
        write_formats = KratosMeshioPlusPlus.MeshioPlusPlusIO.GetSupportedWriteFormats()
        # openfoam round-trips since meshio++ v9.20.0 added the polyMesh writer; it was
        # read-only before, so both directions are pinned.
        self.assertIn("openfoam", read_formats)
        self.assertIn("openfoam", write_formats)
        self.assertIn("gid", read_formats)
        self.assertIn("gid", write_formats)
        for name in ("svg", "tikz"):              # write-only formats
            self.assertIn(name, write_formats)
            self.assertNotIn(name, read_formats)

    def testFormatEnum(self):
        Format = KratosMeshioPlusPlus.MeshioPlusPlusIO.Format
        self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.FormatFromString("vtu"), Format.VTU)
        self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.FormatFromString("auto"), Format.AUTOMATIC)
        self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.FormatFromString("ensight"), Format.ENSIGHT)
        self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.FormatFromString("vtp"), Format.VTP)
        self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.FormatFromString("triangle"), Format.TRIANGLE)
        self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.FormatName(Format.GMSH), "gmsh")
        self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.FormatName(Format.SVG), "svg")
        self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.FormatName(Format.TIKZ), "tikz")
        self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.ResolveFormat("some_mesh.vtu"), Format.VTU)
        self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.ResolveFormat("some_mesh.case"), Format.ENSIGHT)
        self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.ResolveFormat("some_mesh.vtp"), Format.VTP)
        # GiD's compound extension has to beat ".msh"'s own entry, which is gmsh.
        self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.FormatFromString("gid"), Format.GID)
        self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.FormatName(Format.GID), "gid")
        self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.ResolveFormat("results.post.msh"), Format.GID)
        self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.ResolveFormat("results.msh"), Format.GMSH)
        self.assertTrue(KratosMeshioPlusPlus.MeshioPlusPlusIO.IsFormatAvailable(Format.VTU))
        with self.assertRaisesRegex(RuntimeError, "Unknown format"):
            KratosMeshioPlusPlus.MeshioPlusPlusIO.FormatFromString("not_a_format")

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
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name, write_settings).WriteModelPart(write_model_part)
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name).ReadModelPart(read_model_part)

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
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name, write_settings).WriteModelPart(write_model_part)
            self.assertTrue((Path(temp_dir) / "round_trip.geo").exists())
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name).ReadModelPart(read_model_part)

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
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name, settings).WriteModelPart(model_part)
            read_skin = self.model.CreateModelPart("skin_read")
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name).ReadModelPart(read_skin)
            self.assertEqual(read_skin.NumberOfElements(), 6)

            settings = KratosMultiphysics.Parameters("""{"time_series" : "single_file", "skin" : false}""")
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name, settings).WriteModelPart(model_part)
            read_no_skin = self.model.CreateModelPart("no_skin_read")
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name).ReadModelPart(read_no_skin)
            self.assertEqual(read_no_skin.NumberOfElements(), 1)  # only the existing triangle condition

    def testWriteReadRoundTripMed(self):
        Format = KratosMeshioPlusPlus.MeshioPlusPlusIO.Format
        if not KratosMeshioPlusPlus.MeshioPlusPlusIO.IsFormatAvailable(Format.MED):
            self.skipTest("med requires an HDF5-enabled build")
        self._RunWriteReadRoundTrip(".med")

    def _RunWriteReadRoundTrip(self, extension):
        write_model_part = self.model.CreateModelPart("write" + extension.replace(".", "_"))
        read_model_part = self.model.CreateModelPart("read" + extension.replace(".", "_"))
        _PopulateModelPart(write_model_part)

        with tempfile.TemporaryDirectory() as temp_dir:
            file_name = str(Path(temp_dir) / ("round_trip" + extension))
            write_settings = KratosMultiphysics.Parameters("""{"time_series" : "single_file"}""")
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name, write_settings).WriteModelPart(write_model_part)
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name).ReadModelPart(read_model_part)

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
            meshio_io = KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name, settings.Clone())
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
            appending_io = KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name, settings.Clone())
            model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.4
            appending_io.WriteModelPart(model_part)

            with open(file_name) as xdmf_file:
                content = xdmf_file.read()
            self.assertEqual(content.count("<Time Value="), 4)
            self.assertEqual(content.count('Name="mesh"'), 1)

            # Release the series before the temporary directory is removed: a live writer
            # finalizes on destruction and would recreate the file underneath it.
            meshio_io.CloseOutput()
            appending_io.CloseOutput()

    def testMdpaPropertiesRoundTrip(self):
        """Material data survives Kratos -> .mdpa -> Kratos, values included."""
        write_model_part = self.model.CreateModelPart("props_write")
        _PopulateModelPart(write_model_part)
        props = write_model_part.GetProperties()[1]
        props.SetValue(KratosMultiphysics.DENSITY, 7850.0)
        props.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        props.SetValue(KratosMultiphysics.VOLUME_ACCELERATION, [0.0, 0.0, -9.81])

        settings = KratosMultiphysics.Parameters("""{"time_series" : "single_file"}""")

        with tempfile.TemporaryDirectory() as temp_dir:
            file_name = str(Path(temp_dir) / "materials.mdpa")
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name, settings.Clone()).WriteModelPart(write_model_part)

            read_model_part = self.model.CreateModelPart("props_read")
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name).ReadModelPart(read_model_part)

            self.assertTrue(read_model_part.HasProperties(1))
            read_props = read_model_part.GetProperties()[1]
            self.assertAlmostEqual(read_props.GetValue(KratosMultiphysics.DENSITY), 7850.0)
            self.assertEqual(read_props.GetValue(KratosMultiphysics.DOMAIN_SIZE), 3)
            self.assertVectorAlmostEqual(
                read_props.GetValue(KratosMultiphysics.VOLUME_ACCELERATION), [0.0, 0.0, -9.81])

            # In .mdpa a "gmsh:physical" tag stores the properties id, not a physical
            # group, so it must not be turned into a sub model part.
            for sub_model_part in read_model_part.SubModelParts:
                self.assertFalse(sub_model_part.Name.startswith("gmsh_physical_"))

    def testCloseOutput(self):
        """CloseOutput() ends the series, so deleting the output really deletes it."""
        model_part = self.model.CreateModelPart("closing")
        _PopulateModelPart(model_part)

        settings = KratosMultiphysics.Parameters("""{"xdmf_data_format" : "XML"}""")

        with tempfile.TemporaryDirectory() as temp_dir:
            file_name = Path(temp_dir) / "closing.xdmf"
            meshio_io = KratosMeshioPlusPlus.MeshioPlusPlusIO(str(file_name), settings.Clone())
            meshio_io.WriteModelPart(model_part)
            meshio_io.WriteModelPart(model_part)
            meshio_io.CloseOutput()

            self.assertTrue(file_name.is_file())
            self.assertEqual(file_name.read_text().count("<Time Value="), 2)

            # The still-alive IO must not resurrect the deleted output, and a further write
            # starts a fresh series instead of appending to a stale one.
            kratos_utils.DeleteFileIfExisting(str(file_name))
            meshio_io.CloseOutput()  # idempotent
            self.assertFalse(file_name.is_file())

            meshio_io.WriteModelPart(model_part)
            meshio_io.CloseOutput()
            self.assertEqual(file_name.read_text().count("<Time Value="), 1)

    def testReadThroughPythonSolver(self):
        """The simulation-loop entry point: any meshio++ format as "input_type"."""
        from KratosMultiphysics.python_solver import PythonSolver

        write_model_part = self.model.CreateModelPart("solver_write")
        _PopulateModelPart(write_model_part)

        with tempfile.TemporaryDirectory() as temp_dir:
            file_name = str(Path(temp_dir) / "solver_input.vtu")
            write_settings = KratosMultiphysics.Parameters("""{"time_series" : "single_file"}""")
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name, write_settings).WriteModelPart(write_model_part)

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
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name, settings).WriteModelPart(write_model_part)
            read_model_part = self.model.CreateModelPart("entity_type_read")
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name).ReadModelPart(read_model_part)
            self.assertEqual(read_model_part.NumberOfElements(), 2)   # only the tetrahedra were written
            self.assertEqual(read_model_part.NumberOfConditions(), 0)

    def testFileSeriesTimeValues(self):
        """GetNumberOfTimeSteps/GetTimeValues/GetTimeStepIndex must also detect a file
        series (every format but XDMF): mFileName itself is never written in that case."""
        model_part = self.model.CreateModelPart("series_write")
        _PopulateModelPart(model_part)

        with tempfile.TemporaryDirectory() as temp_dir:
            file_name = str(Path(temp_dir) / "series.vtu")
            for step in (1, 2, 5):
                model_part.ProcessInfo[KratosMultiphysics.STEP] = step
                KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name).WriteModelPart(model_part)

            meshio_io = KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name)
            self.assertEqual(meshio_io.GetNumberOfTimeSteps(), 3)
            self.assertVectorAlmostEqual(meshio_io.GetTimeValues(), [1.0, 2.0, 5.0])
            self.assertEqual(meshio_io.GetTimeStepIndex(2.0), 1)
            self.assertEqual(meshio_io.GetTimeStepIndex(3.0), -1)

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
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name, settings).WriteModelPart(model_part)
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
            KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name, settings).WriteModelPart(model_part)
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
        with self.assertRaisesRegex(Exception, 'nsupported "input_type"'):
            solver._ImportModelPart(model_part, import_settings)

    def testReadsGappedMdpaNodeIds(self):
        """Real Kratos decks routinely have gaps; the C++ reader rejected them before v9.13.0."""
        deck = """Begin Properties 0
End Properties

Begin Nodes
    7   0.0   0.0   0.0
   42   1.0   0.0   0.0
   99   0.0   1.0   0.0
End Nodes

Begin Elements Element2D3N
    1   0   7   42   99
End Elements
"""
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "gapped.mdpa"
            path.write_text(deck)

            model_part = self.model.CreateModelPart("gapped")
            KratosMeshioPlusPlus.MeshioPlusPlusIO(str(path)).ReadModelPart(model_part)

        self.assertEqual(model_part.NumberOfNodes(), 3)
        self.assertEqual(model_part.NumberOfElements(), 1)

    def testWriteMdpaIdsPreservesOriginalIds(self):
        """The mdpa writer renumbers to 1..n unless the mesh carries "mdpa:id"."""
        source = self.model.CreateModelPart("gapped_source")
        source.CreateNewNode(7, 0.0, 0.0, 0.0)
        source.CreateNewNode(42, 1.0, 0.0, 0.0)
        source.CreateNewNode(99, 0.0, 1.0, 0.0)
        properties = source.CreateNewProperties(0)
        source.CreateNewElement("Element2D3N", 55, [7, 42, 99], properties)

        with tempfile.TemporaryDirectory() as temp_dir:
            preserved_path = Path(temp_dir) / "preserved.mdpa"
            settings = KratosMultiphysics.Parameters(
                """{"time_series" : "single_file", "write_mdpa_ids" : true}""")
            KratosMeshioPlusPlus.MeshioPlusPlusIO(str(preserved_path), settings).WriteModelPart(source)
            preserved = preserved_path.read_text()

            renumbered_path = Path(temp_dir) / "renumbered.mdpa"
            settings = KratosMultiphysics.Parameters("""{"time_series" : "single_file"}""")
            KratosMeshioPlusPlus.MeshioPlusPlusIO(str(renumbered_path), settings).WriteModelPart(source)
            renumbered = renumbered_path.read_text()

        # The ids cannot collide with the coordinates (0.0 / 1.0), so a substring search is
        # a sound differential oracle here.
        self.assertIn("99", preserved)
        self.assertIn("42", preserved)
        self.assertNotIn("99", renumbered)
        self.assertNotIn("42", renumbered)

    def testSniffFormat(self):
        """Identifies the format from content, not extension."""
        source = self.model.CreateModelPart("to_sniff")
        source.CreateNewNode(1, 0.0, 0.0, 0.0)
        source.CreateNewNode(2, 1.0, 0.0, 0.0)
        source.CreateNewNode(3, 0.0, 1.0, 0.0)
        properties = source.CreateNewProperties(0)
        source.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "mesh.misleading_suffix"
            settings = KratosMultiphysics.Parameters(
                """{"format" : "off", "time_series" : "single_file"}""")
            KratosMeshioPlusPlus.MeshioPlusPlusIO(str(path), settings).WriteModelPart(source)
            self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.SniffFormat(str(path)), "off")

            # Deliberately conservative: "" rather than a guess. mdpa carries no signature
            # the sniffer recognizes.
            opaque = Path(temp_dir) / "opaque.dat"
            opaque.write_text("Begin Nodes\n    1   0.0   0.0   0.0\nEnd Nodes\n")
            self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.SniffFormat(str(opaque)), "")

    def testGidRoundTrip(self):
        Format = KratosMeshioPlusPlus.MeshioPlusPlusIO.Format
        if not KratosMeshioPlusPlus.MeshioPlusPlusIO.IsFormatAvailable(Format.GID):
            self.skipTest("This build has no gidpost (GiD writing needs zlib)")

        write_model_part = self.model.CreateModelPart("Write")
        _PopulateModelPart(write_model_part)

        with tempfile.TemporaryDirectory() as tmp_dir:
            # The ascii flavour writes a .post.msh/.post.res sibling pair.
            file_path = Path(tmp_dir) / "results.post.msh"
            settings = KratosMultiphysics.Parameters(
                '{"time_series" : "single_file", "gid_mode" : "ascii"}')
            io_write = KratosMeshioPlusPlus.MeshioPlusPlusIO(str(file_path), settings)
            io_write.WriteModelPart(write_model_part)
            self.assertTrue(file_path.is_file())

            read_model_part = self.model.CreateModelPart("Read")
            KratosMeshioPlusPlus.MeshioPlusPlusIO(str(file_path)).ReadModelPart(read_model_part)
            self.assertEqual(read_model_part.NumberOfNodes(), write_model_part.NumberOfNodes())
            self.assertGreater(read_model_part.NumberOfElements(), 0)

    def testGidTimeSeries(self):
        Format = KratosMeshioPlusPlus.MeshioPlusPlusIO.Format
        if not KratosMeshioPlusPlus.MeshioPlusPlusIO.IsFormatAvailable(Format.GID):
            self.skipTest("This build has no gidpost (GiD writing needs zlib)")

        model_part = self.model.CreateModelPart("Transient")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        _PopulateModelPart(model_part)

        with tempfile.TemporaryDirectory() as tmp_dir:
            file_path = Path(tmp_dir) / "transient.post.msh"
            settings = KratosMultiphysics.Parameters("""{
                "gid_mode" : "ascii",
                "output_control_type" : "time",
                "nodal_solution_step_data_variables" : ["TEMPERATURE"]
            }""")
            io_write = KratosMeshioPlusPlus.MeshioPlusPlusIO(str(file_path), settings)
            for step in range(3):
                model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.5 * (step + 1)
                for node in model_part.Nodes:
                    node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, step * 10.0 + node.X)
                io_write.WriteModelPart(model_part)
            # Nothing reaches disk until the series is closed: meshio++ pulls every step at once.
            self.assertFalse(file_path.exists())
            io_write.CloseOutput()
            self.assertTrue(file_path.is_file())

            io_read = KratosMeshioPlusPlus.MeshioPlusPlusIO(str(file_path))
            self.assertEqual(io_read.GetNumberOfTimeSteps(), 3)
            self.assertVectorAlmostEqual(io_read.GetTimeValues(), [0.5, 1.0, 1.5])

    def testGidFileSeriesKeepsCompoundExtension(self):
        Format = KratosMeshioPlusPlus.MeshioPlusPlusIO.Format
        if not KratosMeshioPlusPlus.MeshioPlusPlusIO.IsFormatAvailable(Format.GID):
            self.skipTest("This build has no gidpost (GiD writing needs zlib)")

        model_part = self.model.CreateModelPart("Series")
        _PopulateModelPart(model_part)
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 3

        with tempfile.TemporaryDirectory() as tmp_dir:
            file_path = Path(tmp_dir) / "results.post.msh"
            settings = KratosMultiphysics.Parameters(
                '{"time_series" : "file_series", "gid_mode" : "ascii"}')
            KratosMeshioPlusPlus.MeshioPlusPlusIO(str(file_path), settings).WriteModelPart(model_part)

            # The step label goes before the compound extension: a "results.post_3.msh" would no
            # longer resolve to gid at all.
            expected = Path(tmp_dir) / "results_3.post.msh"
            self.assertTrue(expected.is_file())
            self.assertEqual(KratosMeshioPlusPlus.MeshioPlusPlusIO.ResolveFormat(str(expected)),
                             Format.GID)

    def testGidUnknownModeRaises(self):
        with self.assertRaisesRegex(RuntimeError, 'Unknown "gid_mode" setting'):
            KratosMeshioPlusPlus.MeshioPlusPlusIO(
                "unused.post.msh", KratosMultiphysics.Parameters('{"gid_mode" : "not_a_mode"}'))

    def testProvenanceRoundTrip(self):
        model_part = self.model.CreateModelPart("MyStructure")
        _PopulateModelPart(model_part)

        with tempfile.TemporaryDirectory() as tmp_dir:
            file_path = Path(tmp_dir) / "provenance.vtu"
            settings = KratosMultiphysics.Parameters(
                '{"time_series" : "single_file", "provenance" : "best_effort"}')
            KratosMeshioPlusPlus.MeshioPlusPlusIO(str(file_path), settings).WriteModelPart(model_part)

            provenance = KratosMeshioPlusPlus.MeshioPlusPlusIO(str(file_path)).GetProvenance()
            self.assertTrue(provenance["recognised"].GetBool())
            self.assertGreater(provenance["lines"].size(), 0)

    def testUnknownProvenanceModeRaises(self):
        with self.assertRaisesRegex(RuntimeError, 'Unknown "provenance" setting'):
            KratosMeshioPlusPlus.MeshioPlusPlusIO(
                "unused.vtu", KratosMultiphysics.Parameters('{"provenance" : "sometimes"}'))

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()