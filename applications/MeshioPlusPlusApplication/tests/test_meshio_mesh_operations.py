import tempfile
from pathlib import Path

import KratosMultiphysics
import KratosMultiphysics.MeshioPlusPlusApplication as KratosMeshioPlusPlus
import KratosMultiphysics.KratosUnittest as KratosUnittest


def _CreateCubeOfTetrahedra(model_part):
    """A unit cube split into six linear tetrahedra."""
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
    model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
    model_part.CreateNewNode(4, 0.0, 1.0, 0.0)
    model_part.CreateNewNode(5, 0.0, 0.0, 1.0)
    model_part.CreateNewNode(6, 1.0, 0.0, 1.0)
    model_part.CreateNewNode(7, 1.0, 1.0, 1.0)
    model_part.CreateNewNode(8, 0.0, 1.0, 1.0)

    properties = model_part.CreateNewProperties(1)
    connectivities = [
        [1, 2, 3, 7], [1, 2, 7, 6], [1, 6, 7, 5],
        [1, 3, 4, 7], [1, 4, 8, 7], [1, 8, 5, 7],
    ]
    for i, nodes in enumerate(connectivities):
        model_part.CreateNewElement("Element3D4N", i + 1, nodes, properties)


def _CreateTriangulatedSquare(model_part):
    """A single quad (two triangles sharing edge 1-3) on the z=0 plane."""
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
    model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
    model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

    properties = model_part.CreateNewProperties(1)
    model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)
    model_part.CreateNewElement("Element2D3N", 2, [1, 3, 4], properties)


class TestMeshioPlusPlusMeshOperations(KratosUnittest.TestCase):
    def setUp(self):
        self.model = KratosMultiphysics.Model()
        self.source = self.model.CreateModelPart("Source")
        _CreateCubeOfTetrahedra(self.source)

    def _Execute(self, operation, settings=None):
        destination = self.model.CreateModelPart("Destination_" + operation)
        parameters = KratosMultiphysics.Parameters("{}")
        parameters.AddString("operation", operation)
        if settings is not None:
            for key in settings.keys():
                parameters.AddValue(key, settings[key])
        report = KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.Execute(
            self.source, parameters, destination)
        return destination, report

    def test_supported_operations(self):
        operations = KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.GetSupportedOperations()
        for expected in ("clean", "transform", "refine", "partition", "extract_skin", "stats"):
            self.assertIn(expected, operations)

    def test_unknown_operation_raises(self):
        with self.assertRaisesRegex(RuntimeError, "Unknown meshio\\+\\+ operation"):
            self._Execute("this_is_not_an_operation")

    def test_clean_preserves_a_clean_mesh(self):
        destination, _ = self._Execute("clean")
        self.assertEqual(destination.NumberOfNodes(), self.source.NumberOfNodes())
        self.assertEqual(destination.NumberOfElements(), self.source.NumberOfElements())

    def test_transform_translates_the_nodes(self):
        settings = KratosMultiphysics.Parameters("""{"translation" : [1.0, 2.0, 3.0]}""")
        destination, _ = self._Execute("transform", settings)
        self.assertEqual(destination.NumberOfNodes(), self.source.NumberOfNodes())
        # The node created first sits at the origin, so it lands on the translation itself
        node = destination.GetNode(1)
        self.assertAlmostEqual(node.X, 1.0)
        self.assertAlmostEqual(node.Y, 2.0)
        self.assertAlmostEqual(node.Z, 3.0)

    def test_extract_skin_produces_the_boundary(self):
        destination, _ = self._Execute("extract_skin")
        # The cube's surface is 12 triangles, and every node is on the boundary
        self.assertEqual(destination.NumberOfNodes(), 8)
        self.assertGreater(destination.NumberOfElements() + destination.NumberOfConditions(), 0)

    def test_refine_increases_the_cell_count(self):
        destination, _ = self._Execute("refine")
        self.assertGreater(destination.NumberOfElements(), self.source.NumberOfElements())

    # Selective refinement (meshio++ >= 9.5.0 - see MESHIOPLUSPLUS_VERSION_AT_LEAST on the C++
    # side). Two triangles sharing edge 1-3: selecting only cell 0 fully splits it into 4
    # (RedGreen, the default closure); cell 1 sees only its shared edge bisected and splits
    # into 2, so the total is 6 - not the 8 a uniform refine of both would give.
    def test_refine_select_by_explicit_cells(self):
        square = self.model.CreateModelPart("Square")
        _CreateTriangulatedSquare(square)
        destination = self.model.CreateModelPart("RefineByCells")

        settings = KratosMultiphysics.Parameters("""{"operation" : "refine", "cells" : [0]}""")
        KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.Execute(square, settings, destination)

        self.assertEqual(destination.NumberOfElements(), 6)

    def test_refine_select_by_region(self):
        square = self.model.CreateModelPart("SquareForRegion")
        _CreateTriangulatedSquare(square)
        region = square.CreateSubModelPart("target_region")
        region.AddNodes([1, 2, 3])
        region.AddElements([1])
        destination = self.model.CreateModelPart("RefineByRegion")

        settings = KratosMultiphysics.Parameters("""{"operation" : "refine", "region" : "target_region"}""")
        KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.Execute(square, settings, destination)

        self.assertEqual(destination.NumberOfElements(), 6)

    def test_refine_closure_redgreen_vs_propagate(self):
        square = self.model.CreateModelPart("SquareForClosure")
        _CreateTriangulatedSquare(square)

        redgreen = self.model.CreateModelPart("RefineRedGreen")
        redgreen_settings = KratosMultiphysics.Parameters(
            """{"operation" : "refine", "cells" : [0], "closure" : "redgreen"}""")
        KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.Execute(square, redgreen_settings, redgreen)

        propagate = self.model.CreateModelPart("RefinePropagate")
        propagate_settings = KratosMultiphysics.Parameters(
            """{"operation" : "refine", "cells" : [0], "closure" : "propagate"}""")
        KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.Execute(square, propagate_settings, propagate)

        self.assertEqual(redgreen.NumberOfElements(), 6)
        self.assertEqual(propagate.NumberOfElements(), 8)

    def test_refine_two_selectors_raises(self):
        square = self.model.CreateModelPart("SquareForError")
        _CreateTriangulatedSquare(square)
        destination = self.model.CreateModelPart("RefineError")

        settings = KratosMultiphysics.Parameters(
            """{"operation" : "refine", "cells" : [0], "region" : "whatever"}""")
        with self.assertRaises(RuntimeError):
            KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.Execute(square, settings, destination)

    def test_partition_creates_one_model_part_per_piece(self):
        settings = KratosMultiphysics.Parameters("""{"number_of_parts" : 2}""")
        destination, report = self._Execute("partition", settings)
        self.assertEqual(report["number_of_pieces"].GetInt(), 2)
        total_elements = 0
        for i in range(report["pieces"].size()):
            name = report["pieces"][i]["name"].GetString()
            self.assertTrue(self.model.HasModelPart(name))
            total_elements += self.model[name].NumberOfElements()
        # A partition is a partition of unity when there are no ghost layers
        self.assertEqual(total_elements, self.source.NumberOfElements())

    def test_split_by_type(self):
        settings = KratosMultiphysics.Parameters("""{"split_by" : "type"}""")
        destination, report = self._Execute("split", settings)
        # A single cell type, so a single piece
        self.assertEqual(report["number_of_pieces"].GetInt(), 1)
        name = report["pieces"][0]["name"].GetString()
        self.assertTrue(self.model.HasModelPart(name))
        self.assertEqual(self.model[name].NumberOfElements(), self.source.NumberOfElements())

    def test_statistics(self):
        report = KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.ComputeStatistics(self.source)
        self.assertEqual(report["number_of_points"].GetInt(), 8)
        self.assertEqual(report["number_of_cells"].GetInt(), 6)
        self.assertAlmostEqual(report["unsigned_volume"].GetDouble(), 1.0, places=10)

    def test_quality(self):
        report = KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.ComputeQuality(self.source)
        self.assertEqual(report["number_of_cells"].GetInt(), 6)
        self.assertTrue(report.Has("metrics"))

    def test_bandwidth(self):
        self.assertGreater(
            KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.ComputeBandwidth(self.source), 0)

    def test_reorder_reports_the_bandwidth(self):
        settings = KratosMultiphysics.Parameters("""{"method" : "rcm"}""")
        destination, report = self._Execute("reorder", settings)
        self.assertEqual(destination.NumberOfNodes(), self.source.NumberOfNodes())
        self.assertTrue(report.Has("bandwidth_before"))
        self.assertTrue(report.Has("bandwidth_after"))

    def test_diff_of_a_model_part_with_itself(self):
        report = KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.Diff(self.source, self.source)
        self.assertTrue(report["equal"].GetBool())

    def test_merge(self):
        other = self.model.CreateModelPart("Other")
        _CreateCubeOfTetrahedra(other)
        destination = self.model.CreateModelPart("Merged")
        report = KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.Merge(
            [self.source, other], KratosMultiphysics.Parameters("{}"), destination)
        self.assertEqual(report["number_of_sources"].GetInt(), 2)
        self.assertEqual(destination.NumberOfElements(), 12)

    def test_isosurface_needs_a_staged_array(self):
        # Field data is opt-in: with no "*_variables" setting listing it, TEMPERATURE is
        # never staged as point_data even though it is a registered Variable, so the array
        # can never resolve.
        settings = KratosMultiphysics.Parameters("""{"array_name" : "TEMPERATURE"}""")
        with self.assertRaises(RuntimeError):
            self._Execute("isosurface", settings)

    def test_isosurface_with_staged_field_data(self):
        for node in self.source.Nodes:
            node.SetValue(KratosMultiphysics.TEMPERATURE, node.Z)

        settings = KratosMultiphysics.Parameters("""{
            "array_name" : "TEMPERATURE",
            "isovalues" : [0.5],
            "nodal_data_value_variables" : ["TEMPERATURE"]
        }""")
        destination, _ = self._Execute("isosurface", settings)

        self.assertGreater(destination.NumberOfNodes(), 0)
        # "The contoured field is exactly the isovalue on the cut points" - so the write-back
        # path is exercised too.
        for node in destination.Nodes:
            self.assertAlmostEqual(node.GetValue(KratosMultiphysics.TEMPERATURE), 0.5, places=10)

    def test_attach_quality_metrics_do_not_survive_the_write_back(self):
        # attach_quality's own array names ("quality:scaled_jacobian", ...) are not registered
        # Variables, so the metrics are computed but unretrievable from the destination - only
        # the topology survives.
        destination, _ = self._Execute("attach_quality")
        self.assertEqual(destination.NumberOfElements(), self.source.NumberOfElements())

    def test_field_data_round_trips_through_clean(self):
        for node in self.source.Nodes:
            node.SetValue(KratosMultiphysics.TEMPERATURE, 10.0 * node.Id)

        settings = KratosMultiphysics.Parameters("""{"nodal_data_value_variables" : ["TEMPERATURE"]}""")
        destination, _ = self._Execute("clean", settings)

        self.assertEqual(destination.NumberOfNodes(), self.source.NumberOfNodes())
        for node in destination.Nodes:
            self.assertAlmostEqual(node.GetValue(KratosMultiphysics.TEMPERATURE), 10.0 * node.Id, places=10)

    def test_data_calc(self):
        for node in self.source.Nodes:
            node.SetValue(KratosMultiphysics.TEMPERATURE, node.Z)

        settings = KratosMultiphysics.Parameters("""{
            "expression" : "TEMPERATURE * 2.0",
            "output" : "TEMPERATURE",
            "output_overwrite" : true,
            "nodal_data_value_variables" : ["TEMPERATURE"]
        }""")
        destination, _ = self._Execute("data_calc", settings)

        for node in destination.Nodes:
            self.assertAlmostEqual(node.GetValue(KratosMultiphysics.TEMPERATURE), 2.0 * node.Z, places=10)

    def test_data_condition_clamps_values(self):
        for node in self.source.Nodes:
            node.SetValue(KratosMultiphysics.TEMPERATURE, node.Z)

        settings = KratosMultiphysics.Parameters("""{
            "names" : ["TEMPERATURE"],
            "mode" : "clamp",
            "lo" : 0.0,
            "hi" : 0.5,
            "nodal_data_value_variables" : ["TEMPERATURE"]
        }""")
        destination, _ = self._Execute("data_condition", settings)

        for node in destination.Nodes:
            self.assertLessEqual(node.GetValue(KratosMultiphysics.TEMPERATURE), 0.5 + 1e-10)

    def test_data_manage_rename(self):
        # Renaming onto a name that is itself a registered Variable is what lets the result
        # come back through the write-back constraint.
        for node in self.source.Nodes:
            node.SetValue(KratosMultiphysics.TEMPERATURE, 42.0)

        settings = KratosMultiphysics.Parameters("""{
            "rename" : [{"location" : "point", "from" : "TEMPERATURE", "to" : "PRESSURE"}],
            "nodal_data_value_variables" : ["TEMPERATURE"]
        }""")
        destination, report = self._Execute("data_manage", settings)

        self.assertEqual(report["renamed"].size(), 1)
        for node in destination.Nodes:
            self.assertAlmostEqual(node.GetValue(KratosMultiphysics.PRESSURE), 42.0, places=10)

    def test_data_info_reports_array_statistics(self):
        for node in self.source.Nodes:
            node.SetValue(KratosMultiphysics.TEMPERATURE, node.Z)

        settings = KratosMultiphysics.Parameters("""{"nodal_data_value_variables" : ["TEMPERATURE"]}""")
        destination, report = self._Execute("data_info", settings)

        self.assertEqual(report["number_of_point_data"].GetInt(), 1)
        self.assertEqual(report["arrays"][0]["name"].GetString(), "TEMPERATURE")
        self.assertEqual(report["arrays"][0]["num_values"].GetInt(), self.source.NumberOfNodes())
        # Report-only: the destination is never touched
        self.assertEqual(destination.NumberOfNodes(), 0)

    def test_point_data_to_cell_data(self):
        for node in self.source.Nodes:
            node.SetValue(KratosMultiphysics.TEMPERATURE, 5.0)

        settings = KratosMultiphysics.Parameters("""{
            "names" : ["TEMPERATURE"],
            "nodal_data_value_variables" : ["TEMPERATURE"]
        }""")
        destination, _ = self._Execute("point_data_to_cell_data", settings)

        self.assertEqual(destination.NumberOfElements(), self.source.NumberOfElements())
        for element in destination.Elements:
            self.assertAlmostEqual(element.GetValue(KratosMultiphysics.TEMPERATURE), 5.0, places=10)

    def test_cell_data_to_point_data(self):
        for element in self.source.Elements:
            element.SetValue(KratosMultiphysics.DENSITY, 7850.0)

        settings = KratosMultiphysics.Parameters("""{
            "names" : ["DENSITY"],
            "element_data_value_variables" : ["DENSITY"]
        }""")
        destination, _ = self._Execute("cell_data_to_point_data", settings)

        self.assertEqual(destination.NumberOfNodes(), self.source.NumberOfNodes())
        for node in destination.Nodes:
            self.assertAlmostEqual(node.GetValue(KratosMultiphysics.DENSITY), 7850.0, places=10)

    def test_interpolate_nearest(self):
        # Source and target are geometrically identical, so nearest sampling is exact. The
        # target's own (unset, defaulted-to-zero) TEMPERATURE forces the name conflict that
        # "on_conflict" : "overwrite" then resolves.
        for node in self.source.Nodes:
            node.SetValue(KratosMultiphysics.TEMPERATURE, node.Z)
        target = self.model.CreateModelPart("Target")
        _CreateCubeOfTetrahedra(target)
        destination = self.model.CreateModelPart("InterpolatedDestination")

        settings = KratosMultiphysics.Parameters("""{
            "method" : "nearest",
            "on_conflict" : "overwrite",
            "nodal_data_value_variables" : ["TEMPERATURE"]
        }""")
        report = KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.Interpolate(
            self.source, target, settings, destination)

        self.assertEqual(report["number_of_nodes"].GetInt(), target.NumberOfNodes())
        for node in destination.Nodes:
            self.assertAlmostEqual(node.GetValue(KratosMultiphysics.TEMPERATURE), node.Z, places=10)

    def test_interpolate_conflict_raises_by_default(self):
        for node in self.source.Nodes:
            node.SetValue(KratosMultiphysics.TEMPERATURE, node.Z)
        target = self.model.CreateModelPart("TargetConflict")
        _CreateCubeOfTetrahedra(target)
        destination = self.model.CreateModelPart("ConflictDestination")

        settings = KratosMultiphysics.Parameters("""{
            "method" : "nearest",
            "nodal_data_value_variables" : ["TEMPERATURE"]
        }""")
        with self.assertRaises(RuntimeError):
            KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.Interpolate(
                self.source, target, settings, destination)

    def test_entity_names_survive_a_write(self):
        """The registered name must not degrade to the generic cell-type name.

        Written through mdpa, which is the format that carries the Kratos entity name
        verbatim, so the assertion is on what actually reached the file.
        """
        source = self.model.CreateModelPart("Named")
        source.CreateNewNode(1, 0.0, 0.0, 0.0)
        source.CreateNewNode(2, 1.0, 0.0, 0.0)
        source.CreateNewNode(3, 0.0, 1.0, 0.0)
        source.CreateNewNode(4, 0.0, 0.0, 1.0)
        properties = source.CreateNewProperties(1)
        # A concrete element name, not the generic "Element3D4N" the cell type would give
        source.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], properties)

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "named.mdpa"
            # "single_file" so the path is used verbatim, instead of the default file
            # series that appends a step label to the stem
            settings = KratosMultiphysics.Parameters("""{"time_series" : "single_file"}""")
            KratosMeshioPlusPlus.MeshioPlusPlusIO(str(path), settings).WriteModelPart(source)
            written = path.read_text()

        self.assertIn("Begin Elements", written)
        self.assertIn("Element3D4N", written)


if __name__ == "__main__":
    KratosUnittest.main()
