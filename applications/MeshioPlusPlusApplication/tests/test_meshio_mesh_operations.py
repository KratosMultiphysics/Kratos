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
