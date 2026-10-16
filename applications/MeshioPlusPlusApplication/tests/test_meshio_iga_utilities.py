import tempfile
from pathlib import Path

import KratosMultiphysics
import KratosMultiphysics.MeshioPlusPlusApplication as KratosMeshioPlusPlus
import KratosMultiphysics.KratosUnittest as KratosUnittest


def _WriteUnitSquareCadJson(model, file_name):
    """Writes a unit square, a single bilinear quadrilateral, as Kratos CAD JSON."""
    model_part = model.CreateModelPart("UnitSquare")
    props = model_part.CreateNewProperties(1)
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
    model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
    model_part.CreateNewNode(4, 0.0, 1.0, 0.0)
    model_part.CreateNewElement("Element2D4N", 1, [1, 2, 3, 4], props)
    settings = KratosMultiphysics.Parameters('{"time_series" : "single_file"}')
    KratosMeshioPlusPlus.MeshioPlusPlusIO(file_name, settings).WriteModelPart(model_part)


class TestMeshioPlusPlusIgaUtilities(KratosUnittest.TestCase):
    """Tests for the file-based NURBS (isogeometric) utilities."""

    def setUp(self):
        self.model = KratosMultiphysics.Model()
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.square = str(Path(self.temp_dir.name) / "square.cad.json")
        _WriteUnitSquareCadJson(self.model, self.square)

    def test_describe(self):
        info = KratosMeshioPlusPlus.MeshioPlusPlusIgaUtilities.Describe(self.square)
        self.assertEqual(info["kind"].GetString(), "brep")
        self.assertEqual(info["number_of_faces"].GetInt(), 1)
        self.assertEqual(info["number_of_patches"].GetInt(), 1)
        self.assertEqual(info["number_of_trimmed_faces"].GetInt(), 0)
        patch = info["patches"][0]
        self.assertEqual(patch["degrees"].GetVector()[0], 1)
        self.assertEqual(patch["degrees"].GetVector()[1], 1)
        self.assertEqual(patch["number_of_control_points"].GetInt(), 4)
        self.assertFalse(patch["rational"].GetBool())

    def test_sample(self):
        destination = self.model.CreateModelPart("Sampled")
        settings = KratosMultiphysics.Parameters('{"sample_subdivisions" : 2}')
        report = KratosMeshioPlusPlus.MeshioPlusPlusIgaUtilities.Sample(self.square, settings, destination)
        self.assertEqual(report["number_of_nodes"].GetInt(), destination.NumberOfNodes())
        self.assertEqual(destination.NumberOfElements() + destination.NumberOfConditions(), 4)
        self.assertEqual(destination.NumberOfNodes(), 9)
        for node in destination.Nodes:
            self.assertAlmostEqual(node.Z, 0.0, 12)

        with self.assertRaisesRegex(RuntimeError, "must be empty"):
            KratosMeshioPlusPlus.MeshioPlusPlusIgaUtilities.Sample(self.square, settings, destination)

    def test_refine_elevates_and_subdivides_without_moving_the_geometry(self):
        refined_file = str(Path(self.temp_dir.name) / "refined.cad.json")
        settings = KratosMultiphysics.Parameters('{"subdivisions" : [2, 3, 1], "elevate" : [1, 1, 0]}')
        info = KratosMeshioPlusPlus.MeshioPlusPlusIgaUtilities.Refine(self.square, refined_file, settings)
        self.assertTrue(Path(refined_file).exists())
        patch = info["patches"][0]
        self.assertEqual(patch["degrees"].GetVector()[0], 2)
        self.assertEqual(patch["degrees"].GetVector()[1], 2)
        # (2 spans + degree 2) x (3 spans + degree 2) control points
        self.assertEqual(patch["number_of_control_points"].GetInt(), 4 * 5)

        reread = KratosMeshioPlusPlus.MeshioPlusPlusIgaUtilities.Describe(refined_file)
        self.assertEqual(reread["patches"][0]["number_of_control_points"].GetInt(), 20)

        # The geometry is the same unit square
        destination = self.model.CreateModelPart("RefinedSampled")
        KratosMeshioPlusPlus.MeshioPlusPlusIgaUtilities.Sample(
            refined_file, KratosMultiphysics.Parameters("{}"), destination)
        for node in destination.Nodes:
            self.assertAlmostEqual(node.Z, 0.0, 12)
            self.assertTrue(-1e-12 <= node.X <= 1.0 + 1e-12)
            self.assertTrue(-1e-12 <= node.Y <= 1.0 + 1e-12)

        with self.assertRaises(RuntimeError):
            KratosMeshioPlusPlus.MeshioPlusPlusIgaUtilities.Refine(
                self.square, refined_file, KratosMultiphysics.Parameters('{"subdivisions" : [0, 1, 1]}'))

    def test_sampling_a_cubic_model_refuses_what_kratos_has_no_geometry_for(self):
        cubic_file = str(Path(self.temp_dir.name) / "cubic.cad.json")
        KratosMeshioPlusPlus.MeshioPlusPlusIgaUtilities.Refine(
            self.square, cubic_file, KratosMultiphysics.Parameters('{"elevate" : [2, 2, 0]}'))
        # The evaluation and the description need no Kratos geometry
        info = KratosMeshioPlusPlus.MeshioPlusPlusIgaUtilities.Describe(cubic_file)
        self.assertEqual(info["patches"][0]["degrees"].GetVector()[0], 3)
        with self.assertRaisesRegex(RuntimeError, "no registered class"):
            KratosMeshioPlusPlus.MeshioPlusPlusIgaUtilities.Sample(
                cubic_file, KratosMultiphysics.Parameters("{}"), self.model.CreateModelPart("Cubic"))

    def test_evaluate(self):
        for refine in (False, True):
            file_name = self.square
            if refine:
                file_name = str(Path(self.temp_dir.name) / "refined.cad.json")
                KratosMeshioPlusPlus.MeshioPlusPlusIgaUtilities.Refine(
                    self.square, file_name,
                    KratosMultiphysics.Parameters('{"subdivisions" : [2, 2, 1], "elevate" : [1, 1, 0]}'))
            settings = KratosMultiphysics.Parameters('{"parameters" : [[0.0, 0.0], [0.5, 0.5], [1.0, 0.25]], "derivatives" : 1}')
            result = KratosMeshioPlusPlus.MeshioPlusPlusIgaUtilities.Evaluate(file_name, settings)
            points = result["points"].GetMatrix()
            expected = [(0.0, 0.0), (0.5, 0.5), (1.0, 0.25)]
            for i, (x, y) in enumerate(expected):
                self.assertAlmostEqual(points[i, 0], x, 12)
                self.assertAlmostEqual(points[i, 1], y, 12)
                self.assertAlmostEqual(points[i, 2], 0.0, 12)
                normal = result["normals"].GetMatrix()
                self.assertAlmostEqual(abs(normal[i, 2]), 1.0, 12)
            for flag in result["inside"].GetVector():
                self.assertEqual(flag, 1.0)
            # d x / d u = (1, 0, 0), d x / d v = (0, 1, 0), row-major per point
            derivatives = result["derivatives"].GetMatrix()
            self.assertEqual(derivatives.Size2(), 6)
            self.assertAlmostEqual(derivatives[1, 0], 1.0, 12)
            self.assertAlmostEqual(derivatives[1, 4], 1.0, 12)

    def test_evaluate_rejects_a_bad_index_and_a_parameter_out_of_range(self):
        with self.assertRaisesRegex(RuntimeError, "out of range"):
            KratosMeshioPlusPlus.MeshioPlusPlusIgaUtilities.Evaluate(
                self.square, KratosMultiphysics.Parameters('{"index" : 5}'))
        with self.assertRaises(RuntimeError):
            KratosMeshioPlusPlus.MeshioPlusPlusIgaUtilities.Evaluate(
                self.square, KratosMultiphysics.Parameters('{"parameters" : [[7.0, 7.0]]}'))

    def test_a_file_without_a_nurbs_model_is_refused(self):
        with self.assertRaisesRegex(RuntimeError, "carries no NURBS model"):
            KratosMeshioPlusPlus.MeshioPlusPlusIgaUtilities.Describe(str(Path(self.temp_dir.name) / "mesh.vtu"))


if __name__ == "__main__":
    KratosUnittest.main()
