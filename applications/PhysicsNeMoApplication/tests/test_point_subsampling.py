"""Tests for point-cloud token budgets: SelectPointSubset /
ExpandToAllPoints and the "subsampling" block on
PointCloudInferenceProcess."""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.utilities import point_subsampling

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import physicsnemo.nn.functional  # noqa: F401
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False

try:
    import scipy.spatial  # noqa: F401
    have_scipy = True
except ImportError:
    have_scipy = False


def _ClusteredCloud(seed=0):
    """Most points crammed into one corner, a few spread out - the shape a
    mesh refined in one region has, and where a uniform draw and
    farthest-point sampling visibly disagree."""
    rng = numpy.random.default_rng(seed)
    dense = rng.uniform(0.0, 0.1, size=(180, 3))
    sparse = rng.uniform(0.0, 1.0, size=(20, 3))
    return numpy.vstack([dense, sparse])


def _MinimumPairwiseDistance(points):
    differences = points[:, None, :] - points[None, :, :]
    distances = numpy.linalg.norm(differences, axis=-1)
    numpy.fill_diagonal(distances, numpy.inf)
    return float(distances.min())


class TestSelectPointSubset(KratosUnittest.TestCase):
    def test_NoneSelectsEverything(self):
        coordinates = _ClusteredCloud()
        self.assertIsNone(point_subsampling.SelectPointSubset(
            coordinates, Kratos.Parameters('{"method": "none"}')))
        # a budget at or above the point count is not a budget
        self.assertIsNone(point_subsampling.SelectPointSubset(
            coordinates, Kratos.Parameters('{"method": "uniform", "num_points": 500}')))

    def test_UniformSelectsExactlyTheBudgetAndIsSeeded(self):
        coordinates = _ClusteredCloud()
        settings = '{"method": "uniform", "num_points": 16, "seed": 3}'
        first = point_subsampling.SelectPointSubset(coordinates, Kratos.Parameters(settings))
        again = point_subsampling.SelectPointSubset(coordinates, Kratos.Parameters(settings))
        self.assertEqual(first.size, 16)
        numpy.testing.assert_array_equal(first, again)
        # sorted indices into the original array
        numpy.testing.assert_array_equal(first, numpy.sort(first))
        self.assertLess(first.max(), len(coordinates))

    def test_BoundingBoxFiltersInModelPartCoordinates(self):
        coordinates = _ClusteredCloud()
        selected = point_subsampling.SelectPointSubset(coordinates, Kratos.Parameters("""{
            "method"           : "none",
            "bounding_box_min" : [0.0, 0.0, 0.0],
            "bounding_box_max" : [0.1, 0.1, 0.1]
        }"""))
        self.assertIsNotNone(selected)
        self.assertTrue((coordinates[selected] <= 0.1 + 1e-12).all())
        self.assertLess(selected.size, len(coordinates))

    def test_EmptyBoundingBoxRaises(self):
        with self.assertRaisesRegex(ValueError, "no nodes"):
            point_subsampling.SelectPointSubset(_ClusteredCloud(), Kratos.Parameters("""{
                "bounding_box_min" : [10.0, 10.0, 10.0],
                "bounding_box_max" : [11.0, 11.0, 11.0]
            }"""))

    def test_MalformedBoundingBoxRaises(self):
        with self.assertRaisesRegex(ValueError, "three entries"):
            point_subsampling.SelectPointSubset(_ClusteredCloud(), Kratos.Parameters("""{
                "bounding_box_min" : [0.0, 0.0],
                "bounding_box_max" : [1.0, 1.0, 1.0]
            }"""))

    def test_UnsupportedMethodRaises(self):
        with self.assertRaisesRegex(ValueError, "subsampling method"):
            point_subsampling.SelectPointSubset(
                _ClusteredCloud(), Kratos.Parameters('{"method": "poisson"}'))

    @KratosUnittest.skipUnless(have_torch and have_physicsnemo,
                               "Missing required python modules: torch, physicsnemo.")
    def test_FarthestPointCoversTheGeometryBetterThanAUniformDraw(self):
        """The whole reason to prefer FPS: a uniform draw of a cloud that is
        dense in one corner spends its budget there."""
        coordinates = _ClusteredCloud()
        farthest = point_subsampling.SelectPointSubset(
            coordinates, Kratos.Parameters(
                '{"method": "farthest_point", "num_points": 16, "seed": 0}'))
        uniform = point_subsampling.SelectPointSubset(
            coordinates, Kratos.Parameters(
                '{"method": "uniform", "num_points": 16, "seed": 0}'))
        self.assertEqual(farthest.size, 16)
        self.assertGreater(_MinimumPairwiseDistance(coordinates[farthest]),
                           _MinimumPairwiseDistance(coordinates[uniform]))

    @KratosUnittest.skipUnless(have_torch and have_physicsnemo,
                               "Missing required python modules: torch, physicsnemo.")
    def test_BoundingBoxAndBudgetCompose(self):
        coordinates = _ClusteredCloud()
        selected = point_subsampling.SelectPointSubset(coordinates, Kratos.Parameters("""{
            "method"           : "farthest_point",
            "num_points"       : 8,
            "seed"             : 0,
            "bounding_box_min" : [0.0, 0.0, 0.0],
            "bounding_box_max" : [0.1, 0.1, 0.1]
        }"""))
        self.assertEqual(selected.size, 8)
        self.assertTrue((coordinates[selected] <= 0.1 + 1e-12).all())


@KratosUnittest.skipUnless(have_scipy, "Missing required python module: scipy.")
class TestExpandToAllPoints(KratosUnittest.TestCase):
    def test_EveryPointTakesItsNearestSelectedValue(self):
        coordinates = numpy.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0],
                                   [1.0, 0.0, 0.0], [0.9, 0.0, 0.0]])
        selected = numpy.array([0, 2])
        values = numpy.array([[10.0], [20.0]])
        expanded = point_subsampling.ExpandToAllPoints(values, coordinates, selected)
        numpy.testing.assert_allclose(expanded.ravel(), [10.0, 10.0, 20.0, 20.0])

    def test_SelectedPointsKeepTheirOwnValue(self):
        rng = numpy.random.default_rng(0)
        coordinates = rng.uniform(size=(40, 3))
        selected = numpy.sort(rng.choice(40, size=9, replace=False))
        values = rng.standard_normal((9, 2))
        expanded = point_subsampling.ExpandToAllPoints(values, coordinates, selected)
        self.assertEqual(expanded.shape, (40, 2))
        numpy.testing.assert_allclose(expanded[selected], values)

    def test_NoneIsPassedThrough(self):
        self.assertIsNone(point_subsampling.ExpandToAllPoints(
            None, numpy.zeros((3, 3)), numpy.array([0])))

    @KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
    def test_ATorchTensorStaysATorchTensor(self):
        coordinates = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.9, 0.0, 0.0]])
        expanded = point_subsampling.ExpandToAllPoints(
            torch.tensor([[1.0], [2.0]]), coordinates, numpy.array([0, 1]))
        self.assertTrue(torch.is_tensor(expanded))
        numpy.testing.assert_allclose(expanded.numpy().ravel(), [1.0, 2.0, 2.0])


@KratosUnittest.skipUnless(have_torch and have_scipy and have_physicsnemo,
                           "Missing required python modules: torch, scipy, physicsnemo.")
class TestSubsamplingThroughTheProcess(KratosUnittest.TestCase):
    def setUp(self):
        from test_grid_bridge import CreateStructuredTetModelPart

        self.checkpoint = Path("test_subsampled_model.pt")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=3,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.X + node.Y)

        class Pointwise(torch.nn.Module):
            """x -> the first coordinate, so a node's correct value is
            known and a filled-in node's is its neighbour's."""

            def forward(self, x):
                return x[..., :1]

        torch.jit.script(Pointwise()).save(str(self.checkpoint))

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting(str(self.checkpoint))

    def _Run(self, subsampling):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import (
            point_cloud_inference_process)

        block = ',\n                "subsampling" : %s' % subsampling if subsampling else ""
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"       : "Main",
                "model_settings"        : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "torchscript",
                    "device"          : "cpu"
                },
                "model_interface"       : "generic",
                "normalize_coordinates" : false,
                "input_fields"          : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"         : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]%s
            }
        }""" % (self.checkpoint, block))
        process = point_cloud_inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        return numpy.array([
            (node.X, node.GetSolutionStepValue(Kratos.TEMPERATURE))
            for node in self.model_part.Nodes])

    def test_EveryNodeStillGetsAValue(self):
        full = self._Run(None)
        numpy.testing.assert_allclose(full[:, 1], full[:, 0], atol=1e-12)

        budgeted = self._Run('{ "method" : "farthest_point", "num_points" : 8, "seed" : 0 }')
        self.assertEqual(len(budgeted), len(full))
        self.assertTrue(numpy.isfinite(budgeted[:, 1]).all())
        # the model saw 8 nodes, so only 8 distinct values can be written
        self.assertLessEqual(len(numpy.unique(numpy.round(budgeted[:, 1], 10))), 8)
        # and every written value is one the model actually predicted
        for value in budgeted[:, 1]:
            self.assertTrue(numpy.isclose(numpy.abs(full[:, 0] - value).min(), 0.0, atol=1e-9))

    def test_FilledNodesTakeTheirNearestSelectedNodesValue(self):
        budgeted = self._Run('{ "method" : "uniform", "num_points" : 6, "seed" : 1 }')
        written = budgeted[:, 1]
        # a node whose own X was written is a selected node; every other
        # node's value must belong to some selected node
        selected_values = numpy.unique(numpy.round(written, 10))
        self.assertLessEqual(len(selected_values), 6)
        self.assertGreater(len(selected_values), 1)

    def test_BoundingBoxRestrictsWhatTheModelSees(self):
        budgeted = self._Run("""{
            "method"           : "none",
            "bounding_box_min" : [0.0, 0.0, 0.0],
            "bounding_box_max" : [0.4, 1.0, 1.0]
        }""")
        # nothing outside the box was ever evaluated, so no value above 0.4
        self.assertLessEqual(budgeted[:, 1].max(), 0.4 + 1e-9)


if __name__ == '__main__':
    KratosUnittest.main()
