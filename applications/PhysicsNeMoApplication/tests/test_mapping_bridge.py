"""Tests for mapping_bridge: MappingApplication-based transfer onto ML grids
(availability-gated on the compiled MappingApplication)."""

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.bridges import mapping_bridge
from test_grid_bridge import CreateStructuredTetModelPart

have_mapping = kratos_utils.CheckIfApplicationsAvailable("MappingApplication")


def _SetLinearField(model_part, variable):
    for node in model_part.Nodes:
        node.SetSolutionStepValue(variable, 1.0 + 2.0 * node.X + 3.0 * node.Y - node.Z)


class TestBackgroundGrid(KratosUnittest.TestCase):
    """CreateBackgroundGridModelPart and GatherGridArray are pure Kratos."""

    def setUp(self):
        self.model = Kratos.Model()

    def test_GridPartShapes(self):
        grid_part = mapping_bridge.CreateBackgroundGridModelPart(
            self.model, "Grid", (numpy.zeros(3), numpy.ones(3)), divisions=3,
            historical_variables=(Kratos.TEMPERATURE,))
        self.assertEqual(grid_part.NumberOfNodes(), 4 ** 3)
        self.assertGreater(grid_part.NumberOfElements(), 0)

    def test_GatherGridArrayLatticeOrder(self):
        grid_part = mapping_bridge.CreateBackgroundGridModelPart(
            self.model, "Grid", (numpy.zeros(3), numpy.ones(3)), divisions=2,
            historical_variables=(Kratos.TEMPERATURE,))
        for node in grid_part.Nodes:
            node.SetSolutionStepValue(Kratos.TEMPERATURE, node.X + 10.0 * node.Y + 100.0 * node.Z)

        grid = mapping_bridge.GatherGridArray(grid_part, ["TEMPERATURE"], (3, 3, 3))
        self.assertEqual(grid.shape, (1, 3, 3, 3))
        # axes are (x, y, z): stepping each index changes only its coordinate
        self.assertAlmostEqual(grid[0, 1, 0, 0] - grid[0, 0, 0, 0], 0.5, places=12)
        self.assertAlmostEqual(grid[0, 0, 1, 0] - grid[0, 0, 0, 0], 5.0, places=12)
        self.assertAlmostEqual(grid[0, 0, 0, 1] - grid[0, 0, 0, 0], 50.0, places=12)

    def test_InvalidArgumentsRaise(self):
        with self.assertRaisesRegex(ValueError, "divisions"):
            mapping_bridge.CreateBackgroundGridModelPart(
                self.model, "G1", (numpy.zeros(3), numpy.ones(3)), divisions=0)
        with self.assertRaisesRegex(ValueError, "dimension"):
            mapping_bridge.CreateBackgroundGridModelPart(
                self.model, "G2", (numpy.zeros(3), numpy.ones(3)), divisions=2, dimension=4)


@KratosUnittest.skipUnless(have_mapping, "Missing required application: MappingApplication.")
class TestMappingBridge(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.origin_part = CreateStructuredTetModelPart(
            self.model, "Origin", divisions=4,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        _SetLinearField(self.origin_part, Kratos.PRESSURE)
        self.grid_part = mapping_bridge.CreateBackgroundGridModelPart(
            self.model, "Grid", (numpy.zeros(3), numpy.ones(3)), divisions=3,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))

    def test_LinearFieldMapsExactly(self):
        # nearest_element interpolates through the origin's shape functions:
        # exact for a linear field on the tet mesh.
        bridge = mapping_bridge.MappingBridge(self.origin_part, self.grid_part)
        bridge.MapFields([("PRESSURE", "PRESSURE")])

        grid = mapping_bridge.GatherGridArray(self.grid_part, ["PRESSURE"], (4, 4, 4))
        axis = numpy.linspace(0.0, 1.0, 4)
        x, y, z = numpy.meshgrid(axis, axis, axis, indexing="ij")
        numpy.testing.assert_allclose(grid[0], 1.0 + 2.0 * x + 3.0 * y - z, atol=1e-10)

    def test_InverseMapRoundTrip(self):
        bridge = mapping_bridge.MappingBridge(self.origin_part, self.grid_part)
        bridge.MapFields([("PRESSURE", "PRESSURE")])
        # bring the grid values back into a different origin variable
        bridge.InverseMapFields([("TEMPERATURE", "PRESSURE")])

        values = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.origin_part.Nodes])
        reference = numpy.array([
            node.GetSolutionStepValue(Kratos.PRESSURE) for node in self.origin_part.Nodes])
        # both meshes resolve the linear field: the round trip is near-exact
        self.assertLess(numpy.abs(values - reference).max(), 1e-6)

    def test_NearestNeighborMapperOption(self):
        bridge = mapping_bridge.MappingBridge(
            self.origin_part, self.grid_part,
            Kratos.Parameters("""{ "mapper_type" : "nearest_neighbor" }"""))
        bridge.MapFields([("PRESSURE", "PRESSURE")])
        values = numpy.array([
            node.GetSolutionStepValue(Kratos.PRESSURE) for node in self.grid_part.Nodes])
        # Exact, and it DISCRIMINATES: the origin lattice is k/4 and the grid
        # k/3, so each grid node takes its component-wise nearest origin node
        # (no ties). Interpolating instead - i.e. silently running
        # nearest_element - differs by up to 0.5 here.
        coordinates = numpy.array([[node.X, node.Y, node.Z]
                                   for node in self.grid_part.Nodes])
        nearest = numpy.round(4.0 * coordinates) / 4.0
        expected = 1.0 + 2.0 * nearest[:, 0] + 3.0 * nearest[:, 1] - nearest[:, 2]
        numpy.testing.assert_allclose(values, expected, atol=1e-12)
        interpolated = (1.0 + 2.0 * coordinates[:, 0] + 3.0 * coordinates[:, 1]
                        - coordinates[:, 2])
        self.assertGreater(numpy.abs(values - interpolated).max(), 0.1)
        self.assertGreater(numpy.abs(values).max(), 0.0)


if __name__ == '__main__':
    KratosUnittest.main()
