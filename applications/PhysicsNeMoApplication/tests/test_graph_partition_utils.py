"""Serial coverage for the graph partitioning helpers.

The distributed behaviour lives in tests/test_mpi_graph_partition.py; this
covers the contract, the serial degenerate case and the guards.
"""

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
from KratosMultiphysics.PhysicsNeMoApplication.distributed import graph_partition_utils
from test_grid_bridge import CreateStructuredTetModelPart


class TestHaloSubgraphSerial(KratosUnittest.TestCase):
    """On one rank every node is owned and the subgraph is the whole graph."""

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Serial", 3, historical_variables=(Kratos.PRESSURE,))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.X + 2.0 * node.Y)

    def test_MatchesBuildGraphWhenNotDistributed(self):
        features, edge_index, edge_features, node_ids, owned_mask = \
            graph_partition_utils.BuildHaloSubgraph(
                self.model_part, field_specs=(("PRESSURE", "node_historical"),))
        reference = graph_bridge.BuildGraph(
            self.model_part, (("PRESSURE", "node_historical"),))

        numpy.testing.assert_allclose(features, reference[0], atol=1e-12)
        numpy.testing.assert_array_equal(edge_index, reference[1])
        numpy.testing.assert_allclose(edge_features, reference[2], atol=1e-12)
        numpy.testing.assert_array_equal(node_ids, reference[3])
        self.assertTrue(owned_mask.all())   # nothing is a halo on one rank

    def test_HaloRingsValidated(self):
        with self.assertRaisesRegex(ValueError, "num_halo_rings"):
            graph_partition_utils.BuildHaloSubgraph(self.model_part, num_halo_rings=0)

    def test_PartitionUnsafeInterfacesAreRefusedOnlyWhenDistributed(self):
        # serial: every interface is fine, nothing to partition
        for interface in ("meshgraphnet", "bistride", "hybrid"):
            graph_partition_utils.CheckPartitionSafeInterface(self.model_part, interface)

    def test_GatherOwnedPredictionsIsIdentityInSerial(self):
        features, _, _, node_ids, owned_mask = graph_partition_utils.BuildHaloSubgraph(
            self.model_part, field_specs=(("PRESSURE", "node_historical"),))
        ids, values = graph_partition_utils.GatherOwnedPredictionsToRank0(
            self.model_part, node_ids, owned_mask, features)
        numpy.testing.assert_array_equal(ids, numpy.sort(node_ids))
        self.assertEqual(values.shape, (len(node_ids), 1))
        # ...and the VALUES must follow the ids. Checking only that the ids
        # come back sorted lets every prediction attach to the wrong node.
        # Compared against the input features directly, not re-gathered, so
        # a permutation inside the gather cannot cancel out.
        numpy.testing.assert_allclose(
            values[:, 0], features[numpy.argsort(node_ids), 0], rtol=1e-12)


if __name__ == '__main__':
    KratosUnittest.main()
