import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
class TestGraphBridge(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        props = self.model_part.CreateNewProperties(1)
        # Two tetrahedra sharing the face (2, 3, 4).
        for i, xyz in enumerate([
                (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1.0, 1.0, 1.0)]):
            node = self.model_part.CreateNewNode(i + 1, *xyz)
            node.SetSolutionStepValue(Kratos.PRESSURE, float(node.Id))
            node.SetSolutionStepValue(Kratos.VELOCITY, [node.X, node.Y, node.Z])
        self.model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], props)
        self.model_part.CreateNewElement("Element3D4N", 2, [2, 5, 3, 4], props)

    def test_EdgeSet(self):
        _, edge_index, _, node_ids = graph_bridge.BuildGraph(self.model_part)
        # Two tets sharing a face: 6 + 6 - 3 shared = 9 unique edges, 18 directed.
        self.assertEqual(edge_index.shape, (2, 18))
        id_of = dict(enumerate(node_ids))
        directed = {(id_of[i], id_of[j]) for i, j in edge_index.T}
        undirected = {frozenset(edge) for edge in directed}
        self.assertEqual(len(undirected), 9)
        expected = {frozenset(e) for e in
                    [(1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4),  # tet 1
                     (2, 5), (5, 3), (5, 4)]}                          # tet 2 additions
        self.assertEqual(undirected, expected)
        # bidirectional: each directed edge's reverse is present
        self.assertTrue(all((b, a) in directed for a, b in directed))

    def test_EdgeFeatures(self):
        _, edge_index, edge_features, node_ids = graph_bridge.BuildGraph(self.model_part)
        coordinates = {node.Id: numpy.array([node.X, node.Y, node.Z]) for node in self.model_part.Nodes}
        for column in range(edge_index.shape[1]):
            i, j = edge_index[:, column]
            relative = coordinates[node_ids[j]] - coordinates[node_ids[i]]
            self.assertTrue(numpy.allclose(edge_features[column, :3], relative))
            self.assertAlmostEqual(edge_features[column, 3], numpy.linalg.norm(relative), places=12)

    def test_NodeFeatureLayout(self):
        node_features, _, _, node_ids = graph_bridge.BuildGraph(
            self.model_part,
            [("PRESSURE", "node_historical"), ("VELOCITY", "node_historical")])
        self.assertEqual(node_features.shape, (5, 4))  # 1 + 3 channels
        for row, node_id in enumerate(node_ids):
            node = self.model_part.GetNode(int(node_id))
            self.assertAlmostEqual(node_features[row, 0], float(node.Id))
            self.assertTrue(numpy.allclose(node_features[row, 1:], [node.X, node.Y, node.Z]))

    def test_TheSplitHelpersAgreeWithTheMeshItself(self):
        """GatherNodeFeatures/ComputeEdgeFeatures are BuildGraph's per-step half.

        Checked against the mesh, NOT against BuildGraph: BuildGraph now
        delegates to both, so comparing them would move together under any
        bug and prove nothing. The fixture stores PRESSURE = node Id, which
        pins row order independently; edge features are recomputed here from
        the node coordinates.
        """
        _, edge_index, _, node_ids = graph_bridge.BuildGraph(
            self.model_part, (("PRESSURE", "node_historical"),))

        features = graph_bridge.GatherNodeFeatures(
            self.model_part, (("PRESSURE", "node_historical"),), len(node_ids))
        # row k must carry node_ids[k]'s own value, not some other node's
        numpy.testing.assert_allclose(
            features[:, 0], node_ids.astype(float), rtol=0.0, atol=0.0)

        coordinates = numpy.array(
            [[node.X, node.Y, node.Z] for node in self.model_part.Nodes])
        expected_relative = coordinates[edge_index[1]] - coordinates[edge_index[0]]
        edge_features = graph_bridge.ComputeEdgeFeatures(self.model_part, edge_index)
        numpy.testing.assert_allclose(
            edge_features[:, :3], expected_relative, rtol=0.0, atol=0.0)
        numpy.testing.assert_allclose(
            edge_features[:, 3],
            numpy.linalg.norm(expected_relative, axis=1), rtol=0.0, atol=0.0)

    def test_GatherNodeFeaturesHandlesAnEmptySpec(self):
        node_ids = graph_bridge.BuildGraph(self.model_part)[3]
        empty = graph_bridge.GatherNodeFeatures(self.model_part, (), len(node_ids))
        self.assertEqual(empty.shape, (len(node_ids), 0))
        # and without being told the row count
        self.assertEqual(
            graph_bridge.GatherNodeFeatures(self.model_part, ()).shape,
            (self.model_part.NumberOfNodes(), 0))

    def test_PrecomputedScatterRowsMatchTheInternalOnes(self):
        _, _, _, node_ids = graph_bridge.BuildGraph(self.model_part)
        values = numpy.arange(len(node_ids), dtype=float).reshape(-1, 1)

        spec = (("PRESSURE", "node_historical"),)
        graph_bridge.ScatterNodeFeatures(self.model_part, node_ids, values, spec)
        internal = [n.GetSolutionStepValue(Kratos.PRESSURE) for n in self.model_part.Nodes]

        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 0.0)
        rows = graph_bridge.BuildScatterRows(self.model_part, node_ids)
        graph_bridge.ScatterNodeFeatures(
            self.model_part, node_ids, values, spec, rows=rows)
        passed_in = [n.GetSolutionStepValue(Kratos.PRESSURE) for n in self.model_part.Nodes]

        numpy.testing.assert_allclose(passed_in, internal, rtol=0.0, atol=0.0)

    def test_EmptyFieldSpec(self):
        node_features, _, _, _ = graph_bridge.BuildGraph(self.model_part)
        self.assertEqual(node_features.shape, (5, 0))

    def test_ScatterRoundTrip(self):
        node_features, _, _, node_ids = graph_bridge.BuildGraph(
            self.model_part, [("VELOCITY", "node_historical")])
        graph_bridge.ScatterNodeFeatures(
            self.model_part, node_ids, 2.0 * node_features, [("VELOCITY", "node_historical")])
        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(
                node.GetSolutionStepValue(Kratos.VELOCITY), [2.0 * node.X, 2.0 * node.Y, 2.0 * node.Z])

    def test_HexahedronEdges(self):
        model_part = self.model.CreateModelPart("Hex")
        props = model_part.CreateNewProperties(1)
        for i, xyz in enumerate([
                (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]):
            model_part.CreateNewNode(i + 1, *xyz)
        model_part.CreateNewElement("Element3D8N", 1, [1, 2, 3, 4, 5, 6, 7, 8], props)

        _, edge_index, edge_features, _ = graph_bridge.BuildGraph(model_part)
        # True geometry edges only: 12 undirected, 24 directed, no diagonals.
        self.assertEqual(edge_index.shape, (2, 24))
        self.assertTrue(numpy.allclose(edge_features[:, 3], 1.0))  # all unit-length


if __name__ == '__main__':
    KratosUnittest.main()
