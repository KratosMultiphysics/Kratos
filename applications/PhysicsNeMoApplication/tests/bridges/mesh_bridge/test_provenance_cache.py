"""ProvenanceCache: reuse the tessellation, never the geometry."""

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder


class TestProvenanceCache(KratosUnittest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Main")
        properties = self.model_part.CreateNewProperties(1)
        for i, xyz in enumerate([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                                 (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]):
            self.model_part.CreateNewNode(i + 1, *xyz)
        self.model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], properties)
        self.cache = domain_mesh_builder.ProvenanceCache()

    def _CountingBuild(self):
        calls = []
        original = domain_mesh_builder.BuildProvenance

        def Counting(*args, **kwargs):
            calls.append(1)
            return original(*args, **kwargs)
        return calls, original, Counting

    def test_AStaticMeshIsTessellatedOnce(self):
        calls, original, counting = self._CountingBuild()
        domain_mesh_builder.BuildProvenance = counting
        try:
            first = self.cache.Get(self.model_part)
            for _ in range(4):
                self.assertIs(self.cache.Get(self.model_part), first)
        finally:
            domain_mesh_builder.BuildProvenance = original
        self.assertEqual(len(calls), 1)

    def test_MovingANodeRebuildsAndThePointsFollow(self):
        """The map carries simplex_points, so reuse across a moving mesh
        would hand back stale geometry - the reason the guard compares
        coordinates rather than only the entity count."""
        before = self.cache.Get(self.model_part).simplex_points.copy()

        node = self.model_part.GetNode(2)
        node.X = node.X + 0.5
        after = self.cache.Get(self.model_part).simplex_points

        self.assertFalse(numpy.allclose(before, after),
                         "the cache served stale geometry after a node moved")
        coordinates = numpy.array([[n.X, n.Y, n.Z] for n in self.model_part.Nodes])
        for point in coordinates:
            self.assertTrue(numpy.isclose(after, point).all(axis=1).any(),
                            "a current node position is missing from the rebuilt points")

    def test_ATinyMoveStillInvalidates(self):
        before = self.cache.Get(self.model_part)
        self.model_part.GetNode(1).X += 1e-9
        self.assertIsNot(self.cache.Get(self.model_part), before)

    def test_AddedEntitiesRebuild(self):
        before = self.cache.Get(self.model_part)
        properties = self.model_part.GetProperties()[1]
        self.model_part.CreateNewNode(5, 1.0, 1.0, 1.0)
        self.model_part.CreateNewElement("Element3D4N", 2, [2, 3, 4, 5], properties)
        after = self.cache.Get(self.model_part)
        self.assertIsNot(after, before)
        self.assertGreater(after.number_of_cells, before.number_of_cells)


if __name__ == "__main__":
    KratosUnittest.main()
