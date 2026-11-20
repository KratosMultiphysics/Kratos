import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder


class TestProvenanceRoundTrip(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test")
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        props = self.model_part.CreateNewProperties(1)

        # One unit hexahedron + one tetrahedron glued to its x=1 face corner
        # region -> mixed mesh.
        coordinates = [
            (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
            (2.0, 0.0, 0.0),
        ]
        for i, xyz in enumerate(coordinates):
            node = self.model_part.CreateNewNode(i + 1, *xyz)
            node.SetSolutionStepValue(Kratos.VELOCITY, [node.Id * 1.0, node.Id * 2.0, -1.0])
        self.model_part.CreateNewElement("Element3D8N", 1, [1, 2, 3, 4, 5, 6, 7, 8], props)
        self.model_part.CreateNewElement("Element3D4N", 2, [2, 9, 3, 6], props)

    def test_NodalFieldExactRoundTrip(self):
        provenance = domain_mesh_builder.BuildProvenance(self.model_part)
        node_ids = [node.Id for node in self.model_part.Nodes]

        ta = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY)
        ta.CollectData()
        original = numpy.array(ta.data)

        # Kratos -> simplex points -> identity "prediction" -> back to Kratos
        point_field = provenance.GatherNodalField(node_ids, original)
        domain_mesh_builder.ScatterFieldBack(
            provenance, point_field, self.model_part, Kratos.VELOCITY, "node_historical")

        ta.CollectData()
        self.assertTrue(numpy.array_equal(numpy.array(ta.data), original))

    def test_CellFieldAggregation(self):
        provenance = domain_mesh_builder.BuildProvenance(self.model_part)

        # Constant per-source-entity cell field must survive mean-aggregation
        # exactly: assign each simplex cell the id of its source element.
        cell_field = provenance.cell_provenance[:, 0].astype(float)
        domain_mesh_builder.ScatterFieldBack(
            provenance, cell_field, self.model_part, Kratos.PRESSURE, "element")

        for element in self.model_part.Elements:
            self.assertAlmostEqual(element.GetValue(Kratos.PRESSURE), float(element.Id), places=12)

    def test_CellFieldWeightedAggregation(self):
        provenance = domain_mesh_builder.BuildProvenance(self.model_part)

        cell_field = numpy.ones(provenance.number_of_cells)
        entity_ids, values = provenance.AggregateCellField(
            cell_field, "weighted_mean", provenance.ComputeSimplexMeasures())
        self.assertTrue(numpy.allclose(values, 1.0))
        self.assertEqual(sorted(entity_ids.tolist()), [1, 2])

    def test_ScatterShapeMismatchRaises(self):
        provenance = domain_mesh_builder.BuildProvenance(self.model_part)
        with self.assertRaisesRegex(ValueError, "rows but the tessellation"):
            domain_mesh_builder.ScatterFieldBack(
                provenance, numpy.zeros((3, 3)), self.model_part, Kratos.VELOCITY, "node_historical")

    def test_GaussPointScatterBackRefused(self):
        provenance = domain_mesh_builder.BuildProvenance(self.model_part)
        with self.assertRaisesRegex(ValueError, "Gauss-point"):
            domain_mesh_builder.ScatterFieldBack(
                provenance, numpy.zeros(provenance.number_of_cells),
                self.model_part, Kratos.PRESSURE, "element_gauss_point")


if __name__ == '__main__':
    KratosUnittest.main()
