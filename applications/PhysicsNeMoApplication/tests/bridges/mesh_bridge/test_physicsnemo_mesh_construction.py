import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

try:
    import physicsnemo.mesh
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestPhysicsNemoMeshConstruction(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test")
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        props = self.model_part.CreateNewProperties(1)
        for i, xyz in enumerate([
                (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]):
            node = self.model_part.CreateNewNode(i + 1, *xyz)
            node.SetSolutionStepValue(Kratos.PRESSURE, float(i))
        self.model_part.CreateNewElement("Element3D8N", 1, [1, 2, 3, 4, 5, 6, 7, 8], props)

    def test_MeshConstructionFromTessellation(self):
        mesh, provenance = domain_mesh_builder.BuildMesh(
            self.model_part, field_specs=[(Kratos.PRESSURE, "node_historical")])

        self.assertEqual(len(provenance.simplex_points), 8)
        self.assertEqual(len(provenance.simplex_cells), 6)
        self.assertIn("PRESSURE", mesh.point_data)


if __name__ == '__main__':
    KratosUnittest.main()
