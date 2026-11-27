from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

try:
    import physicsnemo.mesh
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestDomainMesh(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        props = self.model_part.CreateNewProperties(1)

        # One unit hex, with its z=0 face exposed as a boundary condition.
        for i, xyz in enumerate([
                (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]):
            node = self.model_part.CreateNewNode(i + 1, *xyz)
            node.SetSolutionStepValue(Kratos.PRESSURE, float(i))
        self.model_part.CreateNewElement("Element3D8N", 1, [1, 2, 3, 4, 5, 6, 7, 8], props)

        bottom = self.model_part.CreateSubModelPart("Bottom")
        bottom.AddNodes([1, 2, 3, 4])
        condition = self.model_part.CreateNewCondition("SurfaceCondition3D4N", 1, [1, 2, 3, 4], props)
        bottom.AddConditions([1])

        self.empty = self.model_part.CreateSubModelPart("Empty")

    def test_BuildDomainMeshWithBoundary(self):
        domain_mesh, provenance_maps = domain_mesh_builder.BuildDomainMesh(
            self.model_part,
            field_specs=[(Kratos.PRESSURE, "node_historical")],
            boundary_sub_model_part_names=["Bottom"])

        self.assertIn("interior", provenance_maps)
        self.assertIn("Bottom", provenance_maps)
        self.assertEqual(provenance_maps["interior"].number_of_cells, 6)   # hex -> 6 tets
        self.assertEqual(provenance_maps["Bottom"].number_of_cells, 2)     # quad -> 2 tris
        self.assertEqual(provenance_maps["Bottom"].number_of_points, 4)
        self.assertIn("PRESSURE", domain_mesh.interior.point_data)
        self.assertIn("PRESSURE", domain_mesh.boundaries["Bottom"].point_data)

    def test_EmptyBoundaryIsSkipped(self):
        domain_mesh, provenance_maps = domain_mesh_builder.BuildDomainMesh(
            self.model_part, boundary_sub_model_part_names=["Empty"])
        self.assertNotIn("Empty", provenance_maps)

    def test_SaveLoadRoundTrip(self):
        mesh, provenance = domain_mesh_builder.BuildMesh(
            self.model_part, field_specs=[(Kratos.PRESSURE, "node_historical")])
        prefix = Path("test_domain_mesh_saved")
        try:
            domain_mesh_builder.SaveMesh(mesh, prefix)
            loaded = domain_mesh_builder.LoadMesh(prefix)
            self.assertEqual(tuple(loaded.points.shape), (provenance.number_of_points, 3))
            self.assertIn("PRESSURE", loaded.point_data)
        finally:
            KratosUtilities.DeleteDirectoryIfExisting(str(prefix))


if __name__ == '__main__':
    KratosUnittest.main()
