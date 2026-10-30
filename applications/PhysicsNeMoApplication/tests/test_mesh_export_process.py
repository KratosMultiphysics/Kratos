from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.processes.export import mesh_export_process
try:
    import physicsnemo.mesh
    from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import CreateMeshDataset
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestMeshExportProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.output_path = Path("test_mesh_export")
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        props = self.model_part.CreateNewProperties(1)
        for i, xyz in enumerate([
                (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]):
            node = self.model_part.CreateNewNode(i + 1, *xyz)
            node.SetSolutionStepValue(Kratos.PRESSURE, float(i))
        self.model_part.CreateNewElement("Element3D8N", 1, [1, 2, 3, 4, 5, 6, 7, 8], props)

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.output_path))

    def _Process(self):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "list_of_fields"  : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
                "output_path"     : "test_mesh_export",
                "output_interval" : 1
            }
        }""")
        process = mesh_export_process.Factory(settings, self.model)
        process.ExecuteInitialize()
        return process

    def test_TheTessellationIsReusedAcrossExports(self):
        """Only the field data changes between steps of a static mesh."""
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

        process = self._Process()
        calls = []
        original = domain_mesh_builder.BuildProvenance

        def Counting(*args, **kwargs):
            calls.append(1)
            return original(*args, **kwargs)

        domain_mesh_builder.BuildProvenance = Counting
        try:
            for step in (1, 2, 3):
                self.model_part.ProcessInfo[Kratos.STEP] = step
                process.ExecuteFinalizeSolutionStep()
        finally:
            domain_mesh_builder.BuildProvenance = original

        self.assertEqual(len(calls), 1, "the mesh was re-tessellated on a static mesh")
        self.assertEqual(
            sorted(path.name for path in self.output_path.glob("*.pmsh")),
            ["mesh_1.pmsh", "mesh_2.pmsh", "mesh_3.pmsh"])

    def test_MovedNodesAreExportedNotTheCachedGeometry(self):
        """A reused provenance map would write the old coordinates."""
        import numpy

        process = self._Process()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        for node in self.model_part.Nodes:
            node.X = node.X * 3.0
        self.model_part.ProcessInfo[Kratos.STEP] = 2
        process.ExecuteFinalizeSolutionStep()

        dataset = CreateMeshDataset(self.output_path)
        self.assertEqual(len(dataset), 2)
        points = [numpy.asarray(dataset[index][0].points) for index in (0, 1)]
        self.assertFalse(
            numpy.allclose(points[0], points[1]),
            "the export wrote cached geometry after the nodes moved")

    def test_InvalidCurvedRefinementRaises(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.export import mesh_export_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"          : "Main",
                "list_of_fields"           : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
                "higher_order_mode"        : "curved",
                "curved_refinement_levels" : 0
            }
        }""")
        with self.assertRaisesRegex(ValueError, "curved_refinement_levels"):
            mesh_export_process.Factory(settings, self.model)

    def test_ExportAndReadBack(self):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "list_of_fields"  : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
                "output_path"     : "test_mesh_export",
                "output_interval" : 1
            }
        }""")
        process = mesh_export_process.Factory(settings, self.model)
        process.ExecuteInitialize()
        for step in (1, 2):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            for node in self.model_part.Nodes:  # stand-in for a solve
                node.SetSolutionStepValue(Kratos.PRESSURE, float(node.Id * step))
            process.ExecuteFinalizeSolutionStep()

        exported = sorted(p.name for p in self.output_path.glob("*.pmsh"))
        self.assertEqual(exported, ["mesh_1.pmsh", "mesh_2.pmsh"])

        dataset = CreateMeshDataset(self.output_path)
        self.assertEqual(len(dataset), 2)
        mesh, _ = dataset[0]
        self.assertEqual(tuple(mesh.points.shape), (8, 3))
        self.assertEqual(tuple(mesh.cells.shape), (6, 4))  # hex -> 6 tets
        self.assertIn("PRESSURE", mesh.point_data)

    def test_IntervalGating(self):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "list_of_fields"  : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
                "output_path"     : "test_mesh_export",
                "output_interval" : 2
            }
        }""")
        process = mesh_export_process.Factory(settings, self.model)
        process.ExecuteInitialize()
        for step in (1, 2, 3):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
        self.assertEqual([p.name for p in sorted(self.output_path.glob("*.pmsh"))], ["mesh_2.pmsh"])

    def test_EmptyDirectoryRaises(self):
        self.output_path.mkdir(exist_ok=True)
        with self.assertRaisesRegex(ValueError, "No paths matching"):
            CreateMeshDataset(self.output_path)


if __name__ == '__main__':
    KratosUnittest.main()
