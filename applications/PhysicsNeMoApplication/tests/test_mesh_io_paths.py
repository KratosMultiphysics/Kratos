"""Tests for the two file-format paths added around the mesh bridge: Zarr
written by physicsnemo itself (no physicsnemo-curator install), and Kratos's
own VTK output read back through physicsnemo's VTKReader."""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import physicsnemo.mesh.io  # noqa: F401
    import zarr  # noqa: F401
    have_zarr_io = True
except ImportError:
    have_zarr_io = False

try:
    import pyvista  # noqa: F401
    from physicsnemo.datapipes.readers import VTKReader  # noqa: F401
    have_vtk = True
except ImportError:
    have_vtk = False


def _ModelPartWithAField(model, name="Main", divisions=3):
    model_part = CreateStructuredTetModelPart(
        model, name, divisions=divisions, historical_variables=(Kratos.PRESSURE,))
    for node in model_part.Nodes:
        node.SetSolutionStepValue(Kratos.PRESSURE, 2.0 * node.X + node.Y)
    return model_part


@KratosUnittest.skipUnless(have_torch and have_zarr_io,
                           "Missing required python modules: torch, physicsnemo (>= 2.2), zarr.")
class TestZarrMeshIo(KratosUnittest.TestCase):
    """physicsnemo 2.2 writes Zarr itself, so the AI-ready store no longer
    needs physicsnemo-curator - a git-only package that builds a Rust
    toolchain at install time and whose sinks reject a DomainMesh."""

    def setUp(self):
        self.store = Path("test_mesh_store.zarr")
        self.model = Kratos.Model()
        self.model_part = _ModelPartWithAField(self.model)

    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting(str(self.store))

    def test_MeshRoundTrip(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
            domain_mesh_builder)

        mesh, _ = domain_mesh_builder.BuildMesh(
            self.model_part, field_specs=[(Kratos.PRESSURE, "node_historical")])
        domain_mesh_builder.SaveMeshZarr(mesh, self.store)
        self.assertTrue(self.store.exists())

        reloaded = domain_mesh_builder.LoadMeshZarr(self.store)
        numpy.testing.assert_allclose(
            reloaded.points.numpy(), mesh.points.numpy(), atol=1e-12)
        numpy.testing.assert_array_equal(reloaded.cells.numpy(), mesh.cells.numpy())
        numpy.testing.assert_allclose(
            reloaded.point_data["PRESSURE"].numpy(),
            mesh.point_data["PRESSURE"].numpy(), atol=1e-12)

    def test_DomainMeshRoundTripKeepsItsBoundaryNames(self):
        """The case curator's own sinks cannot handle: a DomainMesh is
        Mesh-typed there and fails on point_data."""
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
            domain_mesh_builder)

        boundary = self.model_part.CreateSubModelPart("Inlet")
        face_nodes = sorted(node.Id for node in self.model_part.Nodes if node.X < 1e-12)
        boundary.AddNodes(face_nodes)
        properties = self.model_part.GetProperties()[1]
        for index in range(len(face_nodes) - 2):
            boundary.CreateNewCondition(
                "SurfaceCondition3D3N", index + 1,
                [face_nodes[index], face_nodes[index + 1], face_nodes[index + 2]],
                properties)

        domain_mesh, _ = domain_mesh_builder.BuildDomainMesh(
            self.model_part, field_specs=[(Kratos.PRESSURE, "node_historical")],
            boundary_sub_model_part_names=("Inlet",))
        domain_mesh_builder.SaveMeshZarr(domain_mesh, self.store)
        reloaded = domain_mesh_builder.LoadMeshZarr(self.store)
        self.assertIn("Inlet", reloaded.boundaries)
        numpy.testing.assert_allclose(
            reloaded.interior.points.numpy(), domain_mesh.interior.points.numpy(),
            atol=1e-12)


@KratosUnittest.skipUnless(have_torch and have_zarr_io,
                           "Missing required python modules: torch, physicsnemo (>= 2.2), zarr.")
class TestCuratorExportProcessZarrBackend(KratosUnittest.TestCase):
    def setUp(self):
        self.output = Path("test_curator_free_zarr")
        self.model = Kratos.Model()
        self.model_part = _ModelPartWithAField(self.model)

    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting(str(self.output))

    def test_AStorePerStepWithoutCurator(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.export import (
            curator_export_process)
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
            domain_mesh_builder)

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "list_of_fields"  : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
                "output_path"     : "%s",
                "file_prefix"     : "mesh"
            }
        }""" % self.output)
        process = curator_export_process.Factory(settings, self.model)
        process.ExecuteInitialize()
        for step in (1, 2):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
        process.ExecuteFinalize()

        stores = sorted(self.output.glob("mesh_*.zarr"))
        self.assertEqual(len(stores), 2)
        mesh = domain_mesh_builder.LoadMeshZarr(stores[0])
        self.assertIn("PRESSURE", mesh.point_data)
        self.assertEqual(mesh.points.shape[0], self.model_part.NumberOfNodes())

    def test_AnUnknownBackendRaises(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.export import (
            curator_export_process)

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "list_of_fields"  : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
                "zarr_backend"    : "parquet"
            }
        }""")
        with self.assertRaisesRegex(ValueError, "zarr_backend"):
            curator_export_process.Factory(settings, self.model)


@KratosUnittest.skipUnless(have_torch and have_vtk,
                           "Missing required python modules: torch, pyvista, physicsnemo (>= 2.2).")
class TestVtkBridge(KratosUnittest.TestCase):
    """Kratos's own VtkOutputProcess output as training data. Read-only:
    the files carry no provenance, so nothing scatters back."""

    def setUp(self):
        self.output = Path("test_vtk_output")
        self.arranged = Path("test_vtk_samples")
        self.model = Kratos.Model()
        self.model_part = _ModelPartWithAField(self.model)
        self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

    def tearDown(self):
        for directory in (self.output, self.arranged):
            kratos_utils.DeleteDirectoryIfExisting(str(directory))

    def _WriteVtkOutput(self, steps=(1, 2)):
        import KratosMultiphysics.vtk_output_process as vtk_output_process

        vtk_settings = Kratos.Parameters("""{
            "Parameters" : {
                "model_part_name"                    : "Main",
                "output_path"                        : "%s",
                "file_format"                        : "ascii",
                "output_control_type"                : "step",
                "output_interval"                    : 1,
                "nodal_solution_step_data_variables" : ["PRESSURE"]
            }
        }""" % self.output)
        process = vtk_output_process.Factory(vtk_settings, self.model)
        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()
        for step in steps:
            self.model_part.ProcessInfo[Kratos.STEP] = step
            self.model_part.ProcessInfo[Kratos.TIME] = float(step)
            process.ExecuteInitializeSolutionStep()
            process.ExecuteFinalizeSolutionStep()
            process.PrintOutput()
        process.ExecuteFinalize()

    def test_TheSolversOwnOutputBecomesAReaderDataset(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import vtk_bridge

        self._WriteVtkOutput()
        count = vtk_bridge.ArrangeVtkOutputForReader(self.output, self.arranged)
        self.assertGreaterEqual(count, 2)
        self.assertEqual(len(sorted(self.arranged.glob("sample_*"))), count)

        reader = vtk_bridge.CreateVtkReaderDataset(self.arranged)
        self.assertEqual(len(reader), count)

        # ... but it reads NOTHING from them, because its vocabulary is a
        # fixed external-aero one and a field called PRESSURE is in none of
        # its key sets. No error, just an empty sample.
        data, metadata = reader[0]
        self.assertEqual(list(data.keys()), [])
        self.assertIn("source_file", metadata)

    def test_TheGeneralPathReadsTheFieldsByTheirOwnNames(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import vtk_bridge

        self._WriteVtkOutput()
        vtk_bridge.ArrangeVtkOutputForReader(self.output, self.arranged)
        dataset = vtk_bridge.CreateVtkMeshDataset(self.arranged)
        self.assertEqual(len(dataset), 2)
        mesh = dataset[0]
        self.assertIn("PRESSURE", mesh.point_data)
        self.assertEqual(mesh.points.shape[0], self.model_part.NumberOfNodes())

    def test_AFileBecomesAMeshWithItsField(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import vtk_bridge

        self._WriteVtkOutput(steps=(1,))
        files = sorted(self.output.glob("*.vtk"))  # Kratos writes LEGACY vtk
        self.assertTrue(files)
        mesh = vtk_bridge.MeshFromVtkFile(files[0])
        self.assertEqual(mesh.points.shape[0], self.model_part.NumberOfNodes())
        self.assertIn("PRESSURE", mesh.point_data)

        # the field survives the round trip through the file
        written = numpy.array(sorted(
            node.GetSolutionStepValue(Kratos.PRESSURE) for node in self.model_part.Nodes))
        read_back = numpy.sort(numpy.asarray(
            mesh.point_data["PRESSURE"].numpy()).ravel())
        numpy.testing.assert_allclose(read_back, written, atol=1e-6)

    def test_KratosWritesLegacyVtkWhichTheReaderCannotRead(self):
        """The mismatch the arrange step exists for: pointing the reader at
        a Kratos output directory finds nothing, because VtkOutputProcess
        writes .vtk and VTKReader reads .vtu/.vtp/.stl."""
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import vtk_bridge

        self._WriteVtkOutput(steps=(1,))
        self.assertTrue(sorted(self.output.glob("*.vtk")))
        self.assertFalse(sorted(self.output.glob("*.vtu")))

        vtk_bridge.ArrangeVtkOutputForReader(self.output, self.arranged)
        converted = sorted(self.arranged.glob("sample_*/*.vtu"))
        self.assertTrue(converted)

    def test_AnEmptyDirectoryIsReported(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import vtk_bridge

        self.output.mkdir(parents=True, exist_ok=True)
        with self.assertRaisesRegex(ValueError, "No files matching"):
            vtk_bridge.ArrangeVtkOutputForReader(self.output, self.arranged)


if __name__ == '__main__':
    KratosUnittest.main()
