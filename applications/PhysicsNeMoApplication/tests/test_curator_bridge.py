from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.bridges import curator_bridge
from KratosMultiphysics.PhysicsNeMoApplication.processes.export import curator_export_process
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

try:
    import physicsnemo_curator.core.base  # noqa: F401
    have_curator = True
except ImportError:
    have_curator = False

try:
    # The mesh sinks need curator's "mesh" extra (pyarrow), so they fail
    # independently of the core source-side import above.
    import physicsnemo_curator.domains.mesh.sinks.mesh_zarr  # noqa: F401
    have_curator_sinks = True
except ImportError:
    have_curator_sinks = False


def CreateCubeFixture(model, name, z_offset=0.0):
    """Unit cube: 1 hexahedron with a linear nodal PRESSURE field."""
    model_part = model.CreateModelPart(name)
    model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
    props = model_part.CreateNewProperties(1)
    corners = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
               (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]
    for i, (x, y, z) in enumerate(corners):
        node = model_part.CreateNewNode(i + 1, x, y, z + z_offset)
        node.SetSolutionStepValue(Kratos.PRESSURE, 1.0 + node.X + 2.0 * node.Y)
    model_part.CreateNewElement("Element3D8N", 1, list(range(1, 9)), props)
    return model_part


@KratosUnittest.skipUnless(have_curator, "Missing required python module: physicsnemo_curator.")
class TestCuratorSource(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.parts = [CreateCubeFixture(self.model, "A"),
                      CreateCubeFixture(self.model, "B", z_offset=0.5)]
        self.field_specs = [(Kratos.PRESSURE, "node_historical")]

    def test_SourceContract(self):
        source = curator_bridge.CreateKratosMeshSource(self.parts, self.field_specs)
        self.assertEqual(len(source), 2)
        self.assertEqual(type(source).name, "Kratos Mesh")
        self.assertTrue(type(source).description)
        self.assertEqual([p.name for p in type(source).params()],
                         ["source_container", "tessellation_mode",
                          "higher_order_mode", "curved_refinement_levels"])

    def test_ItemMatchesBuildMesh(self):
        source = curator_bridge.CreateKratosMeshSource(self.parts, self.field_specs)
        mesh = next(source[0])
        reference, _ = domain_mesh_builder.BuildMesh(self.parts[0], self.field_specs)
        self.assertVectorAlmostEqual(
            numpy.asarray(mesh.points).ravel(), numpy.asarray(reference.points).ravel(), 12)
        self.assertEqual(numpy.asarray(mesh.cells).tolist(),
                         numpy.asarray(reference.cells).tolist())
        points = numpy.asarray(mesh.points)
        self.assertVectorAlmostEqual(
            numpy.asarray(mesh.point_data["PRESSURE"]).ravel(),
            1.0 + points[:, 0] + 2.0 * points[:, 1], 12)

    def test_DistinctItemsPerIndex(self):
        # The per-source provenance cache must not leak item 0's geometry
        # into item 1 (it invalidates on node coordinates).
        source = curator_bridge.CreateKratosMeshSource(self.parts, self.field_specs)
        self.assertAlmostEqual(float(numpy.asarray(next(source[0]).points)[:, 2].min()), 0.0, 12)
        self.assertAlmostEqual(float(numpy.asarray(next(source[1]).points)[:, 2].min()), 0.5, 12)

    def test_LazyCallableForm(self):
        source = curator_bridge.CreateKratosMeshSource(
            (lambda index: self.parts[index], 2), self.field_specs)
        self.assertEqual(len(source), 2)
        self.assertEqual(tuple(next(source[1]).points.shape), (8, 3))

    def test_OutOfRangeRaises(self):
        source = curator_bridge.CreateKratosMeshSource(self.parts, self.field_specs)
        with self.assertRaisesRegex(IndexError, "out of range"):
            next(source[2])


@KratosUnittest.skipUnless(have_curator_sinks,
                           "Missing required python module: physicsnemo_curator mesh sinks.")
class TestCuratorPipeline(KratosUnittest.TestCase):
    def setUp(self):
        self.output_path = Path("test_curator_pipeline")
        self.model = Kratos.Model()
        self.parts = [CreateCubeFixture(self.model, "A"),
                      CreateCubeFixture(self.model, "B", z_offset=0.5)]
        self.field_specs = [(Kratos.PRESSURE, "node_historical")]

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.output_path))

    def test_ZarrRoundTrip(self):
        import zarr
        source = curator_bridge.CreateKratosMeshSource(self.parts, self.field_specs)
        summary = curator_bridge.RunCuratorPipeline(
            source, curator_bridge.CreateZarrSink(self.output_path))
        self.assertEqual(summary["completed"], 2)
        self.assertEqual(summary["failed"], 0)

        stores = sorted(self.output_path.glob("*.zarr"))
        self.assertEqual(len(stores), 2)
        group = zarr.open(str(stores[0]), mode="r")
        # The sink stores positions as mesh_pos with a leading time axis.
        positions = numpy.asarray(group["mesh_pos"])
        self.assertEqual(positions.shape, (1, 8, 3))
        pressure = numpy.asarray(group["PRESSURE"])
        # written as float32, so compare at single precision
        self.assertVectorAlmostEqual(
            pressure, 1.0 + positions[0][:, 0] + 2.0 * positions[0][:, 1], 5)

    def test_VtuRoundTrip(self):
        import pyvista
        source = curator_bridge.CreateKratosMeshSource(self.parts, self.field_specs)
        curator_bridge.RunCuratorPipeline(
            source, curator_bridge.CreateVtuSink(self.output_path))

        files = sorted(self.output_path.glob("*.vtu"))
        self.assertEqual(len(files), 2)
        grid = pyvista.read(str(files[0]))
        self.assertEqual(grid.n_points, 8)
        self.assertEqual(grid.n_cells, 6)          # hexahedron -> 6 tetrahedra
        points = numpy.asarray(grid.points)
        self.assertVectorAlmostEqual(
            numpy.asarray(grid.point_data["PRESSURE"]),
            1.0 + points[:, 0] + 2.0 * points[:, 1], 5)
        # second item carries the offset geometry, not a repeat of the first
        self.assertAlmostEqual(
            float(numpy.asarray(pyvista.read(str(files[1])).points)[:, 2].min()), 0.5, 5)

    def test_CurvedTessellationReachesSink(self):
        import pyvista
        model = Kratos.Model()
        model_part = model.CreateModelPart("Curved")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        props = model_part.CreateNewProperties(1)
        positions = [(-1, -1), (1, -1), (1, 1), (-1, 1), (0, -1), (1, 0), (0, 1), (-1, 0), (0, 0)]
        for i, (x, y) in enumerate(positions):
            node = model_part.CreateNewNode(i + 1, float(x), float(y),
                                            0.05 * (1 - x * x) * (1 - y * y))
            node.SetSolutionStepValue(Kratos.PRESSURE, 1.0 + node.X)
        model_part.CreateNewCondition("SurfaceCondition3D9N", 1, list(range(1, 10)), props)

        source = curator_bridge.CreateKratosMeshSource(
            [model_part], [(Kratos.PRESSURE, "node_historical")],
            source_container="Conditions", higher_order_mode="curved",
            curved_refinement_levels=2)
        curator_bridge.RunCuratorPipeline(
            source, curator_bridge.CreateVtuSink(self.output_path))
        grid = pyvista.read(str(sorted(self.output_path.glob("*.vtu"))[0]))
        # (2^k + 1)^2 = 25 lattice points, 2 * 4^k = 32 triangles at k = 2
        self.assertEqual(grid.n_points, 25)
        self.assertEqual(grid.n_cells, 32)


@KratosUnittest.skipUnless(have_curator_sinks,
                           "Missing required python module: physicsnemo_curator mesh sinks.")
class TestCuratorExportProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.output_path = Path("test_curator_export")
        self.model = Kratos.Model()
        self.model_part = CreateCubeFixture(self.model, "Main")

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.output_path))

    def _Process(self, sink="zarr", output_interval=1):
        return curator_export_process.Factory(Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "list_of_fields"  : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
                "sink"            : "%s",
                "output_path"     : "test_curator_export",
                "output_interval" : %d
            }
        }""" % (sink, output_interval)), self.model)

    def test_ZarrExportPerStep(self):
        process = self._Process()
        process.ExecuteInitialize()
        for step in (1, 2, 3):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
        self.assertEqual(len(sorted(self.output_path.glob("*.zarr"))), 3)

    def test_VtuExportPerStep(self):
        process = self._Process(sink="vtu")
        process.ExecuteInitialize()
        for step in (1, 2):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
        names = sorted(p.name for p in self.output_path.glob("*.vtu"))
        self.assertEqual(names, ["mesh_0001.vtu", "mesh_0002.vtu"])

    def test_IntervalGating(self):
        process = self._Process(output_interval=2)
        process.ExecuteInitialize()
        for step in (1, 2, 3, 4):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
        names = sorted(p.name for p in self.output_path.glob("*.zarr"))
        self.assertEqual(names, ["mesh_0002.zarr", "mesh_0004.zarr"])

    def test_InvalidSettingsRaise(self):
        with self.assertRaisesRegex(ValueError, "Unsupported \"sink\""):
            self._Process(sink="vtp")
        with self.assertRaisesRegex(ValueError, "output_interval"):
            self._Process(output_interval=0)


if __name__ == '__main__':
    KratosUnittest.main()
