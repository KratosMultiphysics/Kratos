from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.processes.export import cae_dataset_export_process
try:
    import torch  # noqa: F401
    import warp  # noqa: F401 - DoMINO preprocessing needs it (SDF)
    import physicsnemo.datapipes.cae.cae_dataset  # noqa: F401
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


def CreateCaeFixture(model):
    """Unit cube: 1 hexahedron element + 'Skin' sub-part with 6 quad faces."""
    model_part = model.CreateModelPart("Main")
    model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
    model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
    props = model_part.CreateNewProperties(1)
    corners = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
               (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]
    for i, xyz in enumerate(corners):
        node = model_part.CreateNewNode(i + 1, *xyz)
        node.SetSolutionStepValue(Kratos.PRESSURE, 1.0 + node.X + 2.0 * node.Y)
        node.SetSolutionStepValue(Kratos.DISPLACEMENT, [node.X, node.Y, node.Z])
    model_part.CreateNewElement("Element3D8N", 1, list(range(1, 9)), props)

    skin = model_part.CreateSubModelPart("Skin")
    faces = [(1, 4, 3, 2), (5, 6, 7, 8), (1, 2, 6, 5), (2, 3, 7, 6), (3, 4, 8, 7), (4, 1, 5, 8)]
    for i, face in enumerate(faces):
        model_part.CreateNewCondition("SurfaceCondition3D4N", i + 1, list(face), props)
    skin.AddNodes(list(range(1, 9)))
    skin.AddConditions(list(range(1, 7)))
    return model_part


def CreateExportProcess(model, output_path, extra=""):
    settings = Kratos.Parameters("""{
        "Parameters": {
            "model_part_name"         : "Main",
            "surface_model_part_name" : "Main.Skin",
            "surface_fields"          : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" },
                                          { "variable_name" : "DISPLACEMENT", "data_location" : "node_historical" } ],
            "volume_fields"           : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
            "global_params"           : { "stream_velocity" : 30.0, "air_density" : 1.226 },
            "global_params_order"     : ["stream_velocity", "air_density"],
            "output_path"             : "%s"
            %s
        }
    }""" % (output_path, extra))
    return cae_dataset_export_process.Factory(settings, model)


class TestCaeDatasetExportProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.output_path = Path("test_cae_dataset_export")
        self.model = Kratos.Model()
        self.model_part = CreateCaeFixture(self.model)

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.output_path))

    def _Export(self, extra="", step=1):
        process = CreateExportProcess(self.model, self.output_path, extra)
        process.ExecuteInitialize()
        self.model_part.ProcessInfo[Kratos.STEP] = step
        self.model_part.ProcessInfo[Kratos.TIME] = float(step)
        process.ExecuteFinalizeSolutionStep()
        return process

    def test_ExportedSupersetLayout(self):
        self._Export()
        with numpy.load(self.output_path / "case_1.npz") as data:
            arrays = {key: numpy.array(data[key]) for key in data.files}

        for key in ("stl_coordinates", "stl_faces", "stl_centers", "stl_areas",
                    "surface_mesh_centers", "surface_normals", "surface_areas",
                    "surface_fields", "volume_mesh_centers", "volume_fields",
                    "global_params_values", "global_params_reference",
                    "stream_velocity", "air_density", "TIME", "STEP"):
            self.assertIn(key, arrays)

        # 6 quads -> 12 triangles; faces flattened int32
        self.assertEqual(arrays["stl_faces"].shape, (36,))
        self.assertEqual(arrays["stl_faces"].dtype, numpy.int32)
        self.assertEqual(arrays["stl_centers"].shape, (12, 3))
        self.assertAlmostEqual(float(arrays["stl_areas"].sum()), 6.0, places=6)
        numpy.testing.assert_allclose(
            numpy.linalg.norm(arrays["surface_normals"], axis=1), 1.0, atol=1e-6)

        # fields: PRESSURE (1) + DISPLACEMENT (3) on 12 triangles; volume on 8 nodes
        self.assertEqual(arrays["surface_fields"].shape, (12, 4))
        self.assertEqual(arrays["volume_mesh_centers"].shape, (8, 3))
        self.assertEqual(arrays["volume_fields"].shape, (8, 1))
        # nodal linear field -> triangle value equals field at the center
        centers = arrays["stl_centers"]
        numpy.testing.assert_allclose(
            arrays["surface_fields"][:, 0], 1.0 + centers[:, 0] + 2.0 * centers[:, 1], atol=1e-5)

        # globals: both spellings, DoMINO stack in insertion order
        numpy.testing.assert_allclose(arrays["global_params_values"], [[30.0], [1.226]])
        numpy.testing.assert_allclose(arrays["global_params_reference"], [[30.0], [1.226]])
        numpy.testing.assert_allclose(arrays["stream_velocity"], [30.0])

        # regression guard: the npz reader chokes on 0-d arrays
        for key, value in arrays.items():
            self.assertGreaterEqual(value.ndim, 1, f"{key} must not be 0-d")

    def test_TheTessellationIsReusedAcrossExports(self):
        """Tessellation is ~99% of an export and is identical while the mesh is."""
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

        process = CreateExportProcess(self.model, self.output_path, ', "output_interval": 1')
        process.ExecuteInitialize()

        calls = []
        original = domain_mesh_builder.BuildProvenance

        def Counting(*args, **kwargs):
            calls.append(1)
            return original(*args, **kwargs)

        domain_mesh_builder.BuildProvenance = Counting
        try:
            written = []
            for step in (1, 2, 3):
                self.model_part.ProcessInfo[Kratos.STEP] = step
                self.model_part.ProcessInfo[Kratos.TIME] = float(step)
                process.ExecuteFinalizeSolutionStep()
                with numpy.load(self.output_path / f"case_{step}.npz") as data:
                    written.append(numpy.array(data["stl_coordinates"]))
        finally:
            domain_mesh_builder.BuildProvenance = original

        self.assertEqual(len(calls), 1, "the surface was re-tessellated on a static mesh")
        numpy.testing.assert_allclose(written[1], written[0], rtol=0.0, atol=0.0)
        numpy.testing.assert_allclose(written[2], written[0], rtol=0.0, atol=0.0)

    def test_MovedNodesAreExportedNotTheCachedGeometry(self):
        """The provenance map carries simplex_points.

        Reusing one across a moving mesh would write the old coordinates into
        every later case - a silent wrong answer, which is why the cache is
        invalidated by the coordinates and not just by the entity count.
        """
        process = CreateExportProcess(self.model, self.output_path, ', "output_interval": 1')
        process.ExecuteInitialize()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        self.model_part.ProcessInfo[Kratos.TIME] = 1.0
        process.ExecuteFinalizeSolutionStep()
        with numpy.load(self.output_path / "case_1.npz") as data:
            before = numpy.array(data["stl_coordinates"])

        for node in self.model_part.Nodes:
            node.X = node.X * 2.0
        self.model_part.ProcessInfo[Kratos.STEP] = 2
        self.model_part.ProcessInfo[Kratos.TIME] = 2.0
        process.ExecuteFinalizeSolutionStep()
        with numpy.load(self.output_path / "case_2.npz") as data:
            after = numpy.array(data["stl_coordinates"])

        self.assertFalse(
            numpy.allclose(before, after),
            "the export wrote cached geometry after the nodes moved")

    def test_CaseIdNamingAndInterval(self):
        self._Export(extra=', "case_id": 7, "output_interval": 2', step=2)
        self.assertTrue((self.output_path / "case_7.npz").is_file())
        # step 3 is not due
        self.model_part.ProcessInfo[Kratos.STEP] = 3
        CreateExportProcess(self.model, self.output_path, ', "case_id": 8, "output_interval": 2') \
            .ExecuteFinalizeSolutionStep()
        self.assertFalse((self.output_path / "case_8.npz").is_file())

    def test_EmptyFieldListsOmitKeys(self):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"         : "Main",
                "surface_model_part_name" : "Main.Skin",
                "output_path"             : "%s"
            }
        }""" % self.output_path)
        process = cae_dataset_export_process.Factory(settings, self.model)
        process.ExecuteInitialize()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        with numpy.load(self.output_path / "case_1.npz") as data:
            self.assertNotIn("surface_fields", data.files)
            self.assertNotIn("volume_fields", data.files)
            self.assertNotIn("global_params_values", data.files)
            self.assertIn("stl_centers", data.files)

    def test_CurvedSurfaceExport(self):
        # Quadrilateral3D9 skin exported in curved mode: the STL carries the
        # synthetic lattice points (per face: (2^k+1)^2 = 25 points at k=2)
        model = Kratos.Model()
        model_part = model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        props = model_part.CreateNewProperties(1)
        positions = [(-1, -1), (1, -1), (1, 1), (-1, 1), (0, -1), (1, 0), (0, 1), (-1, 0), (0, 0)]
        for i, (x, y) in enumerate(positions):
            node = model_part.CreateNewNode(i + 1, float(x), float(y), 0.05 * (1 - x * x) * (1 - y * y))
            node.SetSolutionStepValue(Kratos.PRESSURE, 1.0 + node.X)
        skin = model_part.CreateSubModelPart("Skin")
        model_part.CreateNewCondition("SurfaceCondition3D9N", 1, list(range(1, 10)), props)
        skin.AddNodes(list(range(1, 10)))
        skin.AddConditions([1])

        process = cae_dataset_export_process.Factory(Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"          : "Main",
                "surface_model_part_name"  : "Main.Skin",
                "surface_fields"           : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
                "higher_order_mode"        : "curved",
                "curved_refinement_levels" : 2,
                "output_path"              : "%s"
            }
        }""" % self.output_path), model)
        process.ExecuteInitialize()
        model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        with numpy.load(self.output_path / "case_1.npz") as data:
            self.assertEqual(numpy.array(data["stl_coordinates"]).shape, (25, 3))
            self.assertEqual(numpy.array(data["stl_faces"]).shape, (32 * 3,))
            self.assertEqual(numpy.array(data["surface_fields"]).shape, (32, 1))
            self.assertTrue(numpy.isfinite(numpy.array(data["surface_fields"])).all())

    def test_InvalidSettingsRaise(self):
        with self.assertRaisesRegex(ValueError, "matching global_params entry"):
            cae_dataset_export_process.Factory(Kratos.Parameters("""{
                "Parameters": {
                    "model_part_name"         : "Main",
                    "surface_model_part_name" : "Main.Skin",
                    "global_params_reference" : { "orphan" : 1.0 }
                }
            }"""), self.model)
        with self.assertRaisesRegex(ValueError, "Volume field"):
            cae_dataset_export_process.Factory(Kratos.Parameters("""{
                "Parameters": {
                    "model_part_name"         : "Main",
                    "surface_model_part_name" : "Main.Skin",
                    "volume_fields"           : [ { "variable_name" : "PRESSURE", "data_location" : "element" } ]
                }
            }"""), self.model)


@KratosUnittest.skipUnless(have_physicsnemo,
                           "Missing required python modules: torch, warp, physicsnemo.")
class TestCaeDatasetThroughDatapipes(KratosUnittest.TestCase):
    """The exported directory must load through the real physicsnemo pipes."""

    def setUp(self):
        self.output_path = Path("test_cae_datapipes")
        model = Kratos.Model()
        model_part = CreateCaeFixture(model)
        process = CreateExportProcess(model, self.output_path)
        process.ExecuteInitialize()
        model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.output_path))

    def test_DoMINOSurfaceAndVolume(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import CreateDoMINODataPipe

        for model_type in ("surface", "volume"):
            pipe = CreateDoMINODataPipe(
                self.output_path, model_type,
                bounding_box=((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
                num_surface_neighbors=5)
            self.assertEqual(len(pipe), 1)
            sample = pipe[0]
            self.assertIn("geometry_coordinates", sample)
            if model_type == "surface":
                self.assertIn("surface_mesh_centers", sample)
                self.assertIn("surface_mesh_neighbors", sample)  # computed internally by kNN
                self.assertIn("surface_fields", sample)
            else:
                self.assertIn("volume_mesh_centers", sample)
                self.assertIn("volume_fields", sample)

    def test_TransolverSurfaceAndVolume(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import CreateTransolverDataPipe

        for model_type in ("surface", "volume"):
            pipe = CreateTransolverDataPipe(self.output_path, model_type)
            self.assertEqual(len(pipe), 1)
            sample = pipe[0]
            self.assertIn("embeddings", sample)
            self.assertIn("fields", sample)
            self.assertIn("fx", sample)  # stream_velocity/air_density globals
            for value in sample.values():
                if hasattr(value, "isfinite"):
                    self.assertTrue(bool(value.isfinite().all()))


if __name__ == '__main__':
    KratosUnittest.main()
