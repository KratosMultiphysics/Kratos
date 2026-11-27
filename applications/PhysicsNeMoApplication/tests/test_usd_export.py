from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from test_curator_bridge import CreateCubeFixture

try:
    import pxr  # noqa: F401 - usd-core
    have_usd = True
except ImportError:
    have_usd = False

try:
    import torch  # noqa: F401
    import physicsnemo.mesh  # noqa: F401
    have_mesh_bridge = True
except ImportError:
    have_mesh_bridge = False


@KratosUnittest.skipUnless(have_usd, "Missing required python module: pxr (usd-core).")
class TestUsdExportCore(KratosUnittest.TestCase):
    """deployment.usd_export against a real pxr round trip - array in, USD out."""

    def setUp(self):
        self.stage_file = Path("test_usd_core.usda")

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.stage_file))

    @staticmethod
    def _Tetrahedron():
        points = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                              [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        triangles = numpy.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
        return points, triangles

    def test_MeshTimeSamplesRoundTrip(self):
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import usd_export
        from pxr import Usd, UsdGeom

        points, triangles = self._Tetrahedron()
        stage = usd_export.CreateUsdStage(self.stage_file, time_codes_per_second=2.0)
        usd_export.WriteMeshTimeSample(
            stage, "/Kratos/Main", points, triangles,
            {"PRESSURE": points[:, 0], "VELOCITY": 1.5 * points,
             "PAIR": numpy.ones((4, 2))}, 1.0)
        usd_export.WriteMeshTimeSample(
            stage, "/Kratos/Main", points + 0.1, triangles,
            {"PRESSURE": 10.0 * points[:, 0]}, 2.0)
        usd_export.SaveStage(stage)

        reopened = Usd.Stage.Open(str(self.stage_file))
        self.assertEqual(reopened.GetStartTimeCode(), 1.0)
        self.assertEqual(reopened.GetEndTimeCode(), 2.0)
        self.assertEqual(reopened.GetTimeCodesPerSecond(), 2.0)
        mesh = UsdGeom.Mesh(reopened.GetPrimAtPath("/Kratos/Main"))
        self.assertTrue(mesh)
        self.assertAlmostEqual(mesh.GetPointsAttr().Get(2.0)[0][0], 0.1, places=6)

        primvars = UsdGeom.PrimvarsAPI(mesh.GetPrim())
        pressure = primvars.GetPrimvar("PRESSURE")
        self.assertEqual(pressure.GetInterpolation(), UsdGeom.Tokens.vertex)
        self.assertVectorAlmostEqual(list(pressure.Get(1.0)), list(points[:, 0]))
        self.assertVectorAlmostEqual(list(pressure.Get(2.0)), list(10.0 * points[:, 0]))
        velocity = primvars.GetPrimvar("VELOCITY")
        self.assertAlmostEqual(velocity.Get(1.0)[1][0], 1.5, places=6)
        pair = primvars.GetPrimvar("PAIR")
        self.assertEqual(pair.GetElementSize(), 2)
        self.assertEqual(len(pair.Get(1.0)), 8)

    def test_UnchangedTopologyIsWrittenOnce(self):
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import usd_export
        from pxr import Usd, UsdGeom

        points, triangles = self._Tetrahedron()
        stage = usd_export.CreateUsdStage(self.stage_file)
        usd_export.WriteMeshTimeSample(stage, "/Kratos/Main", points, triangles, {}, 1.0)
        usd_export.WriteMeshTimeSample(stage, "/Kratos/Main", points + 0.1, triangles, {}, 2.0)
        usd_export.WriteMeshTimeSample(stage, "/Kratos/Main", points, triangles[:2], {}, 3.0)
        usd_export.SaveStage(stage)

        reopened = Usd.Stage.Open(str(self.stage_file))
        mesh = UsdGeom.Mesh(reopened.GetPrimAtPath("/Kratos/Main"))
        # sampled at 1.0 (first) and 3.0 (changed), not at 2.0 (unchanged);
        # an adaptively remeshed series stays valid, a fixed mesh compact
        self.assertEqual(mesh.GetFaceVertexIndicesAttr().GetTimeSamples(), [1.0, 3.0])
        self.assertEqual(len(mesh.GetFaceVertexIndicesAttr().Get(3.0)), 6)

    def test_PointsCloudWithWidths(self):
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import usd_export
        from pxr import Usd, UsdGeom

        points, _ = self._Tetrahedron()
        stage = usd_export.CreateUsdStage(self.stage_file)
        usd_export.WritePointsTimeSample(
            stage, "/Kratos/Cloud", points, {"PRESSURE": points[:, 2]}, 1.0,
            widths=numpy.full(4, 0.05))
        usd_export.SaveStage(stage)

        reopened = Usd.Stage.Open(str(self.stage_file))
        cloud = UsdGeom.Points(reopened.GetPrimAtPath("/Kratos/Cloud"))
        self.assertTrue(cloud)
        self.assertEqual(len(cloud.GetPointsAttr().Get(1.0)), 4)
        self.assertAlmostEqual(cloud.GetWidthsAttr().Get(1.0)[0], 0.05, places=6)

    def test_Validation(self):
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import usd_export

        points, triangles = self._Tetrahedron()
        with self.assertRaisesRegex(ValueError, "up_axis"):
            usd_export.CreateUsdStage(self.stage_file, up_axis="X")
        stage = usd_export.CreateUsdStage(self.stage_file)
        with self.assertRaisesRegex(ValueError, "triangles"):
            usd_export.WriteMeshTimeSample(
                stage, "/Kratos/Main", points, numpy.array([[0, 1, 2, 3]]), {}, 1.0)
        with self.assertRaisesRegex(ValueError, "n_points"):
            usd_export.WriteMeshTimeSample(
                stage, "/Kratos/Main", points, triangles, {"BAD": numpy.ones(3)}, 1.0)
        with self.assertRaisesRegex(ValueError, "widths"):
            usd_export.WritePointsTimeSample(
                stage, "/Kratos/Cloud", points, {}, 1.0, widths=numpy.ones(2))


@KratosUnittest.skipUnless(have_usd, "Missing required python module: pxr (usd-core).")
class TestUsdPointsKind(KratosUnittest.TestCase):
    """kind "points" gathers positions and fields directly - deliberately
    gated on usd-core alone, pinning that the meshless path needs neither
    torch nor physicsnemo."""

    def setUp(self):
        self.stage_file = Path("test_usd_cloud.usda")
        self.model = Kratos.Model()
        self.model_part = CreateCubeFixture(self.model, "Main")

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.stage_file))

    def test_PointsKindNeedsNoTessellation(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.export import usd_export_process
        from pxr import Usd, UsdGeom

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "list_of_fields"  : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
                "output_file"     : "%s",
                "kind"            : "points",
                "prim_path"       : "/Kratos/Cloud"
            }
        }""" % self.stage_file)
        process = usd_export_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteFinalize()

        reopened = Usd.Stage.Open(str(self.stage_file))
        cloud = UsdGeom.Points(reopened.GetPrimAtPath("/Kratos/Cloud"))
        self.assertTrue(cloud)
        self.assertEqual(len(cloud.GetPointsAttr().Get(1.0)), 8)
        pressure = UsdGeom.PrimvarsAPI(cloud.GetPrim()).GetPrimvar("PRESSURE")
        self.assertVectorAlmostEqual(
            list(pressure.Get(1.0)),
            [node.GetSolutionStepValue(Kratos.PRESSURE) for node in self.model_part.Nodes],
            places=5)


@KratosUnittest.skipUnless(have_usd and have_mesh_bridge,
                           "Missing required python modules: pxr (usd-core), torch, physicsnemo.")
class TestUsdExportProcess(KratosUnittest.TestCase):
    """The in-loop digital-twin exporter on a real (tessellated) model part."""

    def setUp(self):
        self.stage_file = Path("test_usd_twin.usda")
        self.model = Kratos.Model()
        self.model_part = CreateCubeFixture(self.model, "Main")

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.stage_file))

    def _Process(self, extra=""):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.export import usd_export_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "list_of_fields"  : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
                "output_file"     : "%s"%s
            }
        }""" % (self.stage_file, extra))
        return usd_export_process.Factory(settings, self.model)

    def test_TimeSampledSurfaceExport(self):
        from pxr import Usd, UsdGeom

        process = self._Process()
        process.ExecuteInitialize()
        for step in (1, 2):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            self.model_part.ProcessInfo[Kratos.TIME] = 0.5 * step
            for node in self.model_part.Nodes:
                node.SetSolutionStepValue(Kratos.PRESSURE, step * (1.0 + node.X))
            process.ExecuteFinalizeSolutionStep()
        process.ExecuteFinalize()

        reopened = Usd.Stage.Open(str(self.stage_file))
        mesh = UsdGeom.Mesh(reopened.GetPrimAtPath("/Kratos/Main"))
        self.assertTrue(mesh)
        # all 8 cube nodes survive; the tessellated hexahedron's outward
        # boundary is 12 triangles
        self.assertEqual(len(mesh.GetPointsAttr().Get(1.0)), 8)
        self.assertEqual(len(mesh.GetFaceVertexCountsAttr().Get(1.0)), 12)
        self.assertEqual(reopened.GetStartTimeCode(), 1.0)
        self.assertEqual(reopened.GetEndTimeCode(), 2.0)

        pressure = UsdGeom.PrimvarsAPI(mesh.GetPrim()).GetPrimvar("PRESSURE")
        expected = {1.0: [1.0 * (1.0 + node.X) for node in self.model_part.Nodes],
                    2.0: [2.0 * (1.0 + node.X) for node in self.model_part.Nodes]}
        for time_code, values in expected.items():
            self.assertVectorAlmostEqual(list(pressure.Get(time_code)), values, places=5)

    def test_IntervalGating(self):
        from pxr import Usd, UsdGeom

        process = self._Process(extra=""",
                "output_interval" : 2""")
        for step in (1, 2, 3):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
        process.ExecuteFinalize()

        reopened = Usd.Stage.Open(str(self.stage_file))
        mesh = UsdGeom.Mesh(reopened.GetPrimAtPath("/Kratos/Main"))
        self.assertEqual(mesh.GetPointsAttr().GetTimeSamples(), [2.0])

    def test_InvalidSettingsRaise(self):
        with self.assertRaisesRegex(ValueError, "kind"):
            self._Process(extra=""",
                "kind" : "curves\"""")
        with self.assertRaisesRegex(ValueError, "nodal"):
            from KratosMultiphysics.PhysicsNeMoApplication.processes.export import usd_export_process
            settings = Kratos.Parameters("""{
                "Parameters": {
                    "model_part_name" : "Main",
                    "list_of_fields"  : [ { "variable_name" : "PRESSURE", "data_location" : "element" } ],
                    "output_file"     : "%s"
                }
            }""" % self.stage_file)
            usd_export_process.Factory(settings, self.model)
        with self.assertRaisesRegex(ValueError, "output_interval"):
            self._Process(extra=""",
                "output_interval" : 0""")
        with self.assertRaisesRegex(ValueError, "prim_path"):
            self._Process(extra=""",
                "prim_path" : "relative/path\"""")


if __name__ == '__main__':
    KratosUnittest.main()
