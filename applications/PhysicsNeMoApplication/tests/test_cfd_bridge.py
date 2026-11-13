import math

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges import cfd_bridge
from KratosMultiphysics.PhysicsNeMoApplication.processes import validation_metrics_process
try:
    import pyvista  # noqa: F401
    have_pyvista = True
except ImportError:
    have_pyvista = False

try:
    import physicsnemo.cfd.evaluation.metrics  # noqa: F401
    import physicsnemo.cfd.hybrid_initialization_tools  # noqa: F401
    have_physicsnemo_cfd = True
except ImportError:
    have_physicsnemo_cfd = False

try:
    import torch  # noqa: F401
    have_torch = True
except ImportError:
    have_torch = False


def CreateCubeFixture(model, name, pressure_fn, velocity):
    """Unit cube: 1 hexahedron + 6 quad skin conditions, nodal p and U."""
    model_part = model.CreateModelPart(name)
    model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
    model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
    model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
    model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
    props = model_part.CreateNewProperties(1)
    corners = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
               (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]
    for i, xyz in enumerate(corners):
        node = model_part.CreateNewNode(i + 1, *xyz)
        node.SetSolutionStepValue(Kratos.PRESSURE, pressure_fn(node))
        node.SetSolutionStepValue(Kratos.VELOCITY, list(velocity))
    model_part.CreateNewElement("Element3D8N", 1, list(range(1, 9)), props)
    faces = [(1, 4, 3, 2), (5, 6, 7, 8), (1, 2, 6, 5), (2, 3, 7, 6), (3, 4, 8, 7), (4, 1, 5, 8)]
    for i, face in enumerate(faces):
        model_part.CreateNewCondition("SurfaceCondition3D4N", i + 1, list(face), props)
    return model_part


@KratosUnittest.skipUnless(have_pyvista, "Missing required python module: pyvista.")
class TestCfdBridgePolyData(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = CreateCubeFixture(
            self.model, "Main", lambda node: 1.0 + node.X + 2.0 * node.Y, [1.0, 0.0, 0.0])

    def test_SurfacePolyData(self):
        polydata, provenance = cfd_bridge.ModelPartToPolyData(
            self.model_part,
            [("p", "PRESSURE", "node_historical"), ("U", "VELOCITY", "node_historical")])
        self.assertEqual(polydata.n_points, 8)
        self.assertEqual(polydata.n_cells, 12)
        points = numpy.asarray(polydata.points)
        expected_p = 1.0 + points[:, 0] + 2.0 * points[:, 1]
        self.assertVectorAlmostEqual(
            numpy.asarray(polydata.point_data["p"]).ravel(), expected_p, 12)
        self.assertMatrixAlmostEqual(
            Kratos.Matrix(numpy.asarray(polydata.point_data["U"])),
            Kratos.Matrix(numpy.tile([1.0, 0.0, 0.0], (8, 1))), 12)
        self.assertEqual(provenance.number_of_points, 8)

    def test_VolumeContainerRaises(self):
        with self.assertRaisesRegex(RuntimeError, "surface tessellation"):
            cfd_bridge.ModelPartToPolyData(self.model_part, (), "Elements")

    def test_NonNodalLocationRaises(self):
        with self.assertRaisesRegex(ValueError, "must be nodal"):
            cfd_bridge.ModelPartToPolyData(
                self.model_part, [("p", "PRESSURE", "condition")])

    def test_NodesToPolyData(self):
        polydata, node_ids = cfd_bridge.NodesToPolyData(
            self.model_part, [("p", "PRESSURE", "node_historical")])
        self.assertEqual(polydata.n_points, 8)
        self.assertEqual(node_ids, list(range(1, 9)))
        points = numpy.asarray(polydata.points)
        expected_p = 1.0 + points[:, 0] + 2.0 * points[:, 1]
        self.assertVectorAlmostEqual(
            numpy.asarray(polydata.point_data["p"]).ravel(), expected_p, 12)

    def test_CurvedSurfacePolyData(self):
        # A curved Quadrilateral3D9 skin: the PolyData carries the synthetic
        # lattice points ((2^k+1)^2 = 25 at k=2) and the gather interpolates
        # the (linear-in-x) pressure exactly at every one of them.
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

        polydata, provenance = cfd_bridge.ModelPartToPolyData(
            model_part, [("p", "PRESSURE", "node_historical")],
            higher_order_mode="curved", curved_refinement_levels=2)
        self.assertEqual(polydata.n_points, 25)
        self.assertEqual(polydata.n_cells, 32)
        self.assertEqual(provenance.number_of_synthetic_points, 25 - 9)
        points = numpy.asarray(polydata.points)
        self.assertVectorAlmostEqual(
            numpy.asarray(polydata.point_data["p"]).ravel(), 1.0 + points[:, 0], 12)


@KratosUnittest.skipUnless(have_physicsnemo_cfd, "Missing required python module: physicsnemo.cfd.")
class TestCfdMetrics(KratosUnittest.TestCase):
    def test_L2PressureHandComputed(self):
        ground_truth = {"pressure": numpy.array([1.0, 2.0, 3.0])}
        predictions = {"pressure": numpy.array([1.1, 1.9, 3.2])}
        values = cfd_bridge.EvaluateCfdMetrics(
            [("l2_pressure", "surface")], ground_truth, predictions)
        expected = (numpy.linalg.norm(predictions["pressure"] - ground_truth["pressure"])
                    / numpy.linalg.norm(ground_truth["pressure"]))
        self.assertAlmostEqual(values["l2_pressure"], expected, places=12)

    def test_UnknownMetricRaises(self):
        with self.assertRaisesRegex(KeyError, "Unknown metric"):
            cfd_bridge.EvaluateCfdMetrics([("no_such_metric", "surface")], {}, {})

    def test_ListCfdMetrics(self):
        names = cfd_bridge.ListCfdMetrics()
        self.assertIn("l2_pressure", names)
        self.assertIn("drag", names)


@KratosUnittest.skipUnless(have_physicsnemo_cfd and have_pyvista,
                           "Missing required python modules: physicsnemo.cfd / pyvista.")
class TestHybridInitializationDelegation(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.part_a = CreateCubeFixture(
            self.model, "A", lambda node: 1.0 + node.X, [1.0, 0.0, 0.0])
        self.part_b = CreateCubeFixture(
            self.model, "B", lambda node: 10.0, [0.0, 2.0, 0.0])
        flowfield_settings = """{
            "velocity_variable" : "VELOCITY",
            "pressure_variable" : "PRESSURE"
        }"""
        self.flowfield_a, self.provenance_a = cfd_bridge.CreateFlowfield(
            self.part_a, Kratos.Parameters(flowfield_settings))
        self.flowfield_b, _ = cfd_bridge.CreateFlowfield(
            self.part_b, Kratos.Parameters(flowfield_settings))

    def test_FlowfieldFillsMissingFields(self):
        for name in ("U", "p", "k", "omega"):
            self.assertIn(name, self.flowfield_a.mesh.point_data)
        self.assertVectorAlmostEqual(
            numpy.asarray(self.flowfield_a.mesh.point_data["k"]), numpy.zeros(8), 12)

    def test_ConstantBlend(self):
        blended = cfd_bridge.CreateHybridInitialization(
            self.flowfield_a, self.flowfield_b,
            Kratos.Parameters("""{ "constant_weight" : 0.25 }"""))
        points = numpy.asarray(blended.mesh.points)
        expected_p = 0.25 * (1.0 + points[:, 0]) + 0.75 * 10.0
        self.assertVectorAlmostEqual(
            numpy.asarray(blended.mesh.point_data["p"]).ravel(), expected_p, 12)
        expected_U = numpy.tile([0.25, 1.5, 0.0], (8, 1))
        self.assertMatrixAlmostEqual(
            Kratos.Matrix(numpy.asarray(blended.mesh.point_data["U"])),
            Kratos.Matrix(expected_U), 12)

    def test_WeightOneReturnsFlowfieldA(self):
        # Pins the shared-mesh path computed locally (upstream 0.0.3a0's
        # no-interpolation branch would return flowfield_b regardless).
        blended = cfd_bridge.CreateHybridInitialization(
            self.flowfield_a, self.flowfield_b,
            Kratos.Parameters("""{ "constant_weight" : 1.0 }"""))
        points = numpy.asarray(blended.mesh.points)
        self.assertVectorAlmostEqual(
            numpy.asarray(blended.mesh.point_data["p"]).ravel(), 1.0 + points[:, 0], 12)

    def test_CallableBlendOverride(self):
        weight_fn = lambda a, b: numpy.asarray(a.mesh.points)[:, 0]  # noqa: E731
        blended = cfd_bridge.CreateHybridInitialization(
            self.flowfield_a, self.flowfield_b, blend_strategy=weight_fn)
        points = numpy.asarray(blended.mesh.points)
        x = points[:, 0]
        expected_p = x * (1.0 + x) + (1.0 - x) * 10.0
        self.assertVectorAlmostEqual(
            numpy.asarray(blended.mesh.point_data["p"]).ravel(), expected_p, 12)

    def test_InvalidBlendSettingsRaise(self):
        with self.assertRaisesRegex(ValueError, "constant_weight"):
            cfd_bridge.CreateHybridInitialization(
                self.flowfield_a, self.flowfield_b,
                Kratos.Parameters("""{ "constant_weight" : 1.5 }"""))
        with self.assertRaisesRegex(ValueError, "blend_strategy"):
            cfd_bridge.CreateHybridInitialization(
                self.flowfield_a, self.flowfield_b,
                Kratos.Parameters("""{ "blend_strategy" : "no_such_strategy" }"""))

    def test_WriteBack(self):
        blended = cfd_bridge.CreateHybridInitialization(
            self.flowfield_a, self.flowfield_b,
            Kratos.Parameters("""{ "constant_weight" : 0.5 }"""))
        cfd_bridge.FlowfieldToModelPart(
            blended, self.part_a, self.provenance_a,
            {"p": "TEMPERATURE", "U": "DISPLACEMENT"})
        for node in self.part_a.Nodes:
            expected_p = 0.5 * (1.0 + node.X) + 0.5 * 10.0
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE), expected_p, places=12)
            self.assertVectorAlmostEqual(
                node.GetSolutionStepValue(Kratos.DISPLACEMENT), [0.5, 1.0, 0.0], 12)


@KratosUnittest.skipUnless(have_physicsnemo_cfd and have_torch,
                           "Missing required python modules: physicsnemo.cfd / torch.")
class TestValidationProcessCfdMetrics(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
        # predicted = reference + 1 -> relative L2 = 2 / sqrt(14)
        for i in range(4):
            node = self.model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.TEMPERATURE, float(i))       # reference
            node.SetSolutionStepValue(Kratos.PRESSURE, float(i) + 1.0)    # predicted

    def test_CfdMetricRecord(self):
        process = validation_metrics_process.Factory(Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"     : "Main",
                "list_of_comparisons" : [],
                "cfd_metrics"         : [
                    {
                        "name"   : "l2_pressure",
                        "domain" : "surface",
                        "fields" : {
                            "pressure" : {
                                "predicted_variable" : "PRESSURE",
                                "reference_variable" : "TEMPERATURE"
                            }
                        }
                    }
                ],
                "output_file"         : "unused_cfd_metrics.json"
            }
        }"""), self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        self.assertEqual(len(process.history), 1)
        self.assertAlmostEqual(
            process.history[0]["cfd_l2_pressure"], 2.0 / math.sqrt(14.0), places=12)

    def test_MissingFieldsRaises(self):
        with self.assertRaisesRegex(ValueError, "fields"):
            validation_metrics_process.Factory(Kratos.Parameters("""{
                "Parameters": {
                    "model_part_name"     : "Main",
                    "list_of_comparisons" : [],
                    "cfd_metrics"         : [ { "name" : "l2_pressure" } ]
                }
            }"""), self.model)


@KratosUnittest.skipUnless(have_physicsnemo_cfd,
                           "Missing required python module: physicsnemo-cfd.")
class TestCfdEvaluationWrappers(KratosUnittest.TestCase):
    def test_GeoTransolverWrapperModulesResolve(self):
        # the checkpoint/NGC-config driven evaluation wrappers must resolve;
        # feeding them Kratos data goes through ModelPartToPolyData (docs recipe)
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # they pull in physicsnemo.experimental
            for name in ("geotransolver", "geotransolver_gp", "geotransolver_drivaerstar"):
                module = cfd_bridge._TryImportCfdEvaluationWrappers(name)
                self.assertIsNotNone(module)

    def test_UnknownWrapperRaisesActionable(self):
        with self.assertRaisesRegex(ImportError, "physicsnemo-cfd"):
            cfd_bridge._TryImportCfdEvaluationWrappers("no_such_wrapper_module")


if __name__ == '__main__':
    KratosUnittest.main()
