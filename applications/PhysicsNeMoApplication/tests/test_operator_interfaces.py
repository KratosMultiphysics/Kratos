"""Tests for the two NON-POINTWISE point-cloud interfaces: xDeepONet, which
maps case PARAMETERS to a field at arbitrary points, and GLOBE, which maps
BOUNDARY data to the interior."""

import sys
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
    import warnings as _warnings
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")
        from physicsnemo.experimental.models.xdeeponet import DeepONet
    have_deeponet = True
except ImportError:
    have_deeponet = False

try:
    import warnings as _warnings
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")
        from physicsnemo.experimental.models.globe import GLOBE
    have_globe = True
except ImportError:
    have_globe = False

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))


def _TinyDeepONet(branch_width=2, out_channels=1, dimension=3, width=16, seed=0):
    torch.manual_seed(seed)
    branch = torch.nn.Sequential(
        torch.nn.Linear(branch_width, 32), torch.nn.Tanh(), torch.nn.Linear(32, width))
    trunk = torch.nn.Sequential(
        torch.nn.Linear(dimension, 32), torch.nn.Tanh(), torch.nn.Linear(32, width))
    return DeepONet(branch, trunk=trunk, dimension=dimension, width=width,
                    out_channels=out_channels)


@KratosUnittest.skipUnless(have_torch and have_deeponet,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestDeepONetContract(KratosUnittest.TestCase):
    def test_CoreModeShapes(self):
        """The contract the bridge is written against: an MLP branch takes
        (B, D) case parameters, the trunk (T, d) query coordinates, and the
        output is (B, T, C_out) - one field per case, at the queried
        points."""
        model = _TinyDeepONet().eval()
        with torch.no_grad():
            single = model(torch.randn(1, 2), torch.rand(37, 3))
            batched = model(torch.randn(5, 2), torch.rand(37, 3))
        self.assertEqual(list(single.shape), [1, 37, 1])
        self.assertEqual(list(batched.shape), [5, 37, 1])

    def test_ItCannotBeCheckpointedAsMdlus(self):
        """Worth pinning, because the error arrives at SAVE time after a
        model has trained: physicsnemo's Module.save refuses plain torch
        submodules, and wrapping them with Module.from_torch instead makes
        the constructor metadata unserializable. Trace it."""
        model = _TinyDeepONet()
        with self.assertRaises(TypeError):
            model.save("unusable.mdlus")
        with self.assertRaises(Exception):
            torch.jit.script(model)

        traced = torch.jit.trace(model.eval(), (torch.randn(1, 2), torch.rand(7, 3)))
        with torch.no_grad():  # the trace holds at a different point count
            self.assertEqual(
                list(traced(torch.randn(1, 2), torch.rand(31, 3)).shape), [1, 31, 1])

    def test_TheSameCaseAtDifferentPointsIsOneField(self):
        """The trunk is queried pointwise, so a subset of the points gives
        the same values as the full set - that is what makes it an
        operator rather than a grid model."""
        model = _TinyDeepONet().eval()
        parameters = torch.randn(1, 2)
        points = torch.rand(20, 3)
        with torch.no_grad():
            everything = model(parameters, points)
            subset = model(parameters, points[:5])
        torch.testing.assert_close(subset[0], everything[0][:5])


@KratosUnittest.skipUnless(have_torch and have_deeponet,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestDeepONetThroughTheProcess(KratosUnittest.TestCase):
    @staticmethod
    def _SaveTraced(model, path, branch_width):
        """An xDeepONet cannot be checkpointed as .mdlus (Module.save
        refuses plain torch submodules, and physicsnemo-wrapped ones make
        its metadata unserializable) and cannot be scripted (its forward
        takes *args). Tracing works, and the trace stays valid at any
        number of query points - which is what a changing mesh needs."""
        traced = torch.jit.trace(
            model.eval(), (torch.randn(1, branch_width), torch.rand(7, 3)))
        torch.jit.save(traced, str(path))

    def setUp(self):
        self.checkpoint = Path("test_deeponet.pt")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        self.model_part.ProcessInfo[Kratos.TIME] = 0.25

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting(str(self.checkpoint))
        kratos_utils.DeleteFileIfExisting(str(self.checkpoint) + ".card.json")

    def _Run(self, branch_block, trunk_dimension=0):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import (
            point_cloud_inference_process)

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"       : "Main",
                "model_settings"        : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "torchscript",
                    "device"          : "cpu"
                },
                "model_interface"       : "deeponet",
                "normalize_coordinates" : false,
                "trunk_dimension"       : %d,
                "branch_input"          : %s,
                "input_fields"          : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"         : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""" % (self.checkpoint, trunk_dimension, branch_block))
        process = point_cloud_inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        return numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])

    def test_CaseParametersFromTheProcessInfoReachTheBranch(self):
        self._SaveTraced(_TinyDeepONet(branch_width=2), self.checkpoint, 2)
        written = self._Run(
            '{ "process_info_variables" : ["TIME"], "constants" : [2.5] }')
        self.assertEqual(len(written), self.model_part.NumberOfNodes())
        self.assertTrue(numpy.isfinite(written).all())
        # an operator's field is not constant over the nodes
        self.assertGreater(written.std(), 0.0)

    def test_ADifferentCaseGivesADifferentField(self):
        """The branch is what makes it an operator: change the case
        parameters and the whole field changes, with the mesh untouched."""
        self._SaveTraced(_TinyDeepONet(branch_width=1), self.checkpoint, 1)
        first = self._Run('{ "constants" : [0.0] }')
        second = self._Run('{ "constants" : [3.0] }')
        self.assertGreater(numpy.abs(first - second).max(), 0.0)

    def test_ProcessInfoValueIsReadAtRunTime(self):
        self._SaveTraced(_TinyDeepONet(branch_width=1), self.checkpoint, 1)
        self.model_part.ProcessInfo[Kratos.TIME] = 0.1
        early = self._Run('{ "process_info_variables" : ["TIME"] }')
        self.model_part.ProcessInfo[Kratos.TIME] = 9.0
        late = self._Run('{ "process_info_variables" : ["TIME"] }')
        self.assertGreater(numpy.abs(early - late).max(), 0.0)

    def test_AnEmptyBranchInputRaises(self):
        self._SaveTraced(_TinyDeepONet(branch_width=1), self.checkpoint, 1)
        with self.assertRaisesRegex(ValueError, "branch_input"):
            self._Run('{ }')

    def test_TrunkDimensionCutsTheCoordinates(self):
        """A 2-D case queries a 2-D trunk; the process follows DOMAIN_SIZE
        unless told otherwise."""
        model = _TinyDeepONet(branch_width=1, dimension=2).eval()
        torch.jit.save(
            torch.jit.trace(model, (torch.randn(1, 1), torch.rand(7, 2))),
            str(self.checkpoint))
        written = self._Run('{ "constants" : [1.0] }', trunk_dimension=2)
        self.assertTrue(numpy.isfinite(written).all())

        self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
        followed = self._Run('{ "constants" : [1.0] }', trunk_dimension=0)
        numpy.testing.assert_allclose(followed, written, atol=1e-12)


@KratosUnittest.skipUnless(have_torch and have_globe,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestGlobeBridge(KratosUnittest.TestCase):
    """GLOBE on CPU: the cluster tree, the forward and the backward all run
    there, so this is not a GPU-only capability."""

    @staticmethod
    def _TinyGlobe(seed=0):
        import warnings

        torch.manual_seed(seed)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from physicsnemo.experimental.models.globe import GLOBE
        return GLOBE(
            n_spatial_dims=2,
            output_field_ranks={"u": 0},
            boundary_source_data_ranks={"wall": {"TEMPERATURE": 0}},
            reference_length_names=["L"], reference_area=1.0,
            n_communication_hyperlayers=1, n_latent_scalars=2, n_latent_vectors=1,
            hidden_layer_sizes=[8], tree_build_device="cpu",
            use_gradient_checkpointing=False)

    @staticmethod
    def _Boundary(seed=0):
        from physicsnemo.mesh import Mesh

        rng = numpy.random.default_rng(seed)
        points = torch.tensor([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
                              dtype=torch.float32)
        cells = torch.tensor([[0, 1], [1, 2], [2, 3], [3, 0]])
        mesh = Mesh(points=points, cells=cells)
        mesh.cell_data["TEMPERATURE"] = torch.tensor(
            rng.uniform(0.0, 1.0, size=4), dtype=torch.float32)
        return {"wall": mesh}

    def test_ForwardAndBackwardRunOnCpu(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import globe_bridge

        model = self._TinyGlobe()
        points = torch.rand(12, 2)
        prediction = globe_bridge.RunGlobeForward(
            model, points, self._Boundary(), {"L": 1.0}, ["u"], enable_grad=True)
        self.assertEqual(list(prediction.shape), [12, 1])
        self.assertTrue(bool(torch.isfinite(prediction).all()))

        prediction.sum().backward()
        self.assertTrue(any(p.grad is not None for p in model.parameters()))

    def test_UnknownOutputNameIsReported(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import globe_bridge

        with self.assertRaisesRegex(ValueError, "output field"):
            globe_bridge.RunGlobeForward(
                self._TinyGlobe(), torch.rand(5, 2), self._Boundary(),
                {"L": 1.0}, ["pressure"])

    def test_RankZeroCellDataMustNotBeAColumn(self):
        """The trap: upstream validates the declared rank against the
        array's own rank, so a rank-0 field shaped (n_cells, 1) is
        rejected while (n_cells,) is accepted."""
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import globe_bridge

        boundary = self._Boundary()
        boundary["wall"].cell_data["TEMPERATURE"] = \
            boundary["wall"].cell_data["TEMPERATURE"].reshape(-1, 1)
        with self.assertRaises(ValueError):
            globe_bridge.RunGlobeForward(
                self._TinyGlobe(), torch.rand(5, 2), boundary, {"L": 1.0}, ["u"])

    def test_TrainingLowersTheLoss(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import globe_training

        model = self._TinyGlobe()
        points = torch.rand(16, 2)
        cases = []
        for seed in range(3):
            boundary = self._Boundary(seed)
            # a target that genuinely depends on the boundary data, so the
            # model has something a boundary-driven kernel can learn
            level = float(boundary["wall"].cell_data["TEMPERATURE"].mean())
            cases.append((boundary, {"L": 1.0}, points,
                          numpy.full((16, 1), level, dtype=numpy.float64)))

        history = globe_training.TrainGlobe(model, cases, Kratos.Parameters("""{
            "epochs"        : 8,
            "learning_rate" : 5e-2,
            "device"        : "cpu",
            "seed"          : 0
        }"""), output_names=["u"])
        self.assertEqual(len(history), 8)
        self.assertTrue(numpy.isfinite(history).all())
        self.assertLess(min(history[1:]), history[0])

    def test_TrainGlobeValidatesItsSettings(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import globe_training

        with self.assertRaisesRegex(ValueError, "at least one case"):
            globe_training.TrainGlobe(self._TinyGlobe(), [], Kratos.Parameters("{}"),
                                      output_names=["u"])


@KratosUnittest.skipUnless(have_torch and have_globe,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestGlobeOnAKratosModelPart(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.TEMPERATURE, Kratos.PRESSURE))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.TEMPERATURE, node.X)
        # a named boundary: the x = 0 face, with real conditions on it -
        # a sub-model-part carrying only NODES has no geometry to tessellate
        boundary = self.model_part.CreateSubModelPart("Inlet")
        face_nodes = sorted(node.Id for node in self.model_part.Nodes if node.X < 1e-12)
        boundary.AddNodes(face_nodes)
        properties = self.model_part.GetProperties()[1]
        condition_id = 1
        for i in range(0, len(face_nodes) - 2, 1):
            boundary.CreateNewCondition(
                "SurfaceCondition3D3N", condition_id,
                [face_nodes[i], face_nodes[i + 1], face_nodes[i + 2]], properties)
            condition_id += 1

    def test_BoundaryMeshesCarryTheirValuesInCellData(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import globe_bridge

        meshes = globe_bridge.BuildGlobeBoundaryMeshes(
            self.model_part, ["Inlet"], [("TEMPERATURE", "node_historical")])
        self.assertIn("Inlet", meshes)
        mesh = meshes["Inlet"]
        self.assertIn("TEMPERATURE", mesh.cell_data)
        # rank 0 means one value per cell, NOT a column
        self.assertEqual(tuple(mesh.cell_data["TEMPERATURE"].shape),
                         (mesh.cells.shape[0],))
        # the inlet sits at x = 0, where the field is 0
        self.assertLess(float(mesh.cell_data["TEMPERATURE"].abs().max()), 1e-6)

    def test_AKratosBuiltBoundaryForwardsThroughGlobe(self):
        """The whole path: Kratos model part -> boundary mesh -> GLOBE.

        Everything else here checks the two halves separately - the bridge
        builds a boundary, and the forward tests run GLOBE on a hand-built
        float32 mesh. Nothing joined them, which is exactly how a dtype
        mismatch survived: the mesh bridge carries a boundary's POINTS as
        Kratos doubles while its cell data is already float32, and GLOBE's
        kernels are float32, so a Kratos-built boundary used to fail with
        "mat1 and mat2 must have the same dtype" from deep inside a Linear,
        naming nothing that would lead you back to the mesh.
        """
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import globe_bridge

        meshes = globe_bridge.BuildGlobeBoundaryMeshes(
            self.model_part, ["Inlet"], [("TEMPERATURE", "node_historical")])
        self.assertEqual(meshes["Inlet"].points.dtype, torch.float64,
                         "the bridge is expected to hand back Kratos doubles here; "
                         "if that changes, this regression no longer covers anything")

        torch.manual_seed(0)
        model = GLOBE(
            n_spatial_dims=3, output_field_ranks={"u": 0},
            boundary_source_data_ranks={"Inlet": {"TEMPERATURE": 0}},
            reference_length_names=["L"], reference_area=1.0,
            n_communication_hyperlayers=1, n_latent_scalars=2, n_latent_vectors=1,
            hidden_layer_sizes=[8], tree_build_device="cpu",
            use_gradient_checkpointing=False)

        points = numpy.array([[node.X, node.Y, node.Z]
                              for node in self.model_part.Nodes], dtype=float)
        prediction = globe_bridge.RunGlobeForward(
            model, points, meshes, {"L": 1.0}, ["u"])

        self.assertEqual(tuple(prediction.shape), (points.shape[0], 1))
        self.assertTrue(bool(torch.isfinite(prediction).all()))
        # the caller's mesh must survive untouched: the inference process
        # caches these per node count and would otherwise be handed a mesh
        # silently converted underneath it
        self.assertEqual(meshes["Inlet"].points.dtype, torch.float64)

    def test_AnUnknownBoundaryIsReported(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import globe_bridge

        with self.assertRaisesRegex(ValueError, "no sub-model-part"):
            globe_bridge.BuildGlobeBoundaryMeshes(
                self.model_part, ["Outlet"], [("TEMPERATURE", "node_historical")])

    def test_ANodesOnlyBoundaryIsReported(self):
        """A sub-model-part carrying only nodes has no geometry to
        tessellate, so it yields no boundary mesh at all."""
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import globe_bridge

        self.model_part.CreateSubModelPart("NodesOnly").AddNodes(
            [node.Id for node in self.model_part.Nodes if node.Y < 1e-12])
        with self.assertRaisesRegex(ValueError, "conditions"):
            globe_bridge.BuildGlobeBoundaryMeshes(
                self.model_part, ["NodesOnly"], [("TEMPERATURE", "node_historical")])


if __name__ == '__main__':
    KratosUnittest.main()
