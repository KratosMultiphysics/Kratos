import copy
from pathlib import Path

import os
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from test_cae_dataset_export_process import CreateCaeFixture

try:
    import torch  # noqa: F401
    import warp  # noqa: F401 - DoMINO preprocessing needs it (SDF)
    import physicsnemo.models.domino  # noqa: F401
    import physicsnemo.datapipes.cae.cae_dataset  # noqa: F401
    have_domino = True
except ImportError:
    have_domino = False


def _ProcessSettings(scratch, checkpoint="", extra=""):
    return Kratos.Parameters("""{
        "Parameters": {
            "volume_model_part_name"  : "Main",
            "surface_model_part_name" : "Main.Skin",
            "model_settings"          : {
                "checkpoint_file" : "%s",
                "checkpoint_type" : "physicsnemo",
                "device"          : "cpu"
            },
            "model_type"           : "surface",
            "bounding_box"         : [-0.1, -0.1, -0.1, 1.1, 1.1, 1.1],
            "global_params"        : { "stream_velocity" : 30.0, "air_density" : 1.226 },
            "global_params_order"  : ["stream_velocity", "air_density"],
            "datapipe_overrides"   : { "grid_resolution" : [8, 8, 8], "num_surface_neighbors" : 4 },
            "output_fields_surface": [ { "variable_name" : "TEMPERATURE", "data_location" : "condition" } ],
            "scratch_directory"    : "%s"
            %s
        }
    }""" % (checkpoint, scratch, extra))


def _FindDominoCheckpointDir():
    """The public domino_drivaerml surface checkpoint, if it is already local.

    Nothing is downloaded during a test run: set
    PHYSICSNEMO_DOMINO_CHECKPOINT_DIR, or have it in the default HuggingFace
    cache, otherwise the checkpoint-backed tests self-skip.
    """
    override = os.environ.get("PHYSICSNEMO_DOMINO_CHECKPOINT_DIR", "")
    candidates = [Path(override)] if override else []
    for root in (Path(os.environ.get("HF_HOME", Path.home() / ".cache" / "huggingface")),
                 Path("/tmp/hf_domino/cache")):
        candidates.extend(root.glob(
            "**/models--nvidia--domino_drivaerml/snapshots/*/domino_drivaerml_surface_checkpoint"))
    for candidate in candidates:
        if (candidate / "DoMINO.0.501.mdlus").is_file() and \
                (candidate / "scaling_factors.pkl").is_file():
            return candidate
    return None


_CHECKPOINT_DIR = _FindDominoCheckpointDir()


class _FakeScalingFactors:
    """Stands in for physicsnemo-cfd's ScalingFactors so the de-normalization
    maths can be tested without a 48 MB checkpoint - or physicsnemo at all."""

    def __init__(self):
        self.max_val = {"surface_fields": numpy.array([2.0, 10.0])}
        self.min_val = {"surface_fields": numpy.array([-2.0, 0.0])}
        self.mean = {"surface_fields": numpy.array([1.0, 5.0])}
        self.std = {"surface_fields": numpy.array([2.0, 4.0])}


def _MakeBareProcess(**attributes):
    """A DominoInferenceProcess with only the de-normalization attributes set,
    bypassing __init__ (which needs a Model and a mesh)."""
    from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import domino_inference_process
    process = domino_inference_process.DominoInferenceProcess.__new__(
        domino_inference_process.DominoInferenceProcess)
    process._scaling = _FakeScalingFactors()
    process.scaling_factors_file = "unused"
    process.normalization = "none"
    process.redimensionalize = False
    process._global_params_values = {}
    for name, value in attributes.items():
        setattr(process, name, value)
    return process


class TestDominoDenormalization(KratosUnittest.TestCase):
    """The de-normalization the process previously did not do.

    A pretrained DoMINO emits dimensionless, min-max normalized values;
    writing them straight onto Kratos entities is wrong by roughly three
    orders of magnitude. Needs no torch and no checkpoint, so it runs
    everywhere - including the torch-free CI.
    """

    def test_NoneIsTheIdentity(self):
        # existing configurations were written against raw output and must
        # keep getting exactly that
        process = _MakeBareProcess(normalization="none", redimensionalize=False)
        raw = numpy.array([[0.5, -0.25]])
        self.assertTrue(numpy.array_equal(process._Denormalize(raw, "surface_fields"), raw))

    def test_MinMaxScalingMatchesTheClosedForm(self):
        # x * (max - min) / 2 + (max + min) / 2
        process = _MakeBareProcess(normalization="min_max_scaling")
        raw = numpy.array([[0.0, 1.0], [-1.0, 0.5]])
        expected = raw * numpy.array([2.0, 5.0]) + numpy.array([0.0, 5.0])
        numpy.testing.assert_allclose(
            process._Denormalize(raw, "surface_fields"), expected, rtol=1e-12)
        # the endpoints map onto the recorded range, which is what makes it
        # an inverse rather than an arbitrary affine map
        endpoints = process._Denormalize(numpy.array([[-1.0, 1.0]]), "surface_fields")
        numpy.testing.assert_allclose(endpoints[0], [-2.0, 10.0], rtol=1e-12)

    def test_MeanStdScaling(self):
        process = _MakeBareProcess(normalization="mean_std_scaling")
        raw = numpy.array([[0.0, 2.0]])
        numpy.testing.assert_allclose(
            process._Denormalize(raw, "surface_fields"), [[1.0, 13.0]], rtol=1e-12)

    def test_RedimensionalizationAppliesDynamicPressure(self):
        process = _MakeBareProcess(
            normalization="none", redimensionalize=True,
            _global_params_values={"stream_velocity": [30.0], "air_density": 1.205})
        self.assertAlmostEqual(process._DynamicPressure(), 30.0 ** 2 * 1.205, places=9)
        raw = numpy.array([[1.0, 2.0]])
        numpy.testing.assert_allclose(
            process._Denormalize(raw, "surface_fields"), raw * 1084.5, rtol=1e-12)

    def test_MissingGlobalParamsIsActionable(self):
        process = _MakeBareProcess(redimensionalize=True,
                                   _global_params_values={"air_density": 1.2})
        with self.assertRaisesRegex(ValueError, "velocity"):
            process._DynamicPressure()

    def test_ChannelCountMismatchIsCaught(self):
        # the commonest way to pair the wrong scaling_factors.pkl with a model
        process = _MakeBareProcess(normalization="min_max_scaling")
        with self.assertRaisesRegex(ValueError, "does not belong to this checkpoint"):
            process._Denormalize(numpy.zeros((3, 5)), "surface_fields")


class TestDominoSurfaceScatter(KratosUnittest.TestCase):
    """The provenance-based triangle -> parent-entity collapse (numpy only,
    no ML dependencies)."""

    def setUp(self):
        self.scratch = Path("test_domino_scatter_scratch")
        self.model = Kratos.Model()
        self.model_part = CreateCaeFixture(self.model)

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.scratch))

    def test_TriangleValuesCollapseToParentConditions(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import domino_inference_process
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

        process = domino_inference_process.Factory(
            _ProcessSettings(self.scratch), self.model)

        # synthetic per-triangle prediction = the parent condition's Id: the
        # collapsed entity value must be exactly that Id (order-exact check)
        provenance = domain_mesh_builder.BuildProvenance(
            self.model_part.GetSubModelPart("Skin"), "Conditions",
            "smallest_id_diagonal", "reduce", 2)
        per_triangle = provenance.cell_provenance[:, 0].astype(float)[:, None]
        process._WriteSurfaceOutputs(per_triangle)

        for condition in self.model_part.GetSubModelPart("Skin").Conditions:
            self.assertAlmostEqual(
                condition.GetValue(Kratos.TEMPERATURE), float(condition.Id), places=12)

    def test_TheProvenanceMapIsBuiltOnceAcrossSteps(self):
        """The map is topological, so a second step must reuse it.

        It is the most expensive per-entity path in benchmark_bridges.py, and
        this process used to rebuild it on every ExecuteFinalizeSolutionStep.
        """
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import domino_inference_process
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

        process = domino_inference_process.Factory(
            _ProcessSettings(self.scratch), self.model)
        skin = self.model_part.GetSubModelPart("Skin")
        provenance = domain_mesh_builder.BuildProvenance(
            skin, "Conditions", "smallest_id_diagonal", "reduce", 2)
        per_triangle = provenance.cell_provenance[:, 0].astype(float)[:, None]

        calls = []
        original = domino_inference_process.domain_mesh_builder.BuildProvenance

        def Counting(*args, **kwargs):
            calls.append(1)
            return original(*args, **kwargs)

        domino_inference_process.domain_mesh_builder.BuildProvenance = Counting
        try:
            process._WriteSurfaceOutputs(per_triangle)
            first = [c.GetValue(Kratos.TEMPERATURE) for c in skin.Conditions]
            self.assertEqual(len(calls), 1)

            for _ in range(3):
                process._WriteSurfaceOutputs(per_triangle)
            # still one build, and the answer is unchanged by caching
            self.assertEqual(len(calls), 1)
            second = [c.GetValue(Kratos.TEMPERATURE) for c in skin.Conditions]
        finally:
            domino_inference_process.domain_mesh_builder.BuildProvenance = original

        numpy.testing.assert_allclose(second, first, rtol=0.0, atol=0.0)

    def test_ChangedTopologyRebuildsTheProvenanceMap(self):
        """AdaptiveRemeshProcess ships here, so the mesh can change mid-run."""
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import domino_inference_process
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

        process = domino_inference_process.Factory(
            _ProcessSettings(self.scratch), self.model)
        skin = self.model_part.GetSubModelPart("Skin")

        def PerTriangle():
            provenance = domain_mesh_builder.BuildProvenance(
                skin, "Conditions", "smallest_id_diagonal", "reduce", 2)
            return provenance.cell_provenance[:, 0].astype(float)[:, None]

        process._WriteSurfaceOutputs(PerTriangle())
        before = len(skin.Conditions)

        # grow the skin: a stale map would either raise or write the old layout
        properties = self.model_part.GetProperties()[1]
        next_node = max(node.Id for node in self.model_part.Nodes) + 1
        skin.CreateNewNode(next_node, 0.5, 0.5, 2.0)
        node_ids = [node.Id for node in skin.Nodes][:2] + [next_node]
        next_condition = max(c.Id for c in self.model_part.Conditions) + 1
        skin.CreateNewCondition(
            "SurfaceCondition3D3N", next_condition, node_ids, properties)

        self.assertEqual(len(skin.Conditions), before + 1)
        process._WriteSurfaceOutputs(PerTriangle())
        # every condition, including the new one, now carries its own Id
        for condition in skin.Conditions:
            self.assertAlmostEqual(
                condition.GetValue(Kratos.TEMPERATURE), float(condition.Id), places=12)

    def test_ValidationErrors(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import domino_inference_process
        settings = _ProcessSettings(self.scratch)
        settings["Parameters"]["bounding_box"].SetVector(Kratos.Vector())
        with self.assertRaisesRegex(ValueError, "bounding_box"):
            domino_inference_process.Factory(settings, self.model)

        settings = _ProcessSettings(self.scratch)
        settings["Parameters"]["output_fields_surface"][0]["data_location"].SetString("node_historical")
        with self.assertRaisesRegex(ValueError, "data locations"):
            domino_inference_process.Factory(settings, self.model)

        settings = _ProcessSettings(self.scratch)
        settings["Parameters"]["model_type"].SetString("voxel")
        with self.assertRaisesRegex(ValueError, "model type"):
            domino_inference_process.Factory(settings, self.model)


@KratosUnittest.skipUnless(have_domino,
                           "Missing required python modules: torch, warp, physicsnemo.")
class TestDominoThroughProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.scratch = Path("test_domino_inference_scratch")
        self.checkpoint = Path("test_domino_tiny.mdlus")
        self.model = Kratos.Model()
        self.model_part = CreateCaeFixture(self.model)

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.scratch))
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def _SaveTinyDomino(self):
        from physicsnemo.models.domino import DoMINO
        from physicsnemo.models.domino.config import DEFAULT_MODEL_PARAMS

        params = copy.deepcopy(DEFAULT_MODEL_PARAMS)
        # shrink everything: neighbor counts below the 12-triangle surface,
        # narrow layers, tiny background grid; interp_res must match the
        # datapipe's grid_resolution and num_neighbors_surface its
        # num_surface_neighbors
        params["model_type"] = "surface"
        params["interp_res"] = [8, 8, 8]
        params["num_neighbors_surface"] = 4
        params["geometry_rep"]["base_filters"] = 4
        params["geometry_rep"]["geo_conv"]["base_neurons"] = 8
        params["geometry_rep"]["geo_conv"]["surface_neighbors_in_radius"] = [4, 4, 4]
        params["geometry_rep"]["geo_conv"]["volume_neighbors_in_radius"] = [4, 4, 4, 4]
        params["geometry_rep"]["geo_processor"]["base_filters"] = 4
        params["geometry_local"]["base_layer"] = 8
        params["geometry_local"]["surface_neighbors_in_radius"] = [4, 4]
        params["geometry_local"]["volume_neighbors_in_radius"] = [4, 4]
        params["nn_basis_functions"]["base_layer"] = 8
        params["aggregation_model"]["base_layer"] = 8
        params["position_encoder"]["base_neurons"] = 8
        params["parameter_model"]["base_layer"] = 8

        torch.manual_seed(0)
        domino = DoMINO(input_features=3, output_features_vol=None,
                        output_features_surf=1, global_features=2,
                        model_parameters=params)
        domino.save(str(self.checkpoint))

    def test_DeNormalizationIsActuallyApplied(self):
        """Guards the CALL SITE, not the maths.

        _Denormalize is exercised in isolation by TestDominoDenormalization,
        which bypasses __init__ entirely, and by the checkpoint-backed test,
        which calls it directly. Both of those pass with the calls in
        RunInference deleted. This runs the same model twice - once with
        normalization off, once on - and asserts the closed-form relation
        between them, so removing either call fails.
        """
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import domino_inference_process
        class OneChannelFactors:
            max_val = {"surface_fields": numpy.array([3.0])}
            min_val = {"surface_fields": numpy.array([-1.0])}

        self._SaveTinyDomino()

        def Run(extra, scaling=None):
            model = Kratos.Model()
            CreateCaeFixture(model)
            process = domino_inference_process.Factory(
                _ProcessSettings(self.scratch, self.checkpoint, extra), model)
            if scaling is not None:
                process._scaling = scaling      # the .pkl load is not what is under test
            process.ExecuteInitialize()
            model.GetModelPart("Main").ProcessInfo[Kratos.STEP] = 1
            process.ExecuteFinalizeSolutionStep()
            skin = model.GetModelPart("Main").GetSubModelPart("Skin")
            return numpy.array([c.GetValue(Kratos.TEMPERATURE) for c in skin.Conditions])

        raw = Run("")
        denormalized = Run(
            ',"scaling_factors_file":"unused","normalization":"min_max_scaling"',
            scaling=OneChannelFactors())

        # min_max inverse over [-1, 1]: x*(max-min)/2 + (max+min)/2 = 2x + 1
        numpy.testing.assert_allclose(denormalized, 2.0 * raw + 1.0, rtol=1e-9, atol=1e-12)
        # and the two really differ, so the assertion is not vacuous
        self.assertGreater(numpy.abs(denormalized - raw).max(), 1e-6)

    def test_SurfaceInferenceEndToEnd(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import domino_inference_process
        self._SaveTinyDomino()
        process = domino_inference_process.Factory(
            _ProcessSettings(self.scratch, self.checkpoint), self.model)
        process.ExecuteInitialize()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        values = numpy.array([
            condition.GetValue(Kratos.TEMPERATURE)
            for condition in self.model_part.GetSubModelPart("Skin").Conditions])
        self.assertEqual(values.shape, (6,))
        self.assertTrue(numpy.isfinite(values).all())
        self.assertGreater(numpy.abs(values).max(), 0.0)


@KratosUnittest.skipUnless(_CHECKPOINT_DIR is not None,
                           "Needs the public nvidia/domino_drivaerml checkpoint locally "
                           "(set PHYSICSNEMO_DOMINO_CHECKPOINT_DIR).")
@KratosUnittest.skipUnless(have_domino, "Missing torch/physicsnemo/warp.")
class TestDominoRealCheckpointDenormalization(KratosUnittest.TestCase):
    """Against the real pretrained checkpoint's own scaling factors.

    Nothing downloads here: the test self-skips unless the checkpoint is
    already on disk. It exists because the synthetic stand-in used by the
    other tests cannot catch a de-normalization error - its output is
    meaningless either way.
    """

    def setUp(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import domino_inference_process
        self.process = domino_inference_process.DominoInferenceProcess.__new__(
            domino_inference_process.DominoInferenceProcess)
        self.process._scaling = None
        self.process.scaling_factors_file = str(_CHECKPOINT_DIR / "scaling_factors.pkl")
        self.process.normalization = "min_max_scaling"
        self.process.redimensionalize = True
        self.process._global_params_values = {"stream_velocity": [30.0],
                                              "air_density": 1.205}

    def test_ScalingFactorsLoadThroughTheRestrictedUnpickler(self):
        # the public ScalingFactors.load is a plain pickle.load and raises
        # ModuleNotFoundError('utils') on these files
        scaling = self.process._GetScalingFactors()
        maximum = numpy.asarray(scaling.max_val["surface_fields"], dtype=float)
        minimum = numpy.asarray(scaling.min_val["surface_fields"], dtype=float)
        self.assertEqual(maximum.shape, (4,))    # p + 3 wall-shear components
        self.assertTrue(numpy.all(maximum > minimum))

    def test_RawOutputBecomesAPhysicalPressure(self):
        # a representative raw prediction from this checkpoint. Raw, it is a
        # dimensionless ~0.1; the value Kratos should receive is hundreds of
        # pascals, which is the whole point of the fix.
        raw = numpy.array([[0.1386, -0.0657, 0.00083, -0.0055]])
        physical = self.process._Denormalize(raw, "surface_fields")

        pressure = float(physical[0, 0])
        self.assertLess(abs(raw[0, 0]), 1.0)               # raw is dimensionless
        self.assertGreater(abs(pressure), 100.0)           # ...the answer is not
        self.assertLess(abs(pressure), 5000.0)
        # |Cp| of order one against the dynamic pressure 0.5 * rho * U^2
        dynamic_pressure = 0.5 * 1.205 * 30.0 ** 2
        self.assertLess(abs(pressure) / dynamic_pressure, 5.0)

    def test_CheckpointDeclaresItsOwnGridResolution(self):
        # what the load-time guard compares against
        import physicsnemo
        model = physicsnemo.Module.from_checkpoint(
            str(_CHECKPOINT_DIR / "DoMINO.0.501.mdlus"))
        self.assertEqual(list(model.grid_resolution), [128, 64, 64])
        self.assertEqual(model.output_features_surf, 4)


if __name__ == '__main__':
    KratosUnittest.main()
