"""Tests for the diffusion bridge on physicsnemo 2.2's PROTOCOL API: the
adapters, MSEDSMLoss training, samplers.sample, DPS guidance (including the
exact discrete FEM residual as the consistency operator) and
multi-diffusion patching."""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.training import diffusion_utils

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import physicsnemo.diffusion.samplers  # noqa: F401
    from physicsnemo.diffusion.noise_schedulers import EDMNoiseScheduler  # noqa: F401
    from physicsnemo.diffusion.guidance import DPSScorePredictor  # noqa: F401
    have_protocol = True
except ImportError:
    have_protocol = False

try:
    import physicsnemo.diffusion.multi_diffusion  # noqa: F401
    have_multi_diffusion = True
except ImportError:
    have_multi_diffusion = False

have_convection_diffusion = kratos_utils.CheckIfApplicationsAvailable(
    "ConvectionDiffusionApplication", "LinearSolversApplication")

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))


def _TinyPreconditioner(seed=0, in_channels=1):
    from physicsnemo.diffusion.preconditioners import EDMPrecondSuperResolution

    torch.manual_seed(seed)
    return EDMPrecondSuperResolution(
        img_resolution=8, img_in_channels=in_channels, img_out_channels=1,
        model_type="SongUNet", model_channels=8, channel_mult=[1, 1],
        num_blocks=1, attn_resolutions=[])


class _ShrinkToCondition(torch.nn.Module if have_torch else object):
    """An analytic x0-predictor: the EDM posterior mean of a latent whose
    prior is the condition. Deterministic, so a guidance term's effect is
    the only thing a test measures."""

    def __init__(self, scheduler, sigma_data=0.5):
        super().__init__()
        self._scheduler = scheduler
        self._sigma_data = sigma_data
        self.img_out_channels = 1
        self.weight = torch.nn.Parameter(torch.zeros(1))  # gives it a device/dtype

    def forward(self, x, t, condition=None, **kwargs):
        sigma = self._scheduler.sigma(t).reshape(-1, *([1] * (x.ndim - 1)))
        share = self._sigma_data ** 2 / (self._sigma_data ** 2 + sigma ** 2)
        return share * x + (1.0 - share) * condition


@KratosUnittest.skipUnless(have_torch and have_protocol,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestWrapDiffusionModel(KratosUnittest.TestCase):
    def test_LegacyAdapterOnlyReordersArguments(self):
        """The whole difference between the two contracts is that the
        condition is positional-and-second in one and a keyword in the
        other - so the adapter must agree with the raw call exactly."""
        from physicsnemo.diffusion import DiffusionModel

        net = _TinyPreconditioner().eval()
        adapted = diffusion_utils.WrapDiffusionModel(net, "edm")
        self.assertIsInstance(adapted, DiffusionModel)
        self.assertEqual(adapted.img_out_channels, 1)

        x = torch.randn(2, 1, 8, 8)
        condition = torch.randn(2, 1, 8, 8)
        t = torch.full((2,), 0.7)
        with torch.no_grad():
            reference = net(x, condition, t)
            through_adapter = adapted(x, t, condition=condition)
        torch.testing.assert_close(through_adapter, reference)

    def test_ConditionIsMandatory(self):
        adapted = diffusion_utils.WrapDiffusionModel(_TinyPreconditioner(), "edm")
        with self.assertRaisesRegex(ValueError, "condition"):
            adapted(torch.randn(1, 1, 8, 8), torch.full((1,), 0.5))

    def test_DoublePreconditioningRefused(self):
        """An EDM-preconditioned denoiser wrapped in EDMPreconditioner
        would be preconditioned twice - silently, with finite garbage."""
        with self.assertRaisesRegex(ValueError, "twice"):
            diffusion_utils.WrapDiffusionModel(
                _TinyPreconditioner(), "edm", preconditioner="edm")

    def test_UnknownInterfaceAndPreconditionerRaise(self):
        with self.assertRaisesRegex(ValueError, "denoiser interface"):
            diffusion_utils.WrapDiffusionModel(_TinyPreconditioner(), "songunet")
        with self.assertRaisesRegex(ValueError, "preconditioner"):
            diffusion_utils.WrapDiffusionModel(
                _TinyPreconditioner(), "edm", preconditioner="vp")

    def test_MissingOutputChannelsRaises(self):
        class NoChannelCount(torch.nn.Module):
            def forward(self, x, t, **kwargs):
                return x

        with self.assertRaisesRegex(ValueError, "out_channels"):
            diffusion_utils.WrapDiffusionModel(NoChannelCount(), "dit")
        # a module that happens to carry the attribute is taken at its word
        self.assertEqual(
            diffusion_utils.WrapDiffusionModel(
                torch.nn.Conv2d(2, 3, 1), "dit").img_out_channels, 3)

    def test_DitAdapterConcatenatesTheCondition(self):
        from physicsnemo.models.dit import DiT

        torch.manual_seed(0)
        dit = DiT(input_size=8, in_channels=2, out_channels=1, patch_size=4,
                  hidden_size=32, depth=1, num_heads=2).eval()
        adapted = diffusion_utils.WrapDiffusionModel(dit, "dit")
        x = torch.randn(2, 1, 8, 8)
        condition = torch.randn(2, 1, 8, 8)
        with torch.no_grad():
            reference = dit(torch.cat([x, condition], dim=1), torch.full((2,), 0.3))
            through_adapter = adapted(x, torch.full((2,), 0.3), condition=condition)
        torch.testing.assert_close(through_adapter, reference)

    def test_PreconditionerWrapsAndKeepsTheChannelCount(self):
        from physicsnemo.models.dit import DiT

        torch.manual_seed(0)
        dit = DiT(input_size=8, in_channels=2, out_channels=1, patch_size=4,
                  hidden_size=32, depth=1, num_heads=2).eval()
        preconditioned = diffusion_utils.WrapDiffusionModel(
            dit, "dit", preconditioner="edm")
        self.assertEqual(preconditioned.img_out_channels, 1)
        with torch.no_grad():
            out = preconditioned(torch.randn(1, 1, 8, 8), torch.full((1,), 0.4),
                                 condition=torch.randn(1, 1, 8, 8))
        self.assertEqual(list(out.shape), [1, 1, 8, 8])

    def test_InterfaceIsNeverSniffedFromTheSignature(self):
        """physicsnemo.models.dit.DiT declares a "condition" parameter, but
        it is a (B, condition_dim) label VECTOR - a conditioning GRID
        belongs in its input channels. A signature probe would call the raw
        DiT a protocol model and feed it the wrong thing, so the interface
        is always explicit."""
        from physicsnemo.models.dit import DiT

        self.assertIn("condition", DiT.forward.__annotations__)
        with self.assertRaisesRegex(ValueError, "denoiser interface"):
            diffusion_utils._AsProtocolModel(_TinyPreconditioner(), "guess")
        # "protocol" is the explicit opt-out, and returns the model itself
        adapted = diffusion_utils.WrapDiffusionModel(_TinyPreconditioner(), "edm")
        self.assertIs(diffusion_utils._AsProtocolModel(adapted, "protocol"), adapted)


@KratosUnittest.skipUnless(have_torch and have_protocol,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestTrainDiffusionModelProtocol(KratosUnittest.TestCase):
    @staticmethod
    def _Dataset(shape=(4, 1, 8, 8), seed=1):
        torch.manual_seed(seed)
        conditions = torch.randn(*shape)
        targets = 2.0 * conditions + 0.1 * torch.randn(*shape)
        return torch.utils.data.TensorDataset(conditions, targets)

    def test_ProtocolTrainingRuns(self):
        history = diffusion_utils.TrainDiffusionModel(
            _TinyPreconditioner(), self._Dataset(), Kratos.Parameters("""{
                "epochs"     : 2,
                "batch_size" : 2,
                "device"     : "cpu",
                "api"        : "protocol",
                "seed"       : 3
            }"""))
        self.assertEqual(len(history), 2)
        self.assertTrue(numpy.isfinite(history).all())

    def test_ApiResolution(self):
        self.assertEqual(diffusion_utils._ResolveDiffusionApi("auto", "edm_sr"), "protocol")
        self.assertEqual(diffusion_utils._ResolveDiffusionApi("auto", "regression"), "legacy")
        self.assertEqual(diffusion_utils._ResolveDiffusionApi("auto", "residual"), "legacy")
        self.assertEqual(diffusion_utils._ResolveDiffusionApi("legacy", "edm_sr"), "legacy")
        with self.assertRaisesRegex(ValueError, "diffusion api"):
            diffusion_utils._ResolveDiffusionApi("protocol_v2", "edm_sr")

    def test_ProtocolWithCorrDiffLossRaises(self):
        """physicsnemo 2.2 has no protocol RegressionLoss/ResidualLoss, so
        asking for one must say that rather than quietly use the legacy."""
        with self.assertRaisesRegex(ValueError, "no protocol-API equivalent"):
            diffusion_utils.TrainDiffusionModel(
                _TinyPreconditioner(), self._Dataset(),
                Kratos.Parameters('{"api": "protocol", "loss": "regression"}'))

    def test_LegacyPathStillTrains(self):
        history = diffusion_utils.TrainDiffusionModel(
            _TinyPreconditioner(), self._Dataset(), Kratos.Parameters("""{
                "epochs"     : 1,
                "batch_size" : 2,
                "device"     : "cpu",
                "api"        : "legacy",
                "seed"       : 3
            }"""))
        self.assertEqual(len(history), 1)
        self.assertTrue(numpy.isfinite(history).all())

    def test_VolumetricTrainsOnTheProtocolAndNotOnTheLegacy(self):
        """MSEDSMLoss is rank-agnostic, so 5-D batches need no local loss
        clone any more; the legacy loss's randn([B,1,1,1]) cannot broadcast
        against them and must say so."""
        class Volumetric(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.convolution = torch.nn.Conv3d(2, 1, 3, padding=1)
                self.img_out_channels = 1

            def forward(self, x, t, condition=None, **kwargs):
                return self.convolution(torch.cat([x, condition], dim=1))

        dataset = self._Dataset(shape=(4, 1, 4, 4, 4))
        history = diffusion_utils.TrainDiffusionModel(
            Volumetric(), dataset, Kratos.Parameters("""{
                "epochs"             : 2,
                "batch_size"         : 2,
                "device"             : "cpu",
                "denoiser_interface" : "protocol",
                "seed"               : 3
            }"""))
        self.assertEqual(len(history), 2)
        self.assertTrue(numpy.isfinite(history).all())

        with self.assertRaisesRegex(ValueError, "rank-agnostic"):
            diffusion_utils.TrainDiffusionModel(
                Volumetric(), dataset,
                Kratos.Parameters('{"api": "legacy", "epochs": 1, "device": "cpu"}'))
        with self.assertRaisesRegex(ValueError, "2D-only"):
            diffusion_utils.TrainDiffusionModel(
                Volumetric(), dataset,
                Kratos.Parameters('{"loss": "regression", "epochs": 1, "device": "cpu"}'))


@KratosUnittest.skipUnless(have_torch and have_protocol,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestGenerateEnsembleProtocol(KratosUnittest.TestCase):
    def test_ShapesAndSeedDeterminism(self):
        model = _TinyPreconditioner().eval()
        settings = """{
            "num_samples" : 2, "num_steps" : 4, "seed" : %d, "api" : "protocol"
        }"""
        first = diffusion_utils.GenerateEnsemble(
            model, numpy.zeros((1, 8, 8)), Kratos.Parameters(settings % 5))
        again = diffusion_utils.GenerateEnsemble(
            model, numpy.zeros((1, 8, 8)), Kratos.Parameters(settings % 5))
        other = diffusion_utils.GenerateEnsemble(
            model, numpy.zeros((1, 8, 8)), Kratos.Parameters(settings % 6))
        self.assertEqual(first.shape, (2, 1, 8, 8))
        self.assertTrue(numpy.isfinite(first).all())
        numpy.testing.assert_allclose(first, again)
        self.assertGreater(numpy.abs(first - other).max(), 0.0)
        self.assertGreater(numpy.abs(first[0] - first[1]).max(), 0.0)

    def test_EverySolverRuns(self):
        model = _TinyPreconditioner().eval()
        for solver in ("heun", "euler", "edm_stochastic_euler", "edm_stochastic_heun"):
            # S_churn belongs to the stochastic solvers; the deterministic
            # ones reject it, so solver_options is per solver
            options = ', "solver_options" : { "S_churn" : 2.0 }' if "stochastic" in solver else ""
            ensemble = diffusion_utils.GenerateEnsemble(
                model, numpy.zeros((1, 8, 8)), Kratos.Parameters("""{
                    "num_samples" : 1,
                    "num_steps"   : 4,
                    "seed"        : 5,
                    "solver"      : "%s"%s
                }""" % (solver, options)))
            self.assertTrue(numpy.isfinite(ensemble).all(), msg=solver)

    def test_UnknownSolverOptionSurfacesUpstreamsError(self):
        with self.assertRaises(TypeError):
            diffusion_utils.GenerateEnsemble(
                _TinyPreconditioner().eval(), numpy.zeros((1, 8, 8)),
                Kratos.Parameters("""{
                    "num_samples"    : 1,
                    "num_steps"      : 4,
                    "solver"         : "heun",
                    "solver_options" : { "S_churn" : 2.0 }
                }"""))

    def test_SingleStepRaises(self):
        """The EDM schedule spans num_steps - 1 intervals, so num_steps = 1
        divides by zero and yields NaN timesteps upstream."""
        with self.assertRaisesRegex(ValueError, "num_steps"):
            diffusion_utils.GenerateEnsemble(
                _TinyPreconditioner().eval(), numpy.zeros((1, 8, 8)),
                Kratos.Parameters('{"num_samples": 1, "num_steps": 1}'))

    def test_LegacyApiRefusesGuidanceAndPatching(self):
        for block in ('"guidance" : { "type" : "data_consistency" }',
                      '"patching" : { "patch_shape" : [4, 4] }'):
            with self.assertRaisesRegex(ValueError, "protocol API"):
                diffusion_utils.GenerateEnsemble(
                    _TinyPreconditioner().eval(), numpy.zeros((1, 8, 8)),
                    Kratos.Parameters('{"api": "legacy", "num_samples": 1, %s}' % block))


@KratosUnittest.skipUnless(have_torch and have_protocol,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestDpsGuidance(KratosUnittest.TestCase):
    def setUp(self):
        from physicsnemo.diffusion.noise_schedulers import EDMNoiseScheduler

        self.scheduler = EDMNoiseScheduler(sigma_data=0.5)
        self.model = _ShrinkToCondition(self.scheduler).eval()
        self.condition = numpy.zeros((1, 8, 8))
        self.observation = numpy.full((1, 8, 8), 3.0)
        self.mask = numpy.zeros((1, 8, 8))
        self.mask[:, :4, :] = 1.0

    def _Sample(self, **kwargs):
        settings = Kratos.Parameters("""{
            "num_samples"        : 1,
            "num_steps"          : 18,
            "seed"               : 0,
            "denoiser_interface" : "protocol",
            "output_channels"    : 1
        }""")
        if "guidance" in kwargs:
            settings.AddValue("guidance", Kratos.Parameters(kwargs["guidance"]))
        return diffusion_utils.GenerateEnsemble(
            self.model, self.condition, settings,
            observation=kwargs.get("observation"), mask=kwargs.get("mask"),
            observation_operator=kwargs.get("observation_operator"))

    def test_DataConsistencyPullsTowardTheObservations(self):
        plain = self._Sample()
        guided = self._Sample(
            guidance='{"type": "data_consistency", "std_y": 0.5}',
            observation=self.observation, mask=self.mask)
        selected = self.mask[None] > 0.5
        plain_error = abs(float(plain[selected].mean()) - 3.0)
        guided_error = abs(float(guided[selected].mean()) - 3.0)
        self.assertLess(guided_error, plain_error)
        # and it is the guidance doing it, not a different random draw
        self.assertGreater(numpy.abs(guided - plain).max(), 0.0)

    def test_GammaRuns(self):
        guided = self._Sample(
            guidance='{"type": "data_consistency", "std_y": 0.5, "gamma": 0.1}',
            observation=self.observation, mask=self.mask)
        self.assertTrue(numpy.isfinite(guided).all())

    def test_TooSmallStdYDivergesAndIsReported(self):
        """The DPS step scales as 1/(2 std_y^2): a std_y far below the
        data's own scale overshoots, and upstream returns the blow-up as
        ordinary finite-then-infinite numbers rather than raising."""
        with self.assertRaisesRegex(ValueError, "std_y"):
            self._Sample(
                guidance='{"type": "data_consistency", "std_y": 1e-4}',
                observation=self.observation, mask=self.mask)

    def test_InvalidGuidanceSettingsRaise(self):
        with self.assertRaisesRegex(ValueError, "std_y"):
            self._Sample(guidance='{"type": "data_consistency", "std_y": 0.0}',
                         observation=self.observation, mask=self.mask)
        with self.assertRaisesRegex(ValueError, "guidance type"):
            self._Sample(guidance='{"type": "l2_projection"}')
        with self.assertRaisesRegex(ValueError, "mask"):
            self._Sample(guidance='{"type": "data_consistency", "std_y": 0.5}',
                         observation=self.observation)
        with self.assertRaisesRegex(ValueError, "mask"):
            self._Sample(guidance='{"type": "data_consistency", "std_y": 0.5}',
                         observation=self.observation,
                         mask=numpy.ones((1, 4, 4)))

    def test_ModelConsistencyPullsTowardTheOperatorsTarget(self):
        def MeanOperator(x0):
            return x0.mean(dim=(2, 3))

        plain = self._Sample()
        guided = self._Sample(
            guidance='{"type": "model_consistency", "std_y": 0.2}',
            observation=numpy.full((1,), 2.0), observation_operator=MeanOperator)
        self.assertLess(abs(float(guided.mean()) - 2.0), abs(float(plain.mean()) - 2.0))

    def test_ModelConsistencyNeedsAnOperator(self):
        with self.assertRaisesRegex(ValueError, "observation_operator"):
            self._Sample(guidance='{"type": "model_consistency", "std_y": 0.2}',
                         observation=numpy.zeros((1,)))


@KratosUnittest.skipUnless(have_torch and have_protocol and have_multi_diffusion,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestMultiDiffusionPatching(KratosUnittest.TestCase):
    """A model trained at patch resolution samples a larger grid."""

    @staticmethod
    def _PatchDataset(samples=4, size=16, seed=0):
        torch.manual_seed(seed)
        conditions = torch.randn(samples, 1, size, size)
        targets = 2.0 * conditions
        return torch.utils.data.TensorDataset(conditions, targets)

    def test_TrainOnPatchesSampleTheWholeGrid(self):
        model = _TinyPreconditioner()  # img_resolution 8, the patch size
        history = diffusion_utils.TrainDiffusionModel(
            model, self._PatchDataset(), Kratos.Parameters("""{
                "epochs"     : 1,
                "batch_size" : 2,
                "device"     : "cpu",
                "seed"       : 0,
                "patching"   : { "patch_shape" : [8, 8], "patch_num" : 2 }
            }"""))
        self.assertTrue(numpy.isfinite(history).all())

        ensemble = diffusion_utils.GenerateEnsemble(
            model.eval(), numpy.zeros((1, 16, 16)), Kratos.Parameters("""{
                "num_samples" : 1,
                "num_steps"   : 4,
                "seed"        : 0,
                "patching"    : { "patch_shape" : [8, 8], "overlap_pix" : 2 }
            }"""))
        self.assertEqual(ensemble.shape, (1, 1, 16, 16))
        self.assertTrue(numpy.isfinite(ensemble).all())

    def test_PatchingRefusedForVolumetricLatents(self):
        with self.assertRaisesRegex(ValueError, "2-D"):
            diffusion_utils.GenerateEnsemble(
                _TinyPreconditioner().eval(), numpy.zeros((1, 4, 4, 4)),
                Kratos.Parameters("""{
                    "num_samples"     : 1,
                    "num_steps"       : 4,
                    "output_channels" : 1,
                    "patching"        : { "patch_shape" : [4, 4] }
                }"""))

    def test_PatchedTrainingRefusesCorrDiffAndVolumetric(self):
        with self.assertRaisesRegex(ValueError, "multi-diffusion"):
            diffusion_utils.TrainDiffusionModel(
                _TinyPreconditioner(), self._PatchDataset(), Kratos.Parameters("""{
                    "epochs"   : 1,
                    "device"   : "cpu",
                    "api"      : "legacy",
                    "patching" : { "patch_shape" : [8, 8], "patch_num" : 2 }
                }"""))

    def test_InvalidPatchGeometrySurfacesUpstreamsError(self):
        with self.assertRaises(Exception):
            diffusion_utils.GenerateEnsemble(
                _TinyPreconditioner().eval(), numpy.zeros((1, 16, 16)),
                Kratos.Parameters("""{
                    "num_samples" : 1,
                    "num_steps"   : 4,
                    "patching"    : { "patch_shape" : [8, 8], "overlap_pix" : 8 }
                }"""))


@KratosUnittest.skipUnless(have_torch and have_protocol and have_convection_diffusion,
                           "Missing torch, physicsnemo (>= 2.2) or "
                           "ConvectionDiffusion/LinearSolvers applications.")
class TestKratosResidualGuidance(KratosUnittest.TestCase):
    """The exact discrete FEM residual as the DPS consistency operator."""

    def setUp(self):
        import thermal_case
        from physicsnemo.diffusion.noise_schedulers import EDMNoiseScheduler

        self.model = Kratos.Model()
        analysis = thermal_case.CreateThermalAnalysis(
            self.model, conductivity=1.0, heat_flux=1.0, divisions=3)
        analysis.Run()
        self.model_part = self.model["ThermalModelPart"]
        self.scheduler = EDMNoiseScheduler(sigma_data=0.5)

    def _Operator(self, **overrides):
        from KratosMultiphysics.PhysicsNeMoApplication.physics import diffusion_residual_operator

        settings = Kratos.Parameters("""{
            "residual_fields" : [ { "variable_name" : "TEMPERATURE" } ],
            "grid_shape"      : [8, 8, 2],
            "squeeze_axis"    : 2
        }""")
        for key, value in overrides.items():
            settings[key] = value
        return diffusion_residual_operator.MakeKratosResidualObservationOperator(
            self.model_part, settings)

    def test_OperatorMatchesTheDirectAssembly(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge
        from KratosMultiphysics.PhysicsNeMoApplication.physics import differentiable_residual

        operator = self._Operator()
        bounding_box = grid_bridge.ComputeBoundingBox(self.model_part)
        grid, _ = grid_bridge.SampleFieldsOnGrid(
            self.model_part, [("TEMPERATURE", "node_historical")], (8, 8, 2), bounding_box)
        squeezed = torch.tensor(grid.mean(axis=3))[None]  # (1, 1, 8, 8)

        residual = operator(squeezed)
        self.assertEqual(tuple(residual.shape), (1, operator._dof_map.n_equations))

        # the same field pushed through the plain autograd Function
        interpolated = grid_bridge.InterpolateGridAtPoints(
            grid, bounding_box,
            numpy.array([[node.X, node.Y, node.Z] for node in self.model_part.Nodes]))
        u = torch.tensor(operator._dof_map.FieldsToDofVector(interpolated))
        u = torch.where(operator._fixed_mask, operator._fixed_values, u)
        reference = differentiable_residual.KratosResidualFunction.Apply(
            u, operator._assembler, operator._dof_map)
        torch.testing.assert_close(residual[0], reference, rtol=1e-10, atol=1e-10)
        operator.Restore()

    def test_TorchInterpolationMatchesNumpyAndCarriesGradients(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge

        bounding_box = grid_bridge.ComputeBoundingBox(self.model_part)
        grid, _ = grid_bridge.SampleFieldsOnGrid(
            self.model_part, [("TEMPERATURE", "node_historical")], (8, 8, 2), bounding_box)
        points = numpy.array([[node.X, node.Y, node.Z] for node in self.model_part.Nodes])

        reference = grid_bridge.InterpolateGridAtPoints(grid, bounding_box, points)
        tensor = torch.tensor(grid, requires_grad=True)
        values = grid_bridge.InterpolateGridAtPointsTorch(tensor, bounding_box, points)
        numpy.testing.assert_allclose(values.detach().numpy(), reference, atol=1e-12)

        values.sum().backward()
        self.assertIsNotNone(tensor.grad)
        # trilinear weights over the whole lattice sum to the point count
        self.assertAlmostEqual(float(tensor.grad.sum()), float(len(points)), places=8)

    def test_TheSolutionHasNoResidualAndAPerturbedFieldDoes(self):
        operator = self._Operator()
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge

        bounding_box = grid_bridge.ComputeBoundingBox(self.model_part)
        grid, _ = grid_bridge.SampleFieldsOnGrid(
            self.model_part, [("TEMPERATURE", "node_historical")], (8, 8, 2), bounding_box)
        solved = torch.tensor(grid.mean(axis=3))[None]
        solved_norm = float(operator(solved).norm())
        perturbed_norm = float(operator(solved + 0.5).norm())
        self.assertLess(solved_norm, perturbed_norm)
        operator.Restore()

    def test_RestorePutsTheModelPartBack(self):
        before = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])
        operator = self._Operator()
        operator(torch.full((1, 1, 8, 8), 123.0, dtype=torch.float64))
        during = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])
        self.assertGreater(numpy.abs(during - before).max(), 1.0)  # assembly wrote the trial
        operator.Restore()
        after = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])
        numpy.testing.assert_allclose(after, before, atol=1e-12)

    def test_EmptyResidualFieldsRaises(self):
        from KratosMultiphysics.PhysicsNeMoApplication.physics import diffusion_residual_operator

        with self.assertRaisesRegex(ValueError, "residual_fields"):
            diffusion_residual_operator.MakeKratosResidualObservationOperator(
                self.model_part, Kratos.Parameters("{}"))

    def test_GuidedSamplingLowersTheResidual(self):
        """The physics grading the generator: guided samples solve the PDE
        better than the same model's unguided ones."""
        operator = self._Operator()
        model = _ShrinkToCondition(self.scheduler).eval()
        condition = numpy.zeros((1, 8, 8))

        def Sample(guidance):
            settings = Kratos.Parameters("""{
                "num_samples"        : 2,
                "num_steps"          : 18,
                "seed"               : 0,
                "denoiser_interface" : "protocol",
                "output_channels"    : 1
            }""")
            settings.AddValue("guidance", Kratos.Parameters(guidance))
            return diffusion_utils.GenerateEnsemble(
                model, condition, settings,
                observation=operator.observation.numpy()[0],
                observation_operator=operator)

        plain = Sample('{"type": "none"}')
        guided = Sample('{"type": "model_consistency", "std_y": 1.0}')
        operator.Restore()

        def ResidualNorm(ensemble):
            return float(operator(torch.tensor(ensemble)).norm(dim=1).mean())

        plain_residual = ResidualNorm(plain)
        guided_residual = ResidualNorm(guided)
        operator.Restore()
        self.assertLess(guided_residual, plain_residual)


@KratosUnittest.skipUnless(have_torch and have_protocol,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestDiffusionInferenceProcessGuidance(KratosUnittest.TestCase):
    """The guidance block reaching the process: observations and the mask
    are read from the model part like any other field."""

    def setUp(self):
        from test_grid_bridge import CreateStructuredTetModelPart

        self.checkpoint = Path("test_diffusion_guidance.mdlus")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE,
                                  Kratos.NODAL_PAUX, Kratos.NODAL_AREA))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.X + node.Y)
            node.SetSolutionStepValue(Kratos.NODAL_PAUX, 3.0)       # the "measurement"
            node.SetSolutionStepValue(
                Kratos.NODAL_AREA, 1.0 if node.X < 0.5 else 0.0)    # where it was taken
        _TinyPreconditioner().save(str(self.checkpoint))

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting(str(self.checkpoint))

    def _Run(self, guidance):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import (
            diffusion_inference_process)

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"  : "Main",
                "model_settings"   : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "input_fields"     : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"    : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "grid_shape"       : [8, 8, 2],
                "squeeze_axis"     : 2,
                "sampler_settings" : {
                    "num_samples" : 2, "num_steps" : 6, "seed" : 7,
                    "guidance"    : %s
                }
            }
        }""" % (self.checkpoint, guidance))
        process = diffusion_inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        return numpy.array([
            (node.GetSolutionStepValue(Kratos.TEMPERATURE),
             node.GetSolutionStepValue(Kratos.NODAL_AREA)) for node in self.model_part.Nodes])

    def test_DataConsistencyReachesTheProcess(self):
        plain = self._Run('{ "type" : "none" }')
        guided = self._Run("""{
            "type"               : "data_consistency",
            "std_y"              : 0.5,
            "observation_fields" : [ { "variable_name" : "NODAL_PAUX", "data_location" : "node_historical" } ],
            "mask_field"         : { "variable_name" : "NODAL_AREA", "data_location" : "node_historical" }
        }""")
        self.assertTrue(numpy.isfinite(guided[:, 0]).all())
        observed = guided[:, 1] > 0.5
        self.assertGreater(observed.sum(), 0)
        plain_error = numpy.abs(plain[observed, 0] - 3.0).mean()
        guided_error = numpy.abs(guided[observed, 0] - 3.0).mean()
        self.assertLess(guided_error, plain_error)

    def test_DataConsistencyWithoutAMaskRaises(self):
        with self.assertRaisesRegex(ValueError, "mask_field"):
            self._Run("""{
                "type"               : "data_consistency",
                "std_y"              : 0.5,
                "observation_fields" : [ { "variable_name" : "NODAL_PAUX", "data_location" : "node_historical" } ]
            }""")

    def test_UnknownGuidanceOperatorRaises(self):
        with self.assertRaisesRegex(ValueError, "guidance operator"):
            self._Run("""{
                "type"     : "model_consistency",
                "std_y"    : 0.5,
                "operator" : "wavelet_projection"
            }""")


@KratosUnittest.skipUnless(have_torch and have_protocol and have_convection_diffusion,
                           "Missing torch, physicsnemo (>= 2.2) or "
                           "ConvectionDiffusion/LinearSolvers applications.")
class TestResidualGuidanceThroughTheProcess(KratosUnittest.TestCase):
    """The whole chain: a process-deployed denoiser steered by the exact
    discrete residual of the very solve it is attached to."""

    def setUp(self):
        import thermal_case

        self.checkpoint = Path("test_diffusion_residual_guidance.mdlus")
        self.model = Kratos.Model()
        analysis = thermal_case.CreateThermalAnalysis(
            self.model, conductivity=1.0, heat_flux=1.0, divisions=3)
        analysis.Run()
        self.model_part = self.model["ThermalModelPart"]
        _TinyPreconditioner().save(str(self.checkpoint))

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting(str(self.checkpoint))

    def test_ModelConsistencyRunsAndLeavesTheSolveIntact(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import (
            diffusion_inference_process)

        before = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"        : "ThermalModelPart",
                "model_settings"         : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "input_fields"           : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "output_fields"          : [ { "variable_name" : "NODAL_PAUX",  "data_location" : "node_non_historical" } ],
                "grid_shape"             : [8, 8, 2],
                "squeeze_axis"           : 2,
                "sampler_settings"       : {
                    "num_samples" : 1, "num_steps" : 6, "seed" : 3,
                    "guidance"    : {
                        "type"            : "model_consistency",
                        "std_y"           : 1.0,
                        "operator"        : "kratos_residual",
                        "residual_fields" : [ { "variable_name" : "TEMPERATURE" } ]
                    }
                }
            }
        }""" % self.checkpoint)
        process = diffusion_inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        written = numpy.array([node.GetValue(Kratos.NODAL_PAUX) for node in self.model_part.Nodes])
        self.assertTrue(numpy.isfinite(written).all())
        # the residual assembly writes trial fields into the DOFs; the
        # process must hand the solve back exactly as it found it
        after = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])
        numpy.testing.assert_allclose(after, before, atol=1e-12)


if __name__ == '__main__':
    KratosUnittest.main()
