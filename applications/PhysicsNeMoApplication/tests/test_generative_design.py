"""Tests for generative topology design on real compliance data: the SIMP
solver case, and TopoDiff trained on density fields conditioned by where the
structure is held and loaded."""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))

have_structural = kratos_utils.CheckIfApplicationsAvailable(
    "StructuralMechanicsApplication", "LinearSolversApplication")

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.models.topodiff import TopoDiff
    have_topodiff = True
except ImportError:
    have_topodiff = False


@KratosUnittest.skipUnless(have_structural,
                           "Missing StructuralMechanics/LinearSolvers applications.")
class TestComplianceCase(KratosUnittest.TestCase):
    """The ground truth a generative design is scored against."""

    def test_RemovingMaterialRaisesCompliance(self):
        import compliance_case

        solid = numpy.ones(2 * 6 * 6)
        stiff, model_part = compliance_case.SolveCompliance(
            Kratos.Model(), solid, divisions=6)
        self.assertEqual(model_part.NumberOfElements(), 2 * 6 * 6)
        self.assertGreater(abs(stiff), 0.0)

        half = numpy.full_like(solid, 0.5)
        softer, _ = compliance_case.SolveCompliance(Kratos.Model(), half, divisions=6)
        self.assertGreater(abs(softer), abs(stiff))

    def test_TheLoadNeedsItsCondition(self):
        """A nodal POINT_LOAD value alone does nothing: without a point-load
        CONDITION the solve returns zero displacement and a compliance of
        exactly zero, which looks like a converged answer."""
        import compliance_case

        model = Kratos.Model()
        model_part = compliance_case.CreateComplianceModelPart(model, divisions=4)
        compliance_case.ApplyDensities(model_part, numpy.ones(model_part.NumberOfElements()))
        compliance_case.ApplyCaseData(model_part)
        self.assertGreater(model_part.NumberOfConditions(), 0)

    def test_SimpPenalizesIntermediateDensities(self):
        import compliance_case

        model = Kratos.Model()
        model_part = compliance_case.CreateComplianceModelPart(model, divisions=4)
        densities = numpy.linspace(0.2, 1.0, model_part.NumberOfElements())
        compliance_case.ApplyDensities(model_part, densities)
        moduli = numpy.array([
            element.Properties.GetValue(Kratos.YOUNG_MODULUS)
            for element in model_part.Elements])
        # E = E_min + rho^3 (E_0 - E_min): a half-dense element is an eighth
        # as stiff, not half - that is what makes designs go black and white
        self.assertLess(moduli[0], moduli[-1])
        ratio = (moduli[-1] - moduli[0]) / moduli[-1]
        self.assertGreater(ratio, 0.9)
        # every element carries its OWN properties
        self.assertEqual(len({element.Properties.Id for element in model_part.Elements}),
                         model_part.NumberOfElements())

    def test_MismatchedDensityCountRaises(self):
        import compliance_case

        model = Kratos.Model()
        model_part = compliance_case.CreateComplianceModelPart(model, divisions=4)
        with self.assertRaisesRegex(ValueError, "densities"):
            compliance_case.ApplyDensities(model_part, numpy.ones(3))

    def test_TheConditioningChannelsDescribeTheProblem(self):
        import compliance_case

        channels = compliance_case.ConstraintChannels(8, 0.4)
        self.assertEqual(channels.shape, (3, 8, 8))
        # the channels must line up with DensityGrid's layout (first index x,
        # second y), not merely carry the right number of ones: the clamp is
        # the x = 0 row and the load sits at the far x, mid-height
        self.assertEqual(channels[0][0, :].sum(), 8.0)      # the clamped edge
        self.assertEqual(channels[0].sum(), 8.0)            # and nothing else
        self.assertEqual(channels[1].sum(), 1.0)            # one loaded cell
        self.assertEqual(channels[1][-1, 4], 1.0)           # far x, mid-height
        numpy.testing.assert_allclose(channels[2], 0.4)     # the volume fraction

    def test_TheDensityGridMatchesTheElementField(self):
        import compliance_case

        densities = compliance_case.SmoothRandomDensities(6, seed=1)
        grid = compliance_case.DensityGrid(densities, 6)
        self.assertEqual(grid.shape, (6, 6))
        self.assertTrue(set(numpy.unique(grid)).issubset({0.2, 1.0}))


@KratosUnittest.skipUnless(have_structural and have_torch and have_topodiff,
                           "Missing StructuralMechanics, torch or physicsnemo.")
class TestTopoDiffOnComplianceData(KratosUnittest.TestCase):
    """TopoDiff on the real thing: designs conditioned on the supports, the
    load and a volume fraction, scored by re-solving their compliance."""

    def setUp(self):
        self.divisions = 16   # TopoDiff's img_resolution
        self.checkpoint = Path("test_topodiff_design.mdlus")

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting(str(self.checkpoint))

    def _Dataset(self, samples=4):
        import compliance_case

        conditions, targets = [], []
        for index in range(samples):
            densities = compliance_case.SmoothRandomDensities(
                self.divisions, seed=index, threshold=0.2)
            grid = compliance_case.DensityGrid(densities, self.divisions)
            channels = compliance_case.ConstraintChannels(self.divisions, float(grid.mean()))
            conditions.append(channels)
            targets.append(grid[None])
        return (torch.tensor(numpy.stack(conditions), dtype=torch.float32),
                torch.tensor(numpy.stack(targets), dtype=torch.float32))

    def test_TrainingOnDensityFieldsAndSampling(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import diffusion_utils

        conditions, targets = self._Dataset()
        dataset = torch.utils.data.TensorDataset(conditions, targets)
        torch.manual_seed(0)
        model = TopoDiff(img_resolution=self.divisions, in_channels=4, out_channels=1,
                         model_channels=64, channel_mult=[1, 1], num_blocks=1,
                         attn_resolutions=[])
        wrapped = diffusion_utils.WrapDiffusionModel(model, "topodiff", out_channels=1)

        history = diffusion_utils.TrainDiffusionModel(
            wrapped, dataset, Kratos.Parameters("""{
                "epochs"             : 2,
                "batch_size"         : 2,
                "device"             : "cpu",
                "seed"               : 0,
                "denoiser_interface" : "protocol"
            }"""))
        self.assertEqual(len(history), 2)
        self.assertTrue(numpy.isfinite(history).all())

        ensemble = diffusion_utils.GenerateEnsemble(
            wrapped, conditions[0].numpy(), Kratos.Parameters("""{
                "num_samples"        : 2,
                "num_steps"          : 4,
                "seed"               : 0,
                "denoiser_interface" : "protocol",
                "output_channels"    : 1
            }"""))
        self.assertEqual(ensemble.shape, (2, 1, self.divisions, self.divisions))
        self.assertTrue(numpy.isfinite(ensemble).all())

    def test_AGeneratedDesignCanBeScoredByResolvingIt(self):
        """What closes the loop: a sampled density field is a real design,
        so Kratos can solve it and report its compliance."""
        import compliance_case

        rng = numpy.random.default_rng(0)
        design = numpy.clip(rng.uniform(size=(self.divisions, self.divisions)), 0.2, 1.0)
        densities = numpy.repeat(design.reshape(-1), 2)
        compliance, model_part = compliance_case.SolveCompliance(
            Kratos.Model(), densities, divisions=self.divisions)
        self.assertTrue(numpy.isfinite(compliance))
        self.assertGreater(abs(compliance), 0.0)
        self.assertEqual(model_part.NumberOfElements(), 2 * self.divisions ** 2)


if __name__ == '__main__':
    KratosUnittest.main()
