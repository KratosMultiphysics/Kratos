"""Fine-tuning a pretrained DoMINO: predictor-corrector and LoRA.

The corrector tests need only torch. The LoRA tests run against a tiny
synthetic DoMINO so they do not need the 48 MB pretrained checkpoint; one
checkpoint-backed test covers the property that actually matters for
deployment - that a merged model is an ordinary .mdlus.
"""

import copy
import os
from pathlib import Path

import numpy

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.training import domino_finetune
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.models.domino import DoMINO  # noqa: F401
    from physicsnemo.experimental.peft import LoRAConfig  # noqa: F401
    have_peft = have_torch
except ImportError:
    have_peft = False

from test_domino_inference_process import _CHECKPOINT_DIR


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestDominoCorrector(KratosUnittest.TestCase):
    """The predictor-corrector decomposition Y = Y_predictor + Y_corrector."""

    def setUp(self):
        torch.manual_seed(0)
        self.rng = numpy.random.default_rng(0)
        self.features = self.rng.random((32, 5))
        self.base = self.rng.random((32, 3))

    def test_UntrainedCorrectorIsExactlyTheIdentityOnThePredictor(self):
        # the last layer is zero-initialized on purpose: fine-tuning starts
        # from the pretrained model's own answer, so it can only improve
        corrector = domino_finetune.CreateCorrector(5, 3, hidden=8, n_layers=2)
        combined = domino_finetune.ApplyCorrector(corrector, self.base, self.features)
        numpy.testing.assert_array_almost_equal(combined, self.base, decimal=12)

    def test_TrainingReducesTheResidual(self):
        corrector = domino_finetune.CreateCorrector(5, 3, hidden=16, n_layers=3)
        # a learnable target: the residual is a fixed linear map of the features
        weights = self.rng.random((5, 3))
        residuals = self.features @ weights
        history = domino_finetune.TrainCorrector(
            corrector, self.features, residuals, epochs=200, learning_rate=1e-2)
        self.assertEqual(len(history), 200)
        self.assertLess(history[-1], 0.2 * history[0])

        combined = domino_finetune.ApplyCorrector(corrector, self.base, self.features)
        # it now moves the prediction toward base + residual
        self.assertLess(numpy.abs(combined - (self.base + residuals)).mean(),
                        numpy.abs(self.base - (self.base + residuals)).mean())

    def test_CorrectorDtypeIsRespectedEndToEnd(self):
        """CreateCorrector takes a dtype, so training and applying must follow it.

        Both helpers used to pin their inputs to float32, so a float64
        corrector - which the dtype argument openly invites - died inside the
        first linear layer with "mat1 and mat2 must have the same dtype".
        """
        for dtype in (torch.float32, torch.float64):
            corrector = domino_finetune.CreateCorrector(
                5, 3, hidden=8, n_layers=2, dtype=dtype)
            self.assertEqual(next(corrector.parameters()).dtype, dtype)

            residuals = self.features @ self.rng.random((5, 3))
            domino_finetune.TrainCorrector(
                corrector, self.features, residuals, epochs=5)
            combined = domino_finetune.ApplyCorrector(
                corrector, self.base, self.features)
            self.assertEqual(combined.shape, self.base.shape)

    def test_MismatchedRowsRejected(self):
        corrector = domino_finetune.CreateCorrector(5, 3)
        with self.assertRaisesRegex(ValueError, "line up"):
            domino_finetune.TrainCorrector(corrector, self.features, self.base[:4])

    def test_MismatchedOutputShapeRejected(self):
        corrector = domino_finetune.CreateCorrector(5, 7)
        with self.assertRaisesRegex(ValueError, "base prediction"):
            domino_finetune.ApplyCorrector(corrector, self.base, self.features)

    def test_CachedPredictionsAreDetachedAndLeaveTrainingModeAlone(self):
        # the predictor must not be perturbed by being cached, and its
        # outputs must not carry a graph into the corrector's optimizer
        class Predictor(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.linear = torch.nn.Linear(3, 2)

            def forward(self, batch):
                return None, self.linear(batch["x"])

        predictor = Predictor()
        predictor.train()
        batches = [{"x": torch.ones(1, 4, 3)}, {"x": torch.zeros(1, 4, 3)}]
        cached = domino_finetune.CacheBasePredictions(predictor, batches)

        self.assertTrue(predictor.training)          # restored
        self.assertEqual(len(cached), 2)
        for volume, surface in cached:
            self.assertIsNone(volume)
            self.assertIsInstance(surface, numpy.ndarray)
            self.assertEqual(surface.shape, (1, 4, 2))


@KratosUnittest.skipUnless(have_peft,
                           "Missing torch/physicsnemo with experimental PEFT.")
class TestDominoLora(KratosUnittest.TestCase):
    """LoRA adapters on a DoMINO, against a tiny synthetic model."""

    @staticmethod
    def _TinyDomino():
        from physicsnemo.models.domino import DoMINO
        from physicsnemo.models.domino.config import DEFAULT_MODEL_PARAMS

        params = copy.deepcopy(DEFAULT_MODEL_PARAMS)
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
        torch.manual_seed(0)
        return DoMINO(input_features=3, output_features_vol=None,
                      output_features_surf=1, global_features=2,
                      model_parameters=params)

    def test_AdaptersWrapLayersAndFreezeTheRest(self):
        model = self._TinyDomino()
        total = sum(p.numel() for p in model.parameters())
        model, wrapped, trainable = domino_finetune.ApplyLora(model, rank=4)
        self.assertGreater(wrapped, 0)
        self.assertGreater(trainable, 0)
        self.assertLess(trainable, total)      # only a fraction is trainable
        frozen = sum(p.numel() for p in model.parameters() if not p.requires_grad)
        self.assertGreater(frozen, trainable)

    def test_APatternThatMatchesNothingRaises(self):
        # the documented trap: solution_calculator_surf re-uses nn_basis_surf
        # and agg_model_surf by reference, so it owns no wrappable layer.
        # Silently training zero parameters looks exactly like success.
        model = self._TinyDomino()
        # physicsnemo raises this itself ("matched 0 wrappable layers"); the
        # bridge keeps a backstop in case that ever moves
        with self.assertRaisesRegex(ValueError, "wrappable"):
            domino_finetune.ApplyLora(model, target_pattern=r"^solution_calculator_surf\..*")

    def test_MergedModelIsAnOrdinaryCheckpoint(self):
        model = self._TinyDomino()
        expected = sum(p.numel() for p in model.parameters())
        model, _, _ = domino_finetune.ApplyLora(model, rank=4)

        path = Path("test_domino_lora_merged.mdlus")
        try:
            domino_finetune.MergeAndSave(model, str(path))
            self.assertTrue(path.is_file())
            import physicsnemo
            reloaded = physicsnemo.Module.from_checkpoint(str(path))
            # merging folds the adapters back in: the parameter count returns
            # to the base model's, which is what makes it deployable unchanged
            self.assertEqual(sum(p.numel() for p in reloaded.parameters()), expected)
            self.assertEqual(list(reloaded.grid_resolution), [8, 8, 8])
        finally:
            KratosUtilities.DeleteFileIfExisting(str(path))


@KratosUnittest.skipUnless(_CHECKPOINT_DIR is not None,
                           "Needs the public nvidia/domino_drivaerml checkpoint locally "
                           "(set PHYSICSNEMO_DOMINO_CHECKPOINT_DIR).")
@KratosUnittest.skipUnless(have_peft, "Missing torch/physicsnemo with experimental PEFT.")
class TestDominoLoraOnRealCheckpoint(KratosUnittest.TestCase):

    def test_LoraIsASmallFractionOfARealCheckpoint(self):
        import physicsnemo
        model = physicsnemo.Module.from_checkpoint(
            str(_CHECKPOINT_DIR / "DoMINO.0.501.mdlus"))
        total = sum(p.numel() for p in model.parameters())
        self.assertGreater(total, 1_000_000)          # a real model, not the stand-in
        model, wrapped, trainable = domino_finetune.ApplyLora(model, rank=8)
        self.assertGreater(wrapped, 0)
        self.assertLess(trainable / total, 0.05)      # a few percent, measured ~1.7 %


if __name__ == '__main__':
    KratosUnittest.main()
