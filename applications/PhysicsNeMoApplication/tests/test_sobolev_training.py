"""Derivative-informed training against stored adjoint gradients.

The claim being tested is the reason the feature exists: a surrogate fitted
on values alone has whatever derivatives the fit left behind, and grading it
on Kratos's exact dJ/dX makes those derivatives right. So the central test is
a comparison - the same architecture, the same seed, the same data, trained
with and without the term - and the assertion is on the ratio of the two
gradient errors at points neither run saw.
"""

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.training import sobolev_training
from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False


def _Field(points):
    """u(x, y, z) = x^2 + 2y + 3xz, and its exact gradient."""
    x, y, z = points[:, 0:1], points[:, 1:2], points[:, 2:3]
    value = x * x + 2.0 * y + 3.0 * x * z
    gradient = numpy.concatenate([2.0 * x + 3.0 * z, 2.0 * numpy.ones_like(y), 3.0 * x], axis=1)
    return value, gradient


def _Sample(n_points, seed):
    generator = numpy.random.default_rng(seed)
    points = generator.uniform(-1.0, 1.0, size=(n_points, 3))
    value, gradient = _Field(points)
    return (torch.tensor(points, dtype=torch.float64),
            torch.tensor(numpy.concatenate([value, gradient], axis=1), dtype=torch.float64))


def _MakeModel(seed=0):
    torch.manual_seed(seed)
    return torch.nn.Sequential(
        torch.nn.Linear(3, 32), torch.nn.Tanh(),
        torch.nn.Linear(32, 32), torch.nn.Tanh(),
        torch.nn.Linear(32, 1)).double()


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestSensitivityGradient(KratosUnittest.TestCase):
    """The surrogate's own dJ/dX, on its own."""

    def test_LinearModelGivesItsWeights(self):
        # J = sum(prediction), so dJ/dx of a linear layer is its weight row -
        # closed form, so the sign and the layout are pinned, not eyeballed
        model = torch.nn.Linear(3, 1, bias=False).double()
        with torch.no_grad():
            model.weight.copy_(torch.tensor([[2.0, -3.0, 5.0]], dtype=torch.float64))
        inputs = torch.rand(7, 3, dtype=torch.float64)

        gradient = sobolev_training.SensitivityGradient(model, inputs)
        self.assertEqual(list(gradient.shape), [7, 3])
        for row in range(7):
            self.assertAlmostEqual(float(gradient[row, 0]), 2.0, places=12)
            self.assertAlmostEqual(float(gradient[row, 1]), -3.0, places=12)
            self.assertAlmostEqual(float(gradient[row, 2]), 5.0, places=12)

    def test_ExtraChannelsAreCarriedThroughUndifferentiated(self):
        # only the leading coordinate_channels columns are the coordinates;
        # a feature column must reach the model but never be differentiated
        model = torch.nn.Linear(4, 1, bias=False).double()
        with torch.no_grad():
            model.weight.copy_(torch.tensor([[1.0, 1.0, 1.0, 9.0]], dtype=torch.float64))
        inputs = torch.rand(5, 4, dtype=torch.float64)

        gradient = sobolev_training.SensitivityGradient(model, inputs, coordinate_channels=3)
        self.assertEqual(list(gradient.shape), [5, 3])
        self.assertAlmostEqual(float(gradient.abs().max()), 1.0, places=12)

    def test_ObjectiveChannelsAndWeightsSelect(self):
        model = torch.nn.Linear(2, 2, bias=False).double()
        with torch.no_grad():
            model.weight.copy_(torch.tensor([[1.0, 0.0], [0.0, 4.0]], dtype=torch.float64))
        inputs = torch.rand(3, 2, dtype=torch.float64)

        both = sobolev_training.SensitivityGradient(model, inputs, coordinate_channels=2)
        self.assertAlmostEqual(float(both[0, 0]), 1.0, places=12)
        self.assertAlmostEqual(float(both[0, 1]), 4.0, places=12)

        second = sobolev_training.SensitivityGradient(
            model, inputs, coordinate_channels=2, objective_channels=(1,))
        self.assertAlmostEqual(float(second[0, 0]), 0.0, places=12)
        self.assertAlmostEqual(float(second[0, 1]), 4.0, places=12)

        weighted = sobolev_training.SensitivityGradient(
            model, inputs, coordinate_channels=2, objective_channels=(1,),
            objective_weights=(0.5,))
        self.assertAlmostEqual(float(weighted[0, 1]), 2.0, places=12)

    def test_BatchedLayoutKeepsSamplesApart(self):
        # (B, N, C): summing the objective over the batch is only legitimate
        # because the model maps samples independently - assert it does not
        # leak, by giving two very different batch entries
        model = torch.nn.Linear(3, 1, bias=False).double()
        with torch.no_grad():
            model.weight.copy_(torch.tensor([[1.0, 0.0, 0.0]], dtype=torch.float64))
        inputs = torch.zeros(2, 4, 3, dtype=torch.float64)
        inputs[1] = 100.0
        gradient = sobolev_training.SensitivityGradient(model, inputs)
        self.assertEqual(list(gradient.shape), [2, 4, 3])
        self.assertAlmostEqual(float(gradient[0, 0, 0]), 1.0, places=12)
        self.assertAlmostEqual(float(gradient[1, 0, 0]), 1.0, places=12)


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestSensitivityLossTerm(KratosUnittest.TestCase):
    """The loss term's own contracts."""

    def _Term(self, json_text='{}'):
        return sobolev_training.MakeSensitivityLossTerm(Kratos.Parameters(json_text))

    def test_TermTakesFourArgumentsAndIsSeenAsSuchByTrainModel(self):
        term = self._Term()
        self.assertTrue(training_utils._WantsTargets(term))

    def test_ThreeArgumentTermsAreUnaffected(self):
        # the back-compat pin: every existing loss term keeps its contract
        def OldStyle(model, inputs, prediction):
            return prediction.sum() * 0.0
        self.assertFalse(training_utils._WantsTargets(OldStyle))

        seen = []

        def Counting(model, inputs, prediction):
            seen.append(3)
            return prediction.sum() * 0.0

        model = torch.nn.Linear(3, 1).double()
        inputs, targets = _Sample(16, seed=1)
        dataset = torch.utils.data.TensorDataset(inputs, targets[:, :1])
        training_utils.TrainModel(model, dataset, Kratos.Parameters("""{
            "epochs" : 2, "batch_size" : 8, "device" : "cpu", "seed" : 0
        }"""), extra_loss_terms=[Counting])
        self.assertEqual(len(seen), 4)  # 2 epochs x 2 batches, none refused

    def test_ZeroTermAtTheExactGradient(self):
        # a model whose gradient IS the reference must score exactly zero
        model = torch.nn.Linear(3, 1, bias=False).double()
        with torch.no_grad():
            model.weight.copy_(torch.tensor([[1.0, 2.0, 3.0]], dtype=torch.float64))
        inputs = torch.rand(6, 3, dtype=torch.float64)
        targets = torch.zeros(6, 4, dtype=torch.float64)
        targets[:, 1:] = torch.tensor([1.0, 2.0, 3.0], dtype=torch.float64)

        term = self._Term('{ "gradient_columns" : [1, 2, 3] }')
        self.assertAlmostEqual(float(term(model, inputs, model(inputs), targets)), 0.0,
                               places=14)

    def test_DefaultGradientColumnsAreTheTrailingOnes(self):
        model = torch.nn.Linear(3, 1, bias=False).double()
        with torch.no_grad():
            model.weight.copy_(torch.tensor([[1.0, 2.0, 3.0]], dtype=torch.float64))
        inputs = torch.rand(6, 3, dtype=torch.float64)
        targets = torch.zeros(6, 4, dtype=torch.float64)
        targets[:, -3:] = torch.tensor([1.0, 2.0, 3.0], dtype=torch.float64)
        term = self._Term()
        self.assertAlmostEqual(float(term(model, inputs, model(inputs), targets)), 0.0,
                               places=14)

    def test_RelativeReductionIsScaleFree(self):
        inputs = torch.rand(4, 3, dtype=torch.float64)

        def Score(scale, reduction):
            model = torch.nn.Linear(3, 1, bias=False).double()
            with torch.no_grad():
                model.weight.copy_(
                    scale * torch.tensor([[1.0, 0.0, 0.0]], dtype=torch.float64))
            targets = torch.zeros(4, 4, dtype=torch.float64)
            targets[:, 1:] = scale * torch.tensor([2.0, 0.0, 0.0], dtype=torch.float64)
            term = self._Term(
                '{ "reduction" : "%s", "gradient_columns" : [1, 2, 3] }' % reduction)
            return float(term(model, inputs, model(inputs), targets))

        # the model is wrong by the SAME factor at both scales, so a
        # scale-free score must not move - and the plain mse must, by
        # scale^2. Shape gradients routinely sit orders of magnitude below
        # the field itself, which is what the relative reduction is for.
        self.assertAlmostEqual(Score(1.0, "relative"), Score(1000.0, "relative"), places=10)
        self.assertAlmostEqual(Score(1000.0, "mse") / Score(1.0, "mse"), 1e6, delta=1.0)

    def test_TermIsDifferentiableThroughTheWeights(self):
        # create_graph=True is what makes the matching term trainable at all
        model = torch.nn.Linear(3, 1, bias=False).double()
        inputs = torch.rand(5, 3, dtype=torch.float64)
        targets = torch.zeros(5, 4, dtype=torch.float64)
        loss = self._Term()(model, inputs, model(inputs), targets)
        loss.backward()
        self.assertIsNotNone(model.weight.grad)
        self.assertGreater(float(model.weight.grad.abs().max()), 0.0)

    def test_SettingsAreValidated(self):
        with self.assertRaisesRegex(ValueError, "reduction"):
            self._Term('{ "reduction" : "vibes" }')
        with self.assertRaisesRegex(ValueError, "coordinate_channels"):
            self._Term('{ "coordinate_channels" : 0 }')
        with self.assertRaisesRegex(ValueError, "gradient_columns"):
            self._Term('{ "gradient_columns" : [0, 1] }')
        with self.assertRaisesRegex(ValueError, "objective_weights"):
            self._Term('{ "objective_channels" : [0, 1], "objective_weights" : [1.0] }')


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestSobolevTrainingImprovesGradients(KratosUnittest.TestCase):
    """The claim: gradient supervision buys gradient accuracy."""

    _TRAINING = """{
        "epochs"          : 400,
        "batch_size"      : 32,
        "learning_rate"   : 0.01,
        "device"          : "cpu",
        "seed"            : 0,
        "target_channels" : [0]
    }"""

    def _Train(self, terms):
        inputs, targets = _Sample(64, seed=7)
        dataset = torch.utils.data.TensorDataset(inputs, targets)
        model = _MakeModel(seed=3)
        training_utils.TrainModel(model, dataset, Kratos.Parameters(self._TRAINING),
                                  extra_loss_terms=terms)
        return model

    @staticmethod
    def _Errors(model):
        inputs, targets = _Sample(256, seed=99)   # points neither run saw
        with torch.no_grad():
            value_error = float((model(inputs) - targets[:, :1]).square().mean().sqrt())
        gradient = sobolev_training.SensitivityGradient(model, inputs)
        gradient_error = float((gradient - targets[:, 1:]).square().mean().sqrt())
        return value_error, gradient_error

    def test_GradientErrorDropsAgainstAValueOnlyFit(self):
        plain_value, plain_gradient = self._Errors(self._Train(None))
        term = sobolev_training.MakeSensitivityLossTerm(Kratos.Parameters("""{
            "gradient_columns" : [1, 2, 3], "weight" : 1.0 }"""))
        sobolev_value, sobolev_gradient = self._Errors(self._Train([term]))

        # both fits must be real fits, not two failures compared
        self.assertLess(plain_value, 0.2)
        self.assertLess(sobolev_value, 0.2)
        # the point of the exercise
        self.assertLess(sobolev_gradient, 0.5 * plain_gradient)

    def test_TargetChannelsIsWhatKeepsTheDataLossOnTheValue(self):
        # The trap "target_channels" exists for: torch does NOT reject a
        # 1-channel prediction against a 4-column target, it BROADCASTS -
        # so leaving the setting out trains the model against the mean of
        # the value and the three gradient columns, silently, and the fit
        # it reports looks like a fit. Pinned as a measured difference.
        inputs, targets = _Sample(64, seed=7)
        dataset = torch.utils.data.TensorDataset(inputs, targets)

        restricted = _MakeModel(seed=3)
        training_utils.TrainModel(restricted, dataset, Kratos.Parameters(self._TRAINING))

        broadcast_settings = Kratos.Parameters(self._TRAINING)
        broadcast_settings.RemoveValue("target_channels")
        broadcast = _MakeModel(seed=3)
        training_utils.TrainModel(broadcast, dataset, broadcast_settings)

        self.assertLess(self._Errors(restricted)[0], 0.2)
        self.assertGreater(self._Errors(broadcast)[0], 2.0 * self._Errors(restricted)[0])


if __name__ == '__main__':
    KratosUnittest.main()
