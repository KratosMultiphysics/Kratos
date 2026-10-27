import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.active_learning import query_strategies

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import physicsnemo.active_learning.protocols
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


class _ListQueue:
    """Minimal AbstractQueue stand-in for exercising sample()."""

    def __init__(self):
        self.items = []

    def put(self, item):
        self.items.append(item)

    def get(self):
        return self.items.pop(0)

    def empty(self):
        return not self.items


def _MakeCandidateSampler(values):
    def sampler(n):
        assert n == len(values)
        return [{"value": float(v)} for v in values]
    return sampler


class TestQueryStrategySelection(KratosUnittest.TestCase):
    """SelectTopCandidates is pure numpy - no ML dependencies."""

    def test_SelectsHighestScoresInOrder(self):
        candidates = ["a", "b", "c", "d"]
        selected = query_strategies.SelectTopCandidates([0.1, 3.0, 2.0, -1.0], candidates, 2)
        self.assertEqual(selected, ["b", "c"])

    def test_StableAmongTies(self):
        selected = query_strategies.SelectTopCandidates([1.0, 1.0, 1.0], ["a", "b", "c"], 2)
        self.assertEqual(selected, ["a", "b"])

    def test_MismatchedLengthsRaise(self):
        with self.assertRaisesRegex(ValueError, "scores for"):
            query_strategies.SelectTopCandidates([1.0], ["a", "b"], 1)


@KratosUnittest.skipUnless(have_torch and have_physicsnemo,
                           "Missing required python modules: torch, physicsnemo.")
class TestEnsembleDisagreementStrategy(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint_dir = Path("test_query_strategy_checkpoints")
        self.checkpoint_dir.mkdir(exist_ok=True)

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.checkpoint_dir))

    def _SaveLinearModel(self, name, slope):
        class Scale(torch.nn.Module):
            def __init__(self, factor):
                super().__init__()
                self.factor = factor

            def forward(self, x):
                return self.factor * x

        path = self.checkpoint_dir / name
        torch.jit.script(Scale(slope)).save(str(path))
        return str(path)

    def _CreateStrategy(self, max_samples=2, pool=4):
        settings = Kratos.Parameters("""{
            "max_samples"          : %d,
            "candidate_pool_size"  : %d,
            "ensemble_checkpoints" : []
        }""" % (max_samples, pool))
        for slope in (1.0, 2.0, 3.0):
            entry = Kratos.Parameters("""{ "checkpoint_file": "", "device": "cpu" }""")
            entry["checkpoint_file"].SetString(self._SaveLinearModel(f"m_{slope}.pt", slope))
            settings["ensemble_checkpoints"].Append(entry)
        return query_strategies.CreateEnsembleDisagreementStrategy(
            settings,
            _MakeCandidateSampler([0.5, -3.0, 1.0, 2.0]),
            lambda candidate: numpy.array([candidate["value"]]))

    def test_HighestDisagreementCandidatesAreEnqueued(self):
        # models predict slope*x -> ensemble variance grows with |x|.
        strategy = self._CreateStrategy()
        queue = _ListQueue()
        strategy.sample(queue)
        self.assertEqual(len(queue.items), 2)
        enqueued_values = [sample.parameters["value"] for sample in queue.items]
        self.assertEqual(enqueued_values, [-3.0, 2.0])
        self.assertTrue(all(sample.sample_id.startswith("query_000_") for sample in queue.items))
        # a second round gets fresh ids
        strategy.sample(queue)
        self.assertTrue(queue.items[-1].sample_id.startswith("query_001_"))

    def test_AttachProtocol(self):
        strategy = self._CreateStrategy()
        self.assertFalse(strategy.is_attached)
        strategy.attach(object())
        self.assertTrue(strategy.is_attached)

    def test_TooFewCheckpointsRaise(self):
        settings = Kratos.Parameters("""{
            "ensemble_checkpoints" : [ { "checkpoint_file": "x.pt" } ]
        }""")
        with self.assertRaisesRegex(ValueError, "at least 2"):
            query_strategies.CreateEnsembleDisagreementStrategy(
                settings, _MakeCandidateSampler([]), lambda c: c)

    def test_PoolSmallerThanMaxSamplesRaises(self):
        settings = Kratos.Parameters("""{
            "max_samples"          : 8,
            "candidate_pool_size"  : 2,
            "ensemble_checkpoints" : [ { "checkpoint_file": "a.pt" }, { "checkpoint_file": "b.pt" } ]
        }""")
        with self.assertRaisesRegex(ValueError, "candidate_pool_size"):
            query_strategies.CreateEnsembleDisagreementStrategy(
                settings, _MakeCandidateSampler([]), lambda c: c)


@KratosUnittest.skipUnless(have_torch and have_physicsnemo,
                           "Missing required python modules: torch, physicsnemo.")
class TestEntropyStrategy(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint_dir = Path("test_entropy_strategy_checkpoints")
        self.checkpoint_dir.mkdir(exist_ok=True)

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.checkpoint_dir))

    def test_LargerActivationsScoreHigherEntropy(self):
        # Dropout on the input: stochastic-pass variance scales with x^2, so
        # the entropy proxy must rank by |x|.
        model = torch.nn.Sequential(torch.nn.Dropout(p=0.5))
        path = self.checkpoint_dir / "dropout.pt"
        torch.jit.script(model).save(str(path))

        settings = Kratos.Parameters("""{
            "max_samples"           : 2,
            "candidate_pool_size"   : 3,
            "model_settings"        : { "checkpoint_file": "%s", "device": "cpu" },
            "num_stochastic_passes" : 64
        }""" % path)
        torch.manual_seed(42)
        strategy = query_strategies.CreateEntropyStrategy(
            settings,
            _MakeCandidateSampler([0.1, 10.0, 1.0]),
            lambda candidate: numpy.full((8,), candidate["value"], dtype=numpy.float32))

        queue = _ListQueue()
        strategy.sample(queue)
        self.assertEqual(len(queue.items), 2)
        self.assertEqual(queue.items[0].parameters["value"], 10.0)


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestSolverResidualStrategy(KratosUnittest.TestCase):
    def test_HighestResidualCandidatesAreEnqueued(self):
        settings = Kratos.Parameters("""{
            "max_samples"         : 2,
            "candidate_pool_size" : 4
        }""")
        strategy = query_strategies.CreateSolverResidualStrategy(
            settings,
            _MakeCandidateSampler([1.0, -4.0, 2.0, 0.5]),
            lambda candidate: candidate["value"] ** 2)

        queue = _ListQueue()
        strategy.sample(queue)
        self.assertEqual([s.parameters["value"] for s in queue.items], [-4.0, 2.0])


if __name__ == '__main__':
    KratosUnittest.main()
