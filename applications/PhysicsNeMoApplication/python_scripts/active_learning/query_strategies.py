"""Query strategies selecting which cases Kratos should label next.

Implements physicsnemo.active_learning's QueryStrategy protocol three ways:

- **Ensemble disagreement**: variance across an ensemble of checkpoints -
  candidates the models disagree on are the ones worth solving.
- **Predictive entropy**: repeated stochastic forward passes of one model
  (dropout kept active), scored by a Gaussian differential-entropy proxy.
- **Solver residual**: the PDE residual of the surrogate's prediction,
  assembled by the real solver (see solver_residuals.py) - the physics
  itself decides where the surrogate is worst.

Each strategy draws a candidate pool from a user-supplied sampler, scores
it, and enqueues the top ``max_samples`` candidates as KratosALSample items
for the label strategy (kratos_label_strategy.py) to solve.

physicsnemo is an optional runtime dependency: the strategy classes derive
from the physicsnemo protocol class, so they are defined lazily inside the
factories - importing this module never requires physicsnemo, and
SelectTopCandidates is usable without it.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.active_learning.sample_io import KratosALSample

_strategy_classes = None


def _TryImportPhysicsNemo():
    try:
        import physicsnemo.active_learning.protocols
        return physicsnemo
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.active_learning.query_strategies requires physicsnemo, "
            "which could not be imported. Install it with e.g. 'pip install nvidia-physicsnemo'.") from e


def SelectTopCandidates(scores, candidates, max_samples):
    """Returns the max_samples highest-scoring candidates, best first.

    Pure numpy (stable order among ties); usable without physicsnemo.
    """
    scores = numpy.asarray(scores, dtype=numpy.float64)
    if scores.shape[0] != len(candidates):
        raise ValueError(
            f"Got {scores.shape[0]} scores for {len(candidates)} candidates.")
    order = numpy.argsort(-scores, kind="stable")[:int(max_samples)]
    return [candidates[int(index)] for index in order]


def _GetStrategyClasses():
    global _strategy_classes
    if _strategy_classes is not None:
        return _strategy_classes

    physicsnemo = _TryImportPhysicsNemo()
    QueryStrategy = physicsnemo.active_learning.protocols.QueryStrategy

    class _KratosQueryStrategyBase(QueryStrategy):
        """Shared plumbing: attachment, candidate scoring loop, enqueueing."""

        __protocol_name__ = "KratosQueryStrategy"

        def __init__(self, settings: Kratos.Parameters, candidate_sampler) -> None:
            settings.ValidateAndAssignDefaults(self._DefaultSettings())
            self.max_samples = settings["max_samples"].GetInt()
            self.candidate_pool_size = settings["candidate_pool_size"].GetInt()
            if self.max_samples < 1 or self.candidate_pool_size < self.max_samples:
                raise ValueError(
                    f"Need candidate_pool_size >= max_samples >= 1 "
                    f"[ candidate_pool_size = {self.candidate_pool_size}, "
                    f"max_samples = {self.max_samples} ].")
            self._settings = settings
            self._candidate_sampler = candidate_sampler
            self._driver = None
            self._iteration = 0

        @classmethod
        def _DefaultSettings(cls) -> Kratos.Parameters:
            return Kratos.Parameters("""{
                "max_samples"         : 4,
                "candidate_pool_size" : 16
            }""")

        def attach(self, other) -> None:
            self._driver = other

        @property
        def is_attached(self) -> bool:
            return self._driver is not None

        def _ScoreCandidates(self, candidates):
            raise NotImplementedError

        def sample(self, query_queue, *args, **kwargs) -> None:
            candidates = list(self._candidate_sampler(self.candidate_pool_size))
            scores = self._ScoreCandidates(candidates)
            selected = SelectTopCandidates(scores, candidates, self.max_samples)
            for rank, candidate in enumerate(selected):
                query_queue.put(KratosALSample(
                    sample_id=f"query_{self._iteration:03d}_{rank:02d}",
                    parameters=dict(candidate)))
            Kratos.Logger.PrintInfo(
                self.__protocol_name__,
                f"Iteration {self._iteration}: enqueued {len(selected)} of "
                f"{len(candidates)} candidates.")
            self._iteration += 1

    class EnsembleDisagreementStrategy(_KratosQueryStrategyBase):
        """Scores candidates by prediction variance across a checkpoint ensemble."""

        __protocol_name__ = "EnsembleDisagreementStrategy"

        def __init__(self, settings: Kratos.Parameters, candidate_sampler, encode_candidate) -> None:
            super().__init__(settings, candidate_sampler)
            self._ensemble_settings = [
                self._settings["ensemble_checkpoints"][i].Clone()
                for i in range(self._settings["ensemble_checkpoints"].size())]
            if len(self._ensemble_settings) < 2:
                raise ValueError("\"ensemble_checkpoints\" needs at least 2 model settings.")
            self._encode_candidate = encode_candidate
            self._models = None

        @classmethod
        def _DefaultSettings(cls) -> Kratos.Parameters:
            return Kratos.Parameters("""{
                "max_samples"          : 4,
                "candidate_pool_size"  : 16,
                "ensemble_checkpoints" : []
            }""")

        def _GetModels(self):
            if self._models is None:
                self._models = [model_registry.LoadModel(s) for s in self._ensemble_settings]
            return self._models

        def _ScoreCandidates(self, candidates):
            import torch  # physicsnemo guarantees torch is present
            models = self._GetModels()
            scores = []
            with torch.no_grad():
                for candidate in candidates:
                    encoded = torch.as_tensor(self._encode_candidate(candidate))
                    predictions = torch.stack([
                        model(encoded.to(device)).cpu().to(torch.float64)
                        for model, device in models])
                    scores.append(float(predictions.var(dim=0, unbiased=False).mean()))
            return scores

    class EntropyStrategy(_KratosQueryStrategyBase):
        """Scores candidates by predictive entropy over stochastic forward
        passes of one model (dropout kept active via model.train())."""

        __protocol_name__ = "EntropyStrategy"

        def __init__(self, settings: Kratos.Parameters, candidate_sampler, encode_candidate) -> None:
            super().__init__(settings, candidate_sampler)
            self._num_passes = self._settings["num_stochastic_passes"].GetInt()
            if self._num_passes < 2:
                raise ValueError("\"num_stochastic_passes\" must be >= 2.")
            self._model_settings = self._settings["model_settings"].Clone()
            self._encode_candidate = encode_candidate
            self._model = None

        @classmethod
        def _DefaultSettings(cls) -> Kratos.Parameters:
            return Kratos.Parameters("""{
                "max_samples"           : 4,
                "candidate_pool_size"   : 16,
                "model_settings"        : {},
                "num_stochastic_passes" : 8
            }""")

        def _GetModel(self):
            if self._model is None:
                self._model = model_registry.LoadModel(self._model_settings)
            return self._model

        def _ScoreCandidates(self, candidates):
            import torch
            model, device = self._GetModel()
            model.train()  # keep dropout layers stochastic
            scores = []
            try:
                with torch.no_grad():
                    for candidate in candidates:
                        encoded = torch.as_tensor(self._encode_candidate(candidate)).to(device)
                        passes = torch.stack(
                            [model(encoded).cpu().to(torch.float64) for _ in range(self._num_passes)])
                        variance = passes.var(dim=0, unbiased=False)
                        # Gaussian differential-entropy proxy, mean over outputs
                        entropy = 0.5 * torch.log(
                            2.0 * torch.pi * torch.e * variance + 1e-30)
                        scores.append(float(entropy.mean()))
            finally:
                model.eval()
            return scores

    class SolverResidualStrategy(_KratosQueryStrategyBase):
        """Scores candidates by a user-supplied residual evaluator - typically
        the surrogate's prediction written into a Kratos model part and
        assembled by solver_residuals.ResidualEvaluator."""

        __protocol_name__ = "SolverResidualStrategy"

        def __init__(self, settings: Kratos.Parameters, candidate_sampler, residual_evaluator) -> None:
            super().__init__(settings, candidate_sampler)
            self._residual_evaluator = residual_evaluator

        def _ScoreCandidates(self, candidates):
            return [float(self._residual_evaluator(candidate)) for candidate in candidates]

    _strategy_classes = {
        "ensemble": EnsembleDisagreementStrategy,
        "entropy": EntropyStrategy,
        "residual": SolverResidualStrategy,
    }
    return _strategy_classes


def CreateEnsembleDisagreementStrategy(settings: Kratos.Parameters, candidate_sampler, encode_candidate):
    """Creates an ensemble-disagreement query strategy.

    Args:
        settings: Parameters with "max_samples", "candidate_pool_size" and
            "ensemble_checkpoints" (array of model_registry.LoadModel
            settings, at least 2 entries).
        candidate_sampler: callable n -> iterable of parameter dicts
            (design-space points, KratosALSample.parameters format).
        encode_candidate: callable dict -> array-like model input for one
            candidate (any shape the ensemble models accept).
    """
    return _GetStrategyClasses()["ensemble"](settings, candidate_sampler, encode_candidate)


def CreateEntropyStrategy(settings: Kratos.Parameters, candidate_sampler, encode_candidate):
    """Creates a predictive-entropy query strategy.

    Args:
        settings: Parameters with "max_samples", "candidate_pool_size",
            "model_settings" (model_registry.LoadModel settings) and
            "num_stochastic_passes" (>= 2). The model should contain dropout
            (or another stochastic layer) for the passes to disagree.
        candidate_sampler / encode_candidate: see
            CreateEnsembleDisagreementStrategy.
    """
    return _GetStrategyClasses()["entropy"](settings, candidate_sampler, encode_candidate)


def CreateSolverResidualStrategy(settings: Kratos.Parameters, candidate_sampler, residual_evaluator):
    """Creates a solver-residual query strategy.

    Args:
        settings: Parameters with "max_samples" and "candidate_pool_size".
        candidate_sampler: callable n -> iterable of parameter dicts.
        residual_evaluator: callable dict -> float scoring one candidate;
            typically: run the surrogate for the candidate, write the
            prediction into the case's model part, and return
            solver_residuals.ResidualEvaluator.ComputeResidualNorm().
    """
    return _GetStrategyClasses()["residual"](settings, candidate_sampler, residual_evaluator)
