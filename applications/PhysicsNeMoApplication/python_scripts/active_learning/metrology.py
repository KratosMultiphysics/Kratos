"""MetrologyStrategy backed by this application's validation metrics.

Implements physicsnemo.active_learning's MetrologyStrategy protocol on top
of validation_metrics_process.ComputeMetricValues: each active-learning
iteration, a user-supplied evaluation callable produces named
(predicted, reference) array pairs - e.g. surrogate predictions against
Kratos solves harvested from KratosALSample.fields - and the strategy
records the configured metrics for each pair. Records serialize to the same
JSON shape as ValidationMetricsProcess reports.

physicsnemo is an optional runtime dependency: the strategy class derives
from the physicsnemo protocol class, so it is defined lazily inside the
factory - importing this module never requires physicsnemo.
"""

import json
from pathlib import Path

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.processes.validation_metrics_process import (
    SUPPORTED_METRICS, ComputeMetricValues)

_strategy_class = None


def _TryImportPhysicsNemo():
    try:
        import physicsnemo.active_learning.protocols
        return physicsnemo
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.active_learning.metrology requires physicsnemo, "
            "which could not be imported. Install it with e.g. 'pip install nvidia-physicsnemo'.") from e


def _GetStrategyClass():
    global _strategy_class
    if _strategy_class is not None:
        return _strategy_class

    physicsnemo = _TryImportPhysicsNemo()
    MetrologyStrategy = physicsnemo.active_learning.protocols.MetrologyStrategy

    class ValidationMetricsMetrology(MetrologyStrategy):
        """Records physicsnemo.metrics values for named prediction/reference
        pairs, once per active-learning iteration."""

        __protocol_name__ = "ValidationMetricsMetrology"

        def __init__(self, settings: Kratos.Parameters, evaluation_callable) -> None:
            default_settings = Kratos.Parameters("""{
                "metrics"     : ["mse", "rmse"],
                "output_file" : "metrology_records.json"
            }""")
            settings.ValidateAndAssignDefaults(default_settings)
            self.metrics = settings["metrics"].GetStringArray()
            for metric in self.metrics:
                if metric not in SUPPORTED_METRICS:
                    raise ValueError(
                        f"Unsupported metric \"{metric}\". Supported: {', '.join(SUPPORTED_METRICS)}.")
            self.output_file = Path(settings["output_file"].GetString())
            self._evaluation_callable = evaluation_callable
            self._driver = None
            self.records = []

        def attach(self, other) -> None:
            self._driver = other

        @property
        def is_attached(self) -> bool:
            return self._driver is not None

        def compute(self, *args, **kwargs) -> None:
            import torch  # physicsnemo guarantees torch is present
            pairs = self._evaluation_callable()
            record = {"iteration": len(self.records)}
            for name, (predicted, reference) in pairs.items():
                record[name] = ComputeMetricValues(
                    torch.as_tensor(predicted), torch.as_tensor(reference), self.metrics)
                Kratos.Logger.PrintInfo(
                    self.__protocol_name__,
                    f"iteration {record['iteration']}: {name}: " +
                    ", ".join(f"{k}={v:.6e}" for k, v in record[name].items()))
            self.append(record)

        def reset(self) -> None:
            self.records = []

        def serialize_records(self, path=None, *args, **kwargs) -> None:
            path = Path(path) if path is not None else self.output_file
            path.parent.mkdir(parents=True, exist_ok=True)
            with open(path, "w") as f:
                json.dump(self.records, f, indent=4)

        def load_records(self, path=None, *args, **kwargs) -> None:
            path = Path(path) if path is not None else self.output_file
            with open(path) as f:
                self.records = json.load(f)

    _strategy_class = ValidationMetricsMetrology
    return _strategy_class


def CreateValidationMetricsMetrology(settings: Kratos.Parameters, evaluation_callable):
    """Creates a ValidationMetricsMetrology strategy (lazily importing physicsnemo).

    Args:
        settings: Parameters with "metrics" (subset of
            validation_metrics_process.SUPPORTED_METRICS) and "output_file"
            (JSON records path, also the serialize_records default).
        evaluation_callable: callable () -> {name: (predicted, reference)}
            of same-shape array-likes; called once per compute(). Typically
            it runs the current surrogate on held-out cases and pairs the
            result with reference fields from KratosALSample.fields (keys
            "<VARIABLE>__<location>") or stored solver runs.

    Returns:
        A MetrologyStrategy instance ready for physicsnemo's
        StrategiesConfig(metrology_strategies=[...]).
    """
    return _GetStrategyClass()(settings, evaluation_callable)
