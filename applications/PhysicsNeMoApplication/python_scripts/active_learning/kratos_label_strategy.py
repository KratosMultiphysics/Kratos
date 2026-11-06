"""Kratos as the ground-truth solver in a physicsnemo active-learning loop.

Implements physicsnemo.active_learning's LabelStrategy protocol: the driver's
query strategies enqueue KratosALSample instances, and this strategy labels
them by running one Kratos solve per sample through an execution backend.

physicsnemo is an optional runtime dependency: the strategy class derives
from the physicsnemo protocol class, so it is defined lazily inside the
factory below — importing this module never requires physicsnemo.
"""

import KratosMultiphysics as Kratos

_strategy_class = None


def _TryImportPhysicsNemo():
    try:
        import physicsnemo.active_learning.protocols
        return physicsnemo
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.active_learning.kratos_label_strategy requires physicsnemo, "
            "which could not be imported. Install it with e.g. 'pip install nvidia-physicsnemo'.") from e


def _GetStrategyClass():
    global _strategy_class
    if _strategy_class is not None:
        return _strategy_class

    physicsnemo = _TryImportPhysicsNemo()
    LabelStrategy = physicsnemo.active_learning.protocols.LabelStrategy

    class KratosLabelStrategy(LabelStrategy):
        """LabelStrategy running one Kratos solve per queued sample.

        Queue item type: KratosALSample. Failed solves are logged and
        skipped so one diverged case cannot abort the whole labeling phase;
        the number of failures is available as `failed_samples`.
        """

        __protocol_name__ = "KratosLabelStrategy"
        __is_external_process__ = True

        def __init__(self, backend, provides_fields):
            self._backend = backend
            self._driver = None
            self.failed_samples = 0
            self.__is_external_process__ = backend.is_external
            self.__provides_fields__ = set(provides_fields)

        def attach(self, other) -> None:
            self._driver = other

        @property
        def is_attached(self) -> bool:
            return self._driver is not None

        def label(self, queue_to_label, serialize_queue, *args, **kwargs) -> None:
            # Drain the queue into a batch so the backend can fan out
            # (SubprocessBackend runs up to max_parallel_jobs concurrently).
            batch = []
            while not queue_to_label.empty():
                batch.append(queue_to_label.get())
            for sample, result in zip(batch, self._backend.RunCases(batch)):
                if isinstance(result, Exception):
                    self.failed_samples += 1
                    Kratos.Logger.PrintWarning(
                        "KratosLabelStrategy",
                        f"Labeling of sample \"{getattr(sample, 'sample_id', '?')}\" failed and "
                        f"is skipped: {result}")
                    continue
                serialize_queue.put(result)

    _strategy_class = KratosLabelStrategy
    return _strategy_class


def CreateKratosLabelStrategy(backend, provides_fields):
    """Creates a KratosLabelStrategy (lazily importing physicsnemo).

    Args:
        backend: A KratosExecutionBackend (InProcessBackend or
            SubprocessBackend) that solves one case per sample.
        provides_fields: The field keys the backend fills on each sample,
            "<VARIABLE>__<location>" (must match the backend's output specs).

    Returns:
        A LabelStrategy instance ready for physicsnemo's StrategiesConfig.
    """
    return _GetStrategyClass()(backend, provides_fields)
