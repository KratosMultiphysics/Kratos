"""Abstract execution backend running one Kratos case per AL sample.

Pure Python: no torch/physicsnemo imports anywhere in the backends.
"""

from abc import ABC, abstractmethod

from KratosMultiphysics.PhysicsNeMoApplication.active_learning.sample_io import KratosALSample


class KratosExecutionBackend(ABC):
    """Runs one Kratos analysis per queued sample and fills its fields."""

    @abstractmethod
    def RunCase(self, sample: KratosALSample) -> KratosALSample:
        """Solves the case described by sample.parameters and returns the
        sample with sample.fields populated."""

    def RunCases(self, samples) -> list:
        """Solves a batch of samples and returns one result per sample, in
        submission order: the labeled sample on success, the raised
        exception on failure (never raises itself, so one diverged case
        cannot abort a batch).

        The base implementation is serial; backends that can fan out
        (e.g. SubprocessBackend with max_parallel_jobs > 1) override it.
        """
        results = []
        for sample in samples:
            try:
                results.append(self.RunCase(sample))
            except Exception as error:  # noqa: BLE001 - reported to the caller
                results.append(error)
        return results

    @property
    @abstractmethod
    def is_external(self) -> bool:
        """True when the solve runs in an external OS process."""
