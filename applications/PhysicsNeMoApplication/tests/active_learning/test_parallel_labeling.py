"""Parallel labeling: SubprocessBackend fan-out and the batched RunCases."""

import os
import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.active_learning.sample_io import KratosALSample
from KratosMultiphysics.PhysicsNeMoApplication.active_learning.execution_backends.base_backend import KratosExecutionBackend
from KratosMultiphysics.PhysicsNeMoApplication.active_learning.execution_backends.subprocess_backend import SubprocessBackend

_AUX_FILES = Path(__file__).parent / "aux_files"


class _FlakyBackend(KratosExecutionBackend):
    """In-memory backend whose RunCase fails for marked samples."""

    def RunCase(self, sample: KratosALSample) -> KratosALSample:
        if sample.parameters.get("fail"):
            raise RuntimeError(f"deliberate failure of {sample.sample_id}")
        sample.fields["X__node_historical"] = numpy.array([1.0])
        return sample

    @property
    def is_external(self) -> bool:
        return False


class TestRunCasesBatch(KratosUnittest.TestCase):
    def test_SerialDefaultReturnsResultsInOrder(self):
        backend = _FlakyBackend()
        samples = [
            KratosALSample("ok_1"),
            KratosALSample("bad", parameters={"fail": True}),
            KratosALSample("ok_2"),
        ]
        results = backend.RunCases(samples)
        self.assertEqual(len(results), 3)
        self.assertTrue(results[0].is_labeled)
        self.assertIsInstance(results[1], RuntimeError)
        self.assertTrue(results[2].is_labeled)


class TestParallelSubprocessBackend(KratosUnittest.TestCase):
    def setUp(self):
        self.working_directory = Path("test_parallel_subprocess_cases")
        self._previous_pythonpath = os.environ.get("PYTHONPATH", "")
        os.environ["PYTHONPATH"] = str(_AUX_FILES) + os.pathsep + self._previous_pythonpath

    def tearDown(self):
        os.environ["PYTHONPATH"] = self._previous_pythonpath
        KratosUtilities.DeleteDirectoryIfExisting(str(self.working_directory))

    def _CreateBackend(self, max_parallel_jobs):
        settings = Kratos.Parameters("""{
            "template_directory" : "",
            "working_directory"  : "test_parallel_subprocess_cases",
            "timeout_seconds"    : 120,
            "max_retries"        : 0,
            "max_parallel_jobs"  : %d
        }""" % max_parallel_jobs)
        settings["template_directory"].SetString(str(_AUX_FILES / "template"))
        settings.AddEmptyArray("run_command")
        for entry in (sys.executable, "MainKratos.py"):
            settings["run_command"].Append(entry)
        return SubprocessBackend(settings)

    def test_ParallelBatchLabelsEverySample(self):
        backend = self._CreateBackend(max_parallel_jobs=2)
        samples = [
            KratosALSample(f"case_{i}", parameters={"dummy_settings/alpha": float(i + 1)})
            for i in range(4)
        ]
        results = backend.RunCases(samples)
        self.assertEqual(len(results), 4)
        for i, result in enumerate(results):
            self.assertNotIsInstance(result, Exception)
            self.assertTrue(numpy.allclose(
                result.fields["PRESSURE__node_historical"], (i + 1.0) * numpy.arange(5.0)))

    def test_ParallelBatchReportsFailuresInPlace(self):
        backend = self._CreateBackend(max_parallel_jobs=2)
        good = KratosALSample("good", parameters={"dummy_settings/alpha": 1.0})
        bad = KratosALSample("bad", parameters={"does/not/exist": 1.0})
        results = backend.RunCases([good, bad])
        self.assertNotIsInstance(results[0], Exception)
        self.assertIsInstance(results[1], Exception)

    def test_InvalidParallelJobCountRaises(self):
        with self.assertRaisesRegex(ValueError, "max_parallel_jobs"):
            self._CreateBackend(max_parallel_jobs=0)


if __name__ == '__main__':
    KratosUnittest.main()
