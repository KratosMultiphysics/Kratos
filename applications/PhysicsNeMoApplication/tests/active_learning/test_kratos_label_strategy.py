import queue

import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.active_learning.sample_io import KratosALSample
from KratosMultiphysics.PhysicsNeMoApplication.active_learning.execution_backends.base_backend import KratosExecutionBackend
from KratosMultiphysics.PhysicsNeMoApplication.active_learning import kratos_label_strategy

try:
    import physicsnemo.active_learning.protocols
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


class _StubBackend(KratosExecutionBackend):
    """Labels instantly; raises for samples whose id contains 'fail'."""

    def __init__(self):
        self.run_ids = []

    def RunCase(self, sample):
        if "fail" in sample.sample_id:
            raise RuntimeError("intentional failure")
        self.run_ids.append(sample.sample_id)
        sample.fields["PRESSURE__node_historical"] = [1.0, 2.0]
        return sample

    @property
    def is_external(self):
        return False


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestKratosLabelStrategy(KratosUnittest.TestCase):
    def _CreateStrategy(self, backend):
        return kratos_label_strategy.CreateKratosLabelStrategy(
            backend, provides_fields={"PRESSURE__node_historical"})

    def test_ProtocolContract(self):
        strategy = self._CreateStrategy(_StubBackend())
        self.assertEqual(strategy.__protocol_name__, "KratosLabelStrategy")
        self.assertFalse(strategy.__is_external_process__)  # stub backend is in-process
        self.assertEqual(strategy.__provides_fields__, {"PRESSURE__node_historical"})
        self.assertFalse(strategy.is_attached)
        strategy.attach(object())
        self.assertTrue(strategy.is_attached)

    def test_LabelDrainsQueue(self):
        backend = _StubBackend()
        strategy = self._CreateStrategy(backend)

        to_label, serialize = queue.Queue(), queue.Queue()
        to_label.put(KratosALSample("s1", parameters={"a": 1.0}))
        to_label.put(KratosALSample("s2", parameters={"a": 2.0}))

        strategy.label(to_label, serialize)

        self.assertTrue(to_label.empty())
        self.assertEqual(serialize.qsize(), 2)
        self.assertEqual(backend.run_ids, ["s1", "s2"])
        while not serialize.empty():
            self.assertTrue(serialize.get().is_labeled)

    def test_FailureIsolation(self):
        strategy = self._CreateStrategy(_StubBackend())

        to_label, serialize = queue.Queue(), queue.Queue()
        to_label.put(KratosALSample("s_ok"))
        to_label.put(KratosALSample("s_fail"))
        to_label.put(KratosALSample("s_ok_too"))

        strategy.label(to_label, serialize)

        self.assertEqual(serialize.qsize(), 2)
        self.assertEqual(strategy.failed_samples, 1)


if __name__ == '__main__':
    KratosUnittest.main()
