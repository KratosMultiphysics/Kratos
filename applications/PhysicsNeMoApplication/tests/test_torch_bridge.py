import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestTorchBridge(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test")
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        for i in range(10):
            node = self.model_part.CreateNewNode(i + 1, i, i + 1.0, i + 2.0)
            node.SetSolutionStepValue(Kratos.VELOCITY, [i + 1.0, 2.0 * i, 3.0 * i])

    def test_KratosTensorToTorchValues(self):
        ta = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY)
        ta.CollectData()
        tensor = torch_bridge.KratosTensorToTorch(ta)

        self.assertEqual(list(tensor.shape), [10, 3])
        self.assertEqual(tensor.dtype, torch.float64)
        self.assertTrue(numpy.allclose(tensor.numpy(), ta.data))

    def test_KratosTensorToTorchZeroCopy(self):
        ta = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY)
        ta.CollectData()
        tensor = torch_bridge.KratosTensorToTorch(ta)

        # The tensor aliases the adaptor's staging buffer.
        tensor[0, 0] = -42.0
        self.assertEqual(ta.data[0, 0], -42.0)

    def test_TorchToKratosTensorRoundTrip(self):
        ta = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY)
        ta.CollectData()

        tensor = torch_bridge.KratosTensorToTorch(ta).clone()
        tensor += 2.0
        torch_bridge.TorchToKratosTensor(tensor, ta)

        for i, node in enumerate(self.model_part.Nodes):
            velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
            self.assertVectorAlmostEqual(velocity, [i + 3.0, 2.0 * i + 2.0, 3.0 * i + 2.0])

    def test_TorchToKratosTensorShapeMismatch(self):
        ta = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY)
        ta.CollectData()

        with self.assertRaisesRegex(ValueError, "Shape mismatch"):
            torch_bridge.TorchToKratosTensor(torch.zeros(5, 3, dtype=torch.float64), ta)


class TestTorchBridgeWithoutTorch(KratosUnittest.TestCase):
    def test_ModuleImportsWithoutTorch(self):
        # Importing the bridge module (already done at the top of this file)
        # must never require torch; only calling its functions may.
        if have_torch:
            self.skipTest("torch is installed; the lazy-import error path is not reachable.")
        with self.assertRaisesRegex(ImportError, "pip install torch"):
            torch_bridge.KratosTensorToTorch(None)


if __name__ == '__main__':
    KratosUnittest.main()
