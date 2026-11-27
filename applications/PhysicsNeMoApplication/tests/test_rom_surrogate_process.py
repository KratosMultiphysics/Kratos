from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from test_rom_bridge import _WriteNumpyBasis, _OrthonormalBasis

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestRomSurrogateProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.basis_folder = Path("test_rom_surrogate_basis")
        self.checkpoint = Path("test_rom_surrogate_model.pt")
        self.model = Kratos.Model()

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.basis_folder))
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint) + ".card.json")

    def _CreateModelPart(self, variables, n_nodes=4):
        model_part = self.model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(Kratos.NODAL_PAUX)
        for variable in variables:
            model_part.AddNodalSolutionStepVariable(variable)
        for i in range(n_nodes):
            node = model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.NODAL_PAUX, 2.0)  # the parameter carrier
        return model_part

    def _SaveLinearModel(self, weight_rows):
        """q = W @ input, saved as TorchScript (weight_rows: (n_modes, C_in))."""
        class ToModes(torch.nn.Module):
            def __init__(self, weight):
                super().__init__()
                self.weight = torch.nn.Parameter(torch.tensor(weight, dtype=torch.float64))

            def forward(self, x):  # (1, C_in) -> (1, n_modes)
                return x @ self.weight.t()

        torch.jit.script(ToModes(weight_rows)).save(str(self.checkpoint))

    def _CreateProcess(self, extra="", card_policy="advisory"):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import rom_surrogate_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"  : "Main",
                "rom_basis_folder" : "%s",
                "model_settings"   : { "checkpoint_file" : "%s", "device" : "cpu",
                                       "model_card_policy" : "%s" },
                "input_fields"     : [ { "variable_name" : "NODAL_PAUX", "data_location" : "node_historical" } ]
                %s
            }
        }""" % (self.basis_folder, self.checkpoint, card_policy, extra))
        return rom_surrogate_process.Factory(settings, self.model)

    def test_SurrogateReconstruction(self):
        model_part = self._CreateModelPart([Kratos.TEMPERATURE])
        phi = _OrthonormalBasis(4, 2)
        _WriteNumpyBasis(self.basis_folder, phi, [1, 2, 3, 4], ["TEMPERATURE"])
        weight = [[1.5], [-0.5]]  # q = [1.5 k, -0.5 k] with carrier k = 2
        self._SaveLinearModel(weight)

        process = self._CreateProcess()
        model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        expected_u = phi @ numpy.array([3.0, -1.0])
        for row, node in enumerate(model_part.Nodes):
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE), expected_u[row], places=10)

    def test_PermutationIsCachedAcrossStepsAndCorrectOnPermutedIds(self):
        # the basis lists the nodes in a different order than the part holds
        # them; every step must still land on the right node, and the
        # permutation (topology) must be built once, not per step
        model_part = self._CreateModelPart([Kratos.TEMPERATURE])
        _WriteNumpyBasis(self.basis_folder, numpy.eye(4)[:, :2], [3, 1, 4, 2], ["TEMPERATURE"])
        self._SaveLinearModel([[1.0], [2.0]])  # q = (2, 4) for the parameter 2
        process = self._CreateProcess()
        for step in (1, 2):
            model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
            if step == 1:
                permutation = process._permutation
            else:
                self.assertIs(process._permutation, permutation)
        # u = phi q = (2, 4, 0, 0) in basis row order [3, 1, 4, 2]
        self.assertAlmostEqual(model_part.GetNode(3).GetSolutionStepValue(Kratos.TEMPERATURE), 2.0)
        self.assertAlmostEqual(model_part.GetNode(1).GetSolutionStepValue(Kratos.TEMPERATURE), 4.0)
        self.assertAlmostEqual(model_part.GetNode(4).GetSolutionStepValue(Kratos.TEMPERATURE), 0.0)
        self.assertAlmostEqual(model_part.GetNode(2).GetSolutionStepValue(Kratos.TEMPERATURE), 0.0)

    def test_ComponentUnknowns(self):
        model_part = self._CreateModelPart([Kratos.DISPLACEMENT])
        phi = _OrthonormalBasis(8, 3)
        _WriteNumpyBasis(self.basis_folder, phi, [1, 2, 3, 4],
                         ["DISPLACEMENT_X", "DISPLACEMENT_Y"])
        self._SaveLinearModel([[1.0], [0.5], [-1.0]])

        process = self._CreateProcess()
        model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        expected_u = (phi @ numpy.array([2.0, 1.0, -2.0])).reshape(4, 2)
        for row, node in enumerate(model_part.Nodes):
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.DISPLACEMENT_X), expected_u[row, 0], places=10)
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y), expected_u[row, 1], places=10)

    def test_WrongModelOutputShapeRaises(self):
        model_part = self._CreateModelPart([Kratos.TEMPERATURE])
        _WriteNumpyBasis(self.basis_folder, _OrthonormalBasis(4, 2), [1, 2, 3, 4], ["TEMPERATURE"])
        self._SaveLinearModel([[1.0], [1.0], [1.0]])  # 3 modes, basis has 2

        process = self._CreateProcess()
        model_part.ProcessInfo[Kratos.STEP] = 1
        with self.assertRaisesRegex(ValueError, r"\(1, 2\) reduced-coordinates"):
            process.ExecuteFinalizeSolutionStep()

    def test_UserOutputFieldsRejected(self):
        self._CreateModelPart([Kratos.TEMPERATURE])
        _WriteNumpyBasis(self.basis_folder, _OrthonormalBasis(4, 2), [1, 2, 3, 4], ["TEMPERATURE"])
        self._SaveLinearModel([[1.0], [1.0]])
        with self.assertRaisesRegex(ValueError, "output_fields"):
            self._CreateProcess(extra=', "output_fields": [ { "variable_name": "TEMPERATURE" } ]')

    def test_BadInputReductionRaises(self):
        self._CreateModelPart([Kratos.TEMPERATURE])
        _WriteNumpyBasis(self.basis_folder, _OrthonormalBasis(4, 2), [1, 2, 3, 4], ["TEMPERATURE"])
        self._SaveLinearModel([[1.0], [1.0]])
        with self.assertRaisesRegex(ValueError, "input reduction"):
            self._CreateProcess(extra=', "input_reduction": "max"')

    def test_ModelCardValidatesAgainstNodalUnknowns(self):
        model_part = self._CreateModelPart([Kratos.TEMPERATURE])
        _WriteNumpyBasis(self.basis_folder, _OrthonormalBasis(4, 2), [1, 2, 3, 4], ["TEMPERATURE"])
        self._SaveLinearModel([[1.0], [1.0]])
        model_registry.SaveModelCard(self.checkpoint, {
            "input_fields": [{"variable_name": "NODAL_PAUX", "data_location": "node_historical"}],
            "output_fields": [{"variable_name": "TEMPERATURE", "data_location": "node_historical"}],
        })
        process = self._CreateProcess(card_policy="strict")
        model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()  # strict passes: card == derived specs
        self.assertNotEqual(
            model_part.GetNode(1).GetSolutionStepValue(Kratos.TEMPERATURE), 0.0)


if __name__ == '__main__':
    KratosUnittest.main()
