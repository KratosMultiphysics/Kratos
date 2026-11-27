"""Tests for rom_bridge: RomApplication basis-format consumption (ungated -
pure numpy + Kratos, no compiled RomApplication needed)."""

import json
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.bridges import rom_bridge
def _WriteNumpyBasis(folder, phi, node_ids, nodal_unknowns, singular_values=None,
                     rom_format="numpy", name="RomParameters"):
    """Writes the exact CalculateRomBasisOutputProcess numpy file layout."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    rom_parameters = {
        "rom_format": rom_format,
        "rom_settings": {
            "nodal_unknowns": list(nodal_unknowns),
            "number_of_rom_dofs": int(phi.shape[1]),
        },
    }
    with open(folder / f"{name}.json", "w") as f:
        json.dump(rom_parameters, f, indent=4)
    if rom_format == "numpy":
        numpy.save(folder / "RightBasisMatrix.npy", phi)
        numpy.save(folder / "NodeIds.npy", numpy.asarray(node_ids))
        if singular_values is not None:
            numpy.save(folder / "SingularValuesVector.npy", singular_values)


def _OrthonormalBasis(n_dofs, n_modes, seed=0):
    rng = numpy.random.default_rng(seed)
    q, _ = numpy.linalg.qr(rng.standard_normal((n_dofs, n_modes)))
    return q


class TestRomBridge(KratosUnittest.TestCase):
    def setUp(self):
        self.folder = Path("test_rom_bridge_basis")
        self.model = Kratos.Model()

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.folder))

    def _CreateNodes(self, variables, node_ids=(1, 2, 3, 4)):
        model_part = self.model.CreateModelPart("Main")
        for variable in variables:
            model_part.AddNodalSolutionStepVariable(variable)
        for i, node_id in enumerate(node_ids):
            model_part.CreateNewNode(int(node_id), float(i), 0.0, 0.0)
        return model_part

    def test_LoadRomBasisNumpyFormat(self):
        phi = _OrthonormalBasis(8, 3)
        sigma = numpy.array([3.0, 2.0, 1.0])
        _WriteNumpyBasis(self.folder, phi, [1, 2, 3, 4], ["DISPLACEMENT_X", "DISPLACEMENT_Y"], sigma)

        basis = rom_bridge.LoadRomBasis(self.folder)
        self.assertEqual(basis.n_nodes, 4)
        self.assertEqual(basis.n_unknowns, 2)
        self.assertEqual(basis.n_modes, 3)
        self.assertEqual(basis.n_dofs, 8)
        self.assertEqual(basis.nodal_unknowns, ("DISPLACEMENT_X", "DISPLACEMENT_Y"))
        numpy.testing.assert_allclose(basis.singular_values, sigma)

        # singular values are optional
        KratosUtilities.DeleteFileIfExisting(str(self.folder / "SingularValuesVector.npy"))
        self.assertIsNone(rom_bridge.LoadRomBasis(self.folder).singular_values)

    def test_JsonFormatBasisRaises(self):
        _WriteNumpyBasis(self.folder, _OrthonormalBasis(4, 2), [1, 2, 3, 4],
                         ["TEMPERATURE"], rom_format="json")
        with self.assertRaisesRegex(ValueError, "numpy"):
            rom_bridge.LoadRomBasis(self.folder)

    def test_MissingFilesRaise(self):
        with self.assertRaisesRegex(FileNotFoundError, "RomParameters"):
            rom_bridge.LoadRomBasis(self.folder)
        _WriteNumpyBasis(self.folder, _OrthonormalBasis(4, 2), [1, 2, 3, 4], ["TEMPERATURE"])
        KratosUtilities.DeleteFileIfExisting(str(self.folder / "NodeIds.npy"))
        with self.assertRaisesRegex(FileNotFoundError, "NodeIds"):
            rom_bridge.LoadRomBasis(self.folder)

    def test_ShapeInconsistenciesRaise(self):
        _WriteNumpyBasis(self.folder, _OrthonormalBasis(6, 2), [1, 2, 3, 4], ["TEMPERATURE"])
        with self.assertRaisesRegex(ValueError, "imply 4 rows"):
            rom_bridge.LoadRomBasis(self.folder)

    def test_GatherScatterRoundTripScalar(self):
        model_part = self._CreateNodes([Kratos.TEMPERATURE])
        _WriteNumpyBasis(self.folder, _OrthonormalBasis(4, 2), [1, 2, 3, 4], ["TEMPERATURE"])
        basis = rom_bridge.LoadRomBasis(self.folder)

        u = numpy.array([10.0, 20.0, 30.0, 40.0])
        rom_bridge.ScatterUnknownsVector(model_part, basis, u)
        for node in model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE), 10.0 * node.Id, places=12)
        numpy.testing.assert_allclose(rom_bridge.GatherUnknownsVector(model_part, basis), u)

    def test_GatherScatterInterleavingComponents(self):
        model_part = self._CreateNodes([Kratos.DISPLACEMENT])
        unknowns = ["DISPLACEMENT_X", "DISPLACEMENT_Y"]
        _WriteNumpyBasis(self.folder, _OrthonormalBasis(8, 2), [1, 2, 3, 4], unknowns)
        basis = rom_bridge.LoadRomBasis(self.folder)

        # value = 10*id + component -> checks node-major, unknown-minor exactly
        for node in model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [10.0 * node.Id, 10.0 * node.Id + 1.0, 0.0])
        u = rom_bridge.GatherUnknownsVector(model_part, basis)
        numpy.testing.assert_allclose(
            u, [10.0, 11.0, 20.0, 21.0, 30.0, 31.0, 40.0, 41.0])

        rom_bridge.ScatterUnknownsVector(model_part, basis, u + 100.0)
        node = model_part.GetNode(3)
        self.assertAlmostEqual(node.GetSolutionStepValue(Kratos.DISPLACEMENT_X), 130.0, places=12)
        self.assertAlmostEqual(node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y), 131.0, places=12)

    def test_ProjectReconstructRoundTrip(self):
        phi = _OrthonormalBasis(12, 4)
        _WriteNumpyBasis(self.folder, phi, list(range(1, 13)), ["TEMPERATURE"])
        basis = rom_bridge.LoadRomBasis(self.folder)

        q_true = numpy.array([1.0, -2.0, 0.5, 3.0])
        u = phi @ q_true  # in the span -> exact round trip
        numpy.testing.assert_allclose(rom_bridge.ProjectToReducedSpace(basis, u), q_true, atol=1e-12)
        numpy.testing.assert_allclose(
            rom_bridge.ReconstructFromReducedSpace(basis, q_true), u, atol=1e-12)

        # snapshot series (n_dofs, T)
        series = phi @ numpy.stack([q_true, 2.0 * q_true], axis=1)
        q_series = rom_bridge.ProjectToReducedSpace(basis, series)
        self.assertEqual(q_series.shape, (4, 2))
        numpy.testing.assert_allclose(
            rom_bridge.ReconstructFromReducedSpace(basis, q_series), series, atol=1e-12)

    def test_PermutedNodeIdsRespected(self):
        model_part = self._CreateNodes([Kratos.TEMPERATURE])
        _WriteNumpyBasis(self.folder, _OrthonormalBasis(4, 2), [3, 1, 4, 2], ["TEMPERATURE"])
        basis = rom_bridge.LoadRomBasis(self.folder)

        rom_bridge.ScatterUnknownsVector(model_part, basis, numpy.array([1.0, 2.0, 3.0, 4.0]))
        # basis row order is [3, 1, 4, 2]
        self.assertAlmostEqual(model_part.GetNode(3).GetSolutionStepValue(Kratos.TEMPERATURE), 1.0)
        self.assertAlmostEqual(model_part.GetNode(1).GetSolutionStepValue(Kratos.TEMPERATURE), 2.0)
        self.assertAlmostEqual(model_part.GetNode(4).GetSolutionStepValue(Kratos.TEMPERATURE), 3.0)
        self.assertAlmostEqual(model_part.GetNode(2).GetSolutionStepValue(Kratos.TEMPERATURE), 4.0)
        numpy.testing.assert_allclose(
            rom_bridge.GatherUnknownsVector(model_part, basis), [1.0, 2.0, 3.0, 4.0])

    def test_PrecomputedPermutationIsHonoured(self):
        model_part = self._CreateNodes([Kratos.TEMPERATURE])
        _WriteNumpyBasis(self.folder, _OrthonormalBasis(4, 2), [3, 1, 4, 2], ["TEMPERATURE"])
        basis = rom_bridge.LoadRomBasis(self.folder)
        permutation = rom_bridge.NodePermutation(model_part, basis)
        numpy.testing.assert_array_equal(permutation, [2, 0, 3, 1])
        rom_bridge.ScatterUnknownsVector(model_part, basis, numpy.array([1.0, 2.0, 3.0, 4.0]),
                                         permutation=permutation)
        numpy.testing.assert_allclose(
            rom_bridge.GatherUnknownsVector(model_part, basis, permutation=permutation),
            [1.0, 2.0, 3.0, 4.0])
        # a wrong permutation is honoured too - it is the caller's contract
        rom_bridge.ScatterUnknownsVector(model_part, basis, numpy.array([1.0, 2.0, 3.0, 4.0]),
                                         permutation=numpy.array([0, 1, 2, 3]))
        self.assertAlmostEqual(model_part.GetNode(1).GetSolutionStepValue(Kratos.TEMPERATURE), 1.0)

    def test_MissingNodeIdRaises(self):
        model_part = self._CreateNodes([Kratos.TEMPERATURE], node_ids=(1, 2, 3))
        _WriteNumpyBasis(self.folder, _OrthonormalBasis(4, 2), [1, 2, 3, 4], ["TEMPERATURE"])
        basis = rom_bridge.LoadRomBasis(self.folder)
        with self.assertRaisesRegex(RuntimeError, "different mesh"):
            rom_bridge.GatherUnknownsVector(model_part, basis)


class TestRowsOfIds(KratosUnittest.TestCase):
    """The shared id -> row lookup (graph scatter, ROM permutation, grid
    sampling). Pure numpy."""

    def test_UnsortedContainerAndShapePreserved(self):
        from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import RowsOfIds
        container = numpy.array([30, 10, 40, 20])
        rows = RowsOfIds(container, [[20, 30], [10, 40]])
        numpy.testing.assert_array_equal(rows, [[3, 0], [1, 2]])
        numpy.testing.assert_array_equal(container[rows], [[20, 30], [10, 40]])
        self.assertEqual(RowsOfIds(container, []).shape, (0,))
        self.assertEqual(RowsOfIds([], []).shape, (0,))

    def test_MissingIdIsNamed(self):
        from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import RowsOfIds
        with self.assertRaises(KeyError) as context:
            RowsOfIds([1, 2, 3], [2, 7])
        self.assertEqual(context.exception.args[0], 7)
        with self.assertRaises(KeyError):
            RowsOfIds([], [1])


if __name__ == '__main__':
    KratosUnittest.main()
