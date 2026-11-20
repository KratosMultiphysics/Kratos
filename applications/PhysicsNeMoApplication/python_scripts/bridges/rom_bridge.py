"""Bridge to RomApplication's POD basis format.

Consumes the numpy-format output of RomApplication's
``CalculateRomBasisOutputProcess``: ``RightBasisMatrix.npy`` (the POD basis
``phi``), ``NodeIds.npy``, an optional ``SingularValuesVector.npy`` and the
``RomParameters.json`` sidecar. Only the FILE FORMAT is consumed - this
module needs numpy and the Kratos core, never the compiled RomApplication
(which is only required to *produce* a basis).

Row-ordering contract (mirrors the producer's snapshot assembly,
``numpy.stack([GetSolutionStepValuesVector(nodes, var, 0) for var in
sorted_vars], axis=1).reshape(-1, 1)``):

    row r of phi  <->  node  ``r // n_unknowns``  in ``NodeIds.npy`` order,
                       unknown ``r % n_unknowns`` in the (alphabetically
                       sorted) ``rom_settings.nodal_unknowns`` order.

Node-major, unknown-minor. ``nodal_unknowns`` entries are plain double
variables - scalars (``TEMPERATURE``) or components (``DISPLACEMENT_X``);
gather/scatter go through ``VariableUtils.Get/SetSolutionStepValuesVector``,
the very calls the producer uses, so orderings match by construction.

Pure numpy + Kratos: this module never imports torch or physicsnemo.
"""

import dataclasses
import json
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.utilities import array_backend_utils

# Measured crossover for the dense basis GEMV. Below ~1e7 basis entries the
# two backends are within noise of each other (the transfer and the kernel
# launch eat the matmul saving); at 1e7 the GPU is ~6x and at 2e7 ~9x.
_GPU_BASIS_THRESHOLD = 10000000


@dataclasses.dataclass(frozen=True)
class RomBasis:
    """A loaded POD basis (see the module docstring for the row ordering).

    Attributes:
        phi: (n_dofs, n_modes) float64 basis matrix (orthonormal columns as
            produced by the SVD).
        node_ids: (n_nodes,) int64, the producer's node ordering.
        nodal_unknowns: variable names in the stored (sorted) order.
        singular_values: (n_modes,) float64 or None.
    """
    phi: numpy.ndarray
    node_ids: numpy.ndarray
    nodal_unknowns: tuple
    singular_values: numpy.ndarray
    _device_phi: object = dataclasses.field(default=None, compare=False, repr=False)

    def _DevicePhi(self, xp):
        """phi in xp's memory, uploaded at most once.

        The basis is the large operand and it never changes between steps,
        so keeping it resident is what makes the GPU projection worth doing:
        only the (small) vector crosses the bus per call. The dataclass is
        frozen, hence object.__setattr__ for the cache slot.
        """
        if xp is numpy:
            return self.phi
        if self._device_phi is None:
            object.__setattr__(self, "_device_phi", xp.asarray(self.phi))
        return self._device_phi

    @property
    def n_nodes(self) -> int:
        return len(self.node_ids)

    @property
    def n_unknowns(self) -> int:
        return len(self.nodal_unknowns)

    @property
    def n_modes(self) -> int:
        return self.phi.shape[1]

    @property
    def n_dofs(self) -> int:
        return self.phi.shape[0]


def LoadRomBasis(folder, rom_parameters_name: str = "RomParameters") -> RomBasis:
    """Loads a numpy-format basis written by CalculateRomBasisOutputProcess.

    Args:
        folder: The rom_basis_output_folder of the producing process.
        rom_parameters_name: The rom_basis_output_name (json stem).

    Returns:
        The RomBasis. Raises with actionable messages for the json-format
        variant, missing files, and shape inconsistencies.
    """
    folder = Path(folder)
    parameters_file = folder / f"{rom_parameters_name}.json"
    if not parameters_file.is_file():
        raise FileNotFoundError(
            f"ROM parameters file \"{parameters_file}\" not found - point rom_basis_folder "
            "at a CalculateRomBasisOutputProcess output folder.")
    with open(parameters_file) as f:
        rom_parameters = json.load(f)

    rom_format = rom_parameters.get("rom_format", "numpy")
    if rom_format != "numpy":
        raise ValueError(
            f"Unsupported rom_format \"{rom_format}\": only the numpy layout is consumed. "
            "Regenerate the basis with \"rom_basis_output_format\": \"numpy\".")

    for name in ("RightBasisMatrix.npy", "NodeIds.npy"):
        if not (folder / name).is_file():
            raise FileNotFoundError(f"Basis file \"{folder / name}\" not found.")
    phi = numpy.load(folder / "RightBasisMatrix.npy").astype(numpy.float64)
    node_ids = numpy.load(folder / "NodeIds.npy").astype(numpy.int64)
    singular_values_file = folder / "SingularValuesVector.npy"
    singular_values = (numpy.load(singular_values_file).astype(numpy.float64)
                       if singular_values_file.is_file() else None)

    rom_settings = rom_parameters.get("rom_settings", {})
    nodal_unknowns = tuple(rom_settings.get("nodal_unknowns", []))
    if not nodal_unknowns:
        raise ValueError(f"\"{parameters_file}\" has no rom_settings.nodal_unknowns.")
    if phi.ndim != 2 or phi.shape[0] != len(node_ids) * len(nodal_unknowns):
        raise ValueError(
            f"RightBasisMatrix has shape {list(phi.shape)} but {len(node_ids)} node(s) x "
            f"{len(nodal_unknowns)} unknown(s) imply {len(node_ids) * len(nodal_unknowns)} rows.")
    number_of_modes = rom_settings.get("number_of_rom_dofs")
    if number_of_modes is not None and phi.shape[1] != number_of_modes:
        raise ValueError(
            f"RightBasisMatrix has {phi.shape[1]} mode(s) but rom_settings.number_of_rom_dofs "
            f"says {number_of_modes}.")

    return RomBasis(phi=phi, node_ids=node_ids, nodal_unknowns=nodal_unknowns,
                    singular_values=singular_values)


def NodePermutation(model_part: Kratos.ModelPart, rom_basis: RomBasis) -> numpy.ndarray:
    """Positions of the basis's node ids inside the model part's Nodes order.

    Topology: a caller stepping a static mesh (RomSurrogateProcess) builds
    this once and hands it to Gather/ScatterUnknownsVector as
    ``permutation=`` instead of paying the O(N) lookup every step. The
    lookup is the shared ``RowsOfIds`` (argsort + searchsorted), not a
    per-call ``{id: row}`` dict plus a generator.
    """
    from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import (
        RowsOfIds)
    part_ids = numpy.fromiter((node.Id for node in model_part.Nodes),
                              dtype=numpy.int64, count=model_part.NumberOfNodes())
    try:
        return RowsOfIds(part_ids, rom_basis.node_ids)
    except KeyError as error:
        raise RuntimeError(
            f"Basis node id {error.args[0]} is not in model part "
            f"\"{model_part.FullName()}\" - the basis belongs to a different mesh.") from None


_NodePermutation = NodePermutation  # the pre-public name


def GatherUnknownsVector(model_part: Kratos.ModelPart, rom_basis: RomBasis,
                         permutation=None) -> numpy.ndarray:
    """Reads the basis's nodal unknowns into a (n_dofs,) vector in the exact
    basis row order (node-major in NodeIds order, unknown-minor).

    Args:
        permutation: Optional precomputed NodePermutation result.
    """
    if permutation is None:
        permutation = NodePermutation(model_part, rom_basis)
    columns = []
    for name in rom_basis.nodal_unknowns:
        variable = Kratos.KratosGlobals.GetVariable(name)
        values = numpy.array(
            Kratos.VariableUtils().GetSolutionStepValuesVector(model_part.Nodes, variable, 0))
        columns.append(values[permutation])
    return numpy.stack(columns, axis=1).reshape(-1)


def ScatterUnknownsVector(model_part: Kratos.ModelPart, rom_basis: RomBasis, u,
                          permutation=None) -> None:
    """Writes a (n_dofs,) vector in basis row order back onto the nodes.

    Args:
        permutation: Optional precomputed NodePermutation result.
    """
    u = numpy.asarray(u, dtype=numpy.float64)
    if u.shape != (rom_basis.n_dofs,):
        raise ValueError(
            f"Expected a ({rom_basis.n_dofs},) unknowns vector, got shape {list(u.shape)}.")
    if permutation is None:
        permutation = NodePermutation(model_part, rom_basis)
    inverse = numpy.empty_like(permutation)
    inverse[permutation] = numpy.arange(len(permutation))
    per_node = u.reshape(rom_basis.n_nodes, rom_basis.n_unknowns)
    for j, name in enumerate(rom_basis.nodal_unknowns):
        variable = Kratos.KratosGlobals.GetVariable(name)
        # Kratos.Vector reads the raw buffer and IGNORES numpy strides -
        # always hand it a contiguous array.
        values = numpy.ascontiguousarray(per_node[inverse, j])
        Kratos.VariableUtils().SetSolutionStepValuesVector(
            model_part.Nodes, variable, Kratos.Vector(values), 0)


def ProjectToReducedSpace(rom_basis: RomBasis, u, backend="numpy") -> numpy.ndarray:
    """q = phi^T u - the L2-optimal reduced coordinates for orthonormal-column
    bases (which the producing SVD guarantees). Accepts a (n_dofs,) vector or
    a (n_dofs, T) snapshot series.

    backend: "numpy" (default), "cupy" or "auto". The crossover measured on
    this dense GEMV sits near phi.size = 1e7: below it the two are within
    noise (a 46875x50 basis measures ~1.0x), at 100000x100 the GPU is ~6.3x
    and at 200000x100 ~8.9x. Below _GPU_BASIS_THRESHOLD a CuPy request falls
    back to numpy on its own.
    """
    u = numpy.asarray(u, dtype=numpy.float64)
    if u.shape[0] != rom_basis.n_dofs:
        raise ValueError(
            f"Expected {rom_basis.n_dofs} row(s) (the basis DOFs), got {u.shape[0]}.")
    xp, _ = array_backend_utils.ResolveArrayModule(
        backend, size_hint=rom_basis.phi.size, threshold=_GPU_BASIS_THRESHOLD)
    return array_backend_utils.ToHost(rom_basis._DevicePhi(xp).T @ xp.asarray(u))


def ReconstructFromReducedSpace(rom_basis: RomBasis, q, backend="numpy") -> numpy.ndarray:
    """u = phi q. Accepts (n_modes,) or (n_modes, T).

    backend: as in ProjectToReducedSpace.
    """
    q = numpy.asarray(q, dtype=numpy.float64)
    if q.shape[0] != rom_basis.n_modes:
        raise ValueError(
            f"Expected {rom_basis.n_modes} reduced coordinate(s), got {q.shape[0]}.")
    xp, _ = array_backend_utils.ResolveArrayModule(
        backend, size_hint=rom_basis.phi.size, threshold=_GPU_BASIS_THRESHOLD)
    return array_backend_utils.ToHost(rom_basis._DevicePhi(xp) @ xp.asarray(q))
