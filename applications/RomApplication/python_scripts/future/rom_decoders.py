import json
import pathlib

import numpy

# Decoders that CreateDecoder builds from a ROM folder, by their 'type_of_decoder' (see RegisterDecoder)
_DECODERS = {}

# Gauss-Newton iterations of RomDecoder.Encode: stop when |dq| <= tolerance |q|
_ENCODE_MAX_ITERATIONS = 50
_ENCODE_RELATIVE_TOLERANCE = 1.0e-12


def RegisterDecoder(type_of_decoder, decoder_class):
    """Makes CreateDecoder build the ROM folders of the given 'type_of_decoder' with decoder_class.FromRomData."""
    _DECODERS[type_of_decoder] = decoder_class


def GetRomParameters(rom_data_folder):
    """The RomParameters.json of a ROM folder, as saved by the RomManager."""
    with open(pathlib.Path(rom_data_folder) / "RomParameters.json", "r") as parameter_file:
        return json.load(parameter_file)


def GetTypeOfDecoder(rom_parameters):
    """The 'type_of_decoder' of the RomParameters.json. If not given, the one of the legacy 'projection_strategy' (e.g. "galerkin_ann"
    for "ann_enhanced" and "galerkin_rbf" for "rbf_enhanced", as written by the RomManager), else "linear"."""
    if "type_of_decoder" in rom_parameters:
        return rom_parameters["type_of_decoder"]
    projection_strategy = rom_parameters.get("projection_strategy", "")
    for suffix, type_of_decoder in (("_ann", "ann_enhanced"), ("_rbf", "rbf_enhanced")):
        if projection_strategy.endswith(suffix):
            return type_of_decoder
    return "linear"


def GetBasisRows(rom_data_folder, dof_set):
    """The bases saved in a ROM folder have one row per node and nodal unknown (NodeIds.npy and the 'nodal_unknowns' of the
    RomParameters.json). Returns their rows in the order of the DOFs of dof_set (e.g. the effective DOF set of the Future scheme)."""
    folder = pathlib.Path(rom_data_folder)
    nodal_unknowns = GetRomParameters(folder)["rom_settings"]["nodal_unknowns"]
    node_id_to_index = {node_id: index for index, node_id in enumerate(numpy.load(folder / "NodeIds.npy"))}
    unknown_to_index = {name: index for index, name in enumerate(nodal_unknowns)}
    return [node_id_to_index[dof.Id()] * len(nodal_unknowns) + unknown_to_index[dof.GetVariable().Name()] for dof in dof_set]


def CreateDecoder(rom_data_folder, dof_set):
    """The decoder saved in a ROM folder, with its rows in the order of the DOFs of dof_set: the class registered for the type of
    decoder of its RomParameters.json (see GetTypeOfDecoder), built with its FromRomData."""
    type_of_decoder = GetTypeOfDecoder(GetRomParameters(rom_data_folder))
    if type_of_decoder not in _DECODERS:
        err_msg = "Unknown 'type_of_decoder' \'{}\' in {}. Available options are: {}.".format(type_of_decoder, rom_data_folder, ", ".join("\'{}\'".format(name) for name in _DECODERS))
        raise Exception(err_msg)
    return _DECODERS[type_of_decoder].FromRomData(rom_data_folder, dof_set)


class RomDecoder:
    """Base class of the decoders u = D(q) used by the Future ROM solver.

    The rows of the decoded vector and of the Jacobian follow the effective DOF set
    of the Future scheme (row i belongs to the DOF with effective equation id i).
    A decoder may treat the fixed DOFs in its own way (e.g. FixedDofsDecoder), reading
    the values imposed in the current step in InitializeSolutionStep.
    A decoder that can be saved in a ROM folder implements FromRomData and is registered (RegisterDecoder),
    so that CreateDecoder builds it from the folder alone.
    """

    @classmethod
    def FromRomData(cls, rom_data_folder, dof_set):
        """Creates the decoder from the files of a ROM folder, with its rows in the order of the DOFs of dof_set."""
        raise NotImplementedError("Calling base class 'FromRomData'.")

    def NumberOfRomDofs(self):
        raise NotImplementedError("Calling base class 'NumberOfRomDofs'.")

    def Decode(self, q):
        """Returns the full order solution vector D(q)."""
        raise NotImplementedError("Calling base class 'Decode'.")

    def Jacobian(self, q):
        """Returns the dense matrix dD/dq evaluated at q."""
        raise NotImplementedError("Calling base class 'Jacobian'.")

    def Encode(self, u, q=None):
        """Returns the reduced coordinates whose decoding is closest to u (least squares, min |D(q) - u|): q itself for u = D(q).
        Gauss-Newton iterations with Decode and Jacobian from q (e.g. the reduced coordinates of the previous step, else zero).
        Decoders with a cheaper or exact inverse override it."""
        q = numpy.zeros(self.NumberOfRomDofs()) if q is None else numpy.array(q, dtype=numpy.float64)
        for _ in range(_ENCODE_MAX_ITERATIONS):
            dq = numpy.linalg.lstsq(self.Jacobian(q), u - self.Decode(q), rcond=None)[0]
            q += dq
            if numpy.linalg.norm(dq) <= _ENCODE_RELATIVE_TOLERANCE * numpy.linalg.norm(q):
                return q
        raise Exception("Encode did not converge in {} Gauss-Newton iterations.".format(_ENCODE_MAX_ITERATIONS))

    def InitializeSolutionStep(self):
        """Called once the conditions of the step are applied (fixed DOFs and their values)."""
        pass


class LinearDecoder(RomDecoder):
    """Affine decoder u = u0 + Phi q (standard POD basis)."""

    def __init__(self, phi, u0=None):
        # C-ordered so that it can be passed to the C++ projector without copies
        self.phi = numpy.ascontiguousarray(phi, dtype=numpy.float64)
        if self.phi.ndim != 2:
            raise Exception("The basis 'phi' must be a matrix with one row per DOF and one column per mode.")
        self.u0 = numpy.zeros(self.phi.shape[0]) if u0 is None else numpy.array(u0, dtype=numpy.float64)
        if self.u0.shape != (self.phi.shape[0],):
            raise Exception("Size of 'u0' ({}) does not match the rows of 'phi' ({}).".format(self.u0.shape, self.phi.shape[0]))

    @classmethod
    def FromRomData(cls, rom_data_folder, dof_set):
        """From the basis saved by the RomManager (the first 'number_of_rom_dofs' columns of RightBasisMatrix.npy), taking the
        current values of the DOFs (the initial state) as origin."""
        folder = pathlib.Path(rom_data_folder)
        phi = numpy.load(folder / "RightBasisMatrix.npy")
        if phi.ndim == 1:
            phi = phi.reshape(-1, 1)
        phi = phi[:, :GetRomParameters(folder)["rom_settings"]["number_of_rom_dofs"]]
        return cls(phi[GetBasisRows(folder, dof_set), :], numpy.array(dof_set.GetValues()))

    def NumberOfRomDofs(self):
        return self.phi.shape[1]

    def Decode(self, q):
        return self.u0 + self.phi @ q

    def Jacobian(self, q):
        return self.phi

    def Encode(self, u, q=None):
        """Least squares on the basis, q = Phi^+ (u - u0), through its thin QR factorization (computed once)."""
        if not hasattr(self, "_phi_qr"):
            self._phi_qr = numpy.linalg.qr(self.phi)
        phi_q, phi_r = self._phi_qr
        return numpy.linalg.solve(phi_r, phi_q.T @ (u - self.u0))


RegisterDecoder("linear", LinearDecoder)


class FixedDofsDecoder(RomDecoder):
    """Decoder with two DOF sets: the free DOFs are given by the provided decoder and the fixed DOFs by their prescribed values.

    The fixed DOFs (e.g. Dirichlet conditions) are not approximated: their values are the ones imposed in the model, read in
    InitializeSolutionStep. Their rows of the Jacobian are zero. Any decoder can be wrapped (linear or nonlinear).
    Encode is the Gauss-Newton one of the base class: as the Jacobian is zero on the fixed DOFs, it fits the free DOFs only.
    """

    def __init__(self, decoder, dof_set):
        self.decoder = decoder
        self.dof_set = dof_set
        self.InitializeSolutionStep()

    def InitializeSolutionStep(self):
        """Reads which DOFs are fixed and their prescribed values"""
        self.decoder.InitializeSolutionStep()
        self.fixed_dofs = numpy.array([dof.IsFixed() for dof in self.dof_set], dtype=bool)
        self.fixed_values = numpy.array(self.dof_set.GetValues())[self.fixed_dofs]

    def NumberOfRomDofs(self):
        return self.decoder.NumberOfRomDofs()

    def Decode(self, q):
        u = numpy.array(self.decoder.Decode(q), dtype=numpy.float64)
        u[self.fixed_dofs] = self.fixed_values
        return u

    def Jacobian(self, q):
        # A new C-ordered array, so that it can be passed to the C++ projector without copies
        jacobian = numpy.array(self.decoder.Jacobian(q), dtype=numpy.float64, order="C")
        jacobian[self.fixed_dofs, :] = 0.0
        return jacobian
