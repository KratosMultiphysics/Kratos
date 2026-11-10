import json
import pathlib

import numpy

# Decoders that CreateDecoder builds from a ROM folder, by their 'type_of_decoder' (see RegisterDecoder)
_DECODERS = {}


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


RegisterDecoder("linear", LinearDecoder)


class TorchDecoder(RomDecoder):
    """Decoder u = g(q) given by a function of PyTorch tensors (e.g. a neural network). Its Jacobian is obtained by automatic
    differentiation (forward mode), so only g is implemented. Requires PyTorch, which is imported here only.

    The rows of g(q) follow the effective DOF set (see GetBasisRows, and 'decoder_factory' in CreateFutureRomAnalysisInstance).
    It is not registered for CreateDecoder, as it is given by a function, not by the files of a ROM folder.
    """

    def __init__(self, function, number_of_rom_dofs):
        import torch
        self.torch = torch
        self.function = function
        self.number_of_rom_dofs = number_of_rom_dofs
        self.jacobian_function = torch.func.jacfwd(function)

    def NumberOfRomDofs(self):
        return self.number_of_rom_dofs

    def Decode(self, q):
        with self.torch.no_grad():
            return self.function(self.torch.as_tensor(q, dtype=self.torch.float64)).numpy()

    def Jacobian(self, q):
        # C-ordered so that it can be passed to the C++ projector without copies
        return numpy.ascontiguousarray(self.jacobian_function(self.torch.as_tensor(q, dtype=self.torch.float64)).detach().numpy())


class RbfEnhancedDecoder(RomDecoder):
    """RBF-enhanced decoder of the legacy 'galerkin_rbf' ROM, with the same data:
    u = Phi_inf q + Phi_sup Sigma_sup W^T k(s), s = Sigma_inf^-1 q, with k_i the kernel at the distance from s to the center c_i,
    gaussian exp(-(eps r)^2) or inverse multiquadric 1 / sqrt(1 + (eps r)^2).
    phi: basis with at least n_sup columns, sigma: the singular values used to scale it in the training, weights: the RBF weights
    (number of centers x (n_sup - n_inf)), centers: the scaled inferior reduced coordinates of the training snapshots (one per row)."""

    def __init__(self, phi, sigma, weights, centers, kernel, epsilon):
        self.weights = numpy.asarray(weights, dtype=numpy.float64)
        self.centers = numpy.asarray(centers, dtype=numpy.float64)
        n_inf = self.centers.shape[1]
        n_sup = n_inf + self.weights.shape[1]
        sigma = numpy.asarray(sigma, dtype=numpy.float64)
        self.phi_inf = numpy.ascontiguousarray(phi[:, :n_inf], dtype=numpy.float64)
        self.phi_sup_weighted = phi[:, n_inf:n_sup] * sigma[n_inf:n_sup]
        self.sigma_inf = sigma[:n_inf]
        if kernel not in ("gaussian", "imq"):
            raise Exception("Unknown RBF kernel '{}'. Available options are 'gaussian' and 'imq'.".format(kernel))
        self.kernel = kernel
        self.epsilon = float(epsilon)

    def NumberOfRomDofs(self):
        return self.phi_inf.shape[1]

    def _Kernels(self, q):
        """Kernels, and the factors f_i such that the gradient of k_i with respect to s is f_i (s - c_i)"""
        differences = q / self.sigma_inf - self.centers
        squared_distances = self.epsilon**2 * numpy.sum(differences**2, axis=1)
        if self.kernel == "gaussian":
            kernels = numpy.exp(-squared_distances)
            factors = -2.0 * self.epsilon**2 * kernels
        else:
            kernels = 1.0 / numpy.sqrt(1.0 + squared_distances)
            factors = -self.epsilon**2 * kernels**3
        return kernels, factors, differences

    def Decode(self, q):
        kernels, _, _ = self._Kernels(q)
        return self.phi_inf @ q + self.phi_sup_weighted @ (self.weights.T @ kernels)

    def Jacobian(self, q):
        _, factors, differences = self._Kernels(q)
        rbf_gradient = self.weights.T @ (factors[:, None] * differences) / self.sigma_inf
        return numpy.ascontiguousarray(self.phi_inf + self.phi_sup_weighted @ rbf_gradient)

    @classmethod
    def FromNumpyFiles(cls, rom_data_folder, kernel, epsilon, rows):
        """From the files of a standalone legacy RBF-enhanced ROM (RightBasisMatrix.npy, SingularValues.npy, rbf_weights.npy and
        rbf_centers.npy, see RBF_ROM_Interface.FromNumpyFiles), with the rows of the basis in the order of the effective DOF set"""
        folder = pathlib.Path(rom_data_folder)
        return cls(numpy.load(folder / "RightBasisMatrix.npy")[rows, :], numpy.load(folder / "SingularValues.npy"),
                   numpy.load(folder / "rbf_weights.npy"), numpy.load(folder / "rbf_centers.npy"), kernel, epsilon)


class FixedDofsDecoder(RomDecoder):
    """Decoder with two DOF sets: the free DOFs are given by the provided decoder and the fixed DOFs by their prescribed values.

    The fixed DOFs (e.g. Dirichlet conditions) are not approximated: their values are the ones imposed in the model, read in
    InitializeSolutionStep. Their rows of the Jacobian are zero. Any decoder can be wrapped (linear or nonlinear).
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
