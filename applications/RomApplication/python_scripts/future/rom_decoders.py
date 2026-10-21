import pathlib

import numpy


class RomDecoder:
    """Base class of the decoders u = D(q) used by the Future ROM solver.

    The rows of the decoded vector and of the Jacobian follow the effective DOF set
    of the Future scheme (row i belongs to the DOF with effective equation id i).
    """

    def NumberOfRomDofs(self):
        raise NotImplementedError("Calling base class 'NumberOfRomDofs'.")

    def Decode(self, q):
        """Returns the full order solution vector D(q)."""
        raise NotImplementedError("Calling base class 'Decode'.")

    def Jacobian(self, q):
        """Returns the dense matrix dD/dq evaluated at q."""
        raise NotImplementedError("Calling base class 'Jacobian'.")


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

    def NumberOfRomDofs(self):
        return self.phi.shape[1]

    def Decode(self, q):
        return self.u0 + self.phi @ q

    def Jacobian(self, q):
        return self.phi


class TorchDecoder(RomDecoder):
    """Decoder u = g(q) given by a function of PyTorch tensors (e.g. a neural network). Its Jacobian is obtained by automatic
    differentiation (forward mode), so only g is implemented. Requires PyTorch, which is imported here only.

    The rows of g(q) follow the effective DOF set (see 'decoder_factory' in CreateFutureRomAnalysisInstance).
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


class AnnEnhancedDecoder(TorchDecoder):
    """ANN-enhanced decoder of the legacy 'galerkin_ann' ROM (Barnett, Farhat and Maday 2023), with the same data:
    u = Phi_inf q + Phi_sup Sigma_sup N(Sigma_inf^-1 q), with N a network of dense layers without bias, ELU activation in all
    of them but the last (h <- W^T h). The Jacobian is obtained by automatic differentiation (see TorchDecoder).
    phi: basis with at least n_sup columns, sigma: the singular values used to scale it in the training, layers: the weight
    matrices of the network, modes: [n_inf, n_sup]."""

    def __init__(self, phi, sigma, layers, modes):
        import torch
        n_inf, n_sup = int(modes[0]), int(modes[1])
        phi = torch.as_tensor(numpy.asarray(phi, dtype=numpy.float64))
        sigma = torch.as_tensor(numpy.asarray(sigma, dtype=numpy.float64))
        phi_inf = phi[:, :n_inf]
        phi_sup_weighted = phi[:, n_inf:n_sup] * sigma[n_inf:n_sup]
        sigma_inf = sigma[:n_inf]
        layers = [torch.as_tensor(numpy.asarray(layer, dtype=numpy.float64)) for layer in layers]

        def Function(q):
            h = q / sigma_inf
            for i, layer in enumerate(layers):
                h = layer.T @ h
                if i < len(layers) - 1:
                    h = torch.nn.functional.elu(h)
            return phi_inf @ q + phi_sup_weighted @ h

        super().__init__(Function, n_inf)

    @classmethod
    def FromNumpyFiles(cls, rom_data_folder, modes, rows):
        """From the files of a standalone legacy ANN-enhanced ROM (RightBasisMatrix.npy, SingularValues.npy and model_weights.npy,
        see NN_ROM_Interface.FromNumpyFiles), with the rows of the basis in the order of the effective DOF set"""
        folder = pathlib.Path(rom_data_folder)
        return cls(numpy.load(folder / "RightBasisMatrix.npy")[rows, :], numpy.load(folder / "SingularValues.npy"),
                   numpy.load(folder / "model_weights.npy", allow_pickle=True), modes)
