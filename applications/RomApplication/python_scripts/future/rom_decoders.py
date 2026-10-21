import numpy


class RomDecoder:
    """Base class of the decoders u = D(q) used by the Future ROM solver.

    The rows of the decoded vector and of the Jacobian follow the effective DOF set
    of the Future scheme (row i belongs to the DOF with effective equation id i).
    A decoder may treat the fixed DOFs in its own way (e.g. FixedDofsDecoder), reading
    the values imposed in the current step in InitializeSolutionStep.
    """

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
