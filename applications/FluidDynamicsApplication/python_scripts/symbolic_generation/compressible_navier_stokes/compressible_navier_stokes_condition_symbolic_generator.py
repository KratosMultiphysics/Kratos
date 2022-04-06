import sympy

import KratosMultiphysics
import KratosMultiphysics.sympy_fe_utilities as KratosSympy

from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src import generate_convective_flux
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src import generate_diffusive_flux

from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src \
    .base_generator import CompressibleNavierStokesSymbolicGeneratorBase

from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src \
    .defines import CompressibleNavierStokesDefines as defs

from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src \
    .symbolic_parameters import FormulationParameters, ShockCapturingParameters, ShockCapturingNodalParameters


class CompressibleNavierStokesConditionSymbolicGenerator(CompressibleNavierStokesSymbolicGeneratorBase):
    """
    This class is in charge of generating a condition to solve the compressible
    Navier-Stokes using the conservative magnitudes (density, momentum, energy)
    as solution.

    Formulation
    -----------
    This formulation is a variant of:
    Bayona (2017). Adaptive Mesh Simulations of Compressible Flows using Stabilized
    Formulations. Chapter 3.
    """

    def __init__(self, settings):
        super().__init__(settings)
        self.is_neumann = settings["neumann"].GetBool()

    @classmethod
    def GetDefaultParameters(cls):
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "neumann": true
        }""")
        default_parameters.RecursivelyAddMissingParameters(super().GetDefaultParameters())
        return default_parameters

    def TouchFiles(self):
        super()._TouchFiles(__file__)

    def _ComputeFluxes(self, n, Ug, H, params, sc_params):
        KratosMultiphysics.Logger.Print(" - Compute Euler Jacobian matrix")
        E = generate_convective_flux.ComputeConvectiveFlux(Ug, params)

        if self.shock_capturing:
            sc_params = ShockCapturingParameters()
            KratosMultiphysics.Logger.Print(" - Compute diffusive flux (shock capturing ON)")
            G = generate_diffusive_flux.ComputeDiffusiveFluxWithShockCapturing(Ug, H, params, sc_params)
        else:
            sc_params = None
            KratosMultiphysics.Logger.Print(" - Compute diffusive flux (shock capturing OFF)")
            G = generate_diffusive_flux.ComputeDiffusiveFlux(Ug, H, params)

        return (E + G) * n

    def _ComputeVariationalFormulation(self, flux, V):
        KratosMultiphysics.Logger.Print(" - Compute variational formulation")
        return -V.transpose()*flux

    def _ComputeResidualAtGaussPointNeumann(self, F, F_gauss, i_gauss, rv_gauss, V, w):
        KratosMultiphysics.Logger.Print("    Gauss point: {}".format(i_gauss))

        # Get Gauss point geometry data
        Ng = self.geometry.N_gauss(i_gauss)

        # Data interpolation at the gauss point
        w_gauss = w.transpose() * Ng

        fluxes_gauss = F_gauss[i_gauss]
        KratosMultiphysics.Logger.Print("    - Substitution in the variational formulation")
        KratosSympy.SubstituteMatrixValue(rv_gauss, F, fluxes_gauss)
        KratosSympy.SubstituteMatrixValue(rv_gauss, V, w_gauss)

        return rv_gauss

    def _ComputeResidualAtGaussPointFree(self, grad_U, H, i_gauss, rv_gauss, sc_nodes, sc_params, U, Ug, V, w):
        KratosMultiphysics.Logger.Print("    Gauss point: {}".format(i_gauss))

        # Get Gauss point geometry data
        Ng = self.geometry.N_gauss(i_gauss)

        # Data interpolation at the gauss point
        w_gauss = w.transpose() * Ng

        U_gauss = U.transpose() * Ng
        grad_U_gauss = grad_U[i_gauss]

        alpha_sc_gauss = (sc_nodes.alpha.transpose()*Ng)[0]
        mu_sc_gauss    = (sc_nodes.mu.transpose()*Ng)[0]
        beta_sc_gauss  = (sc_nodes.beta.transpose()*Ng)[0]
        lamb_sc_gauss  = (sc_nodes.lamb.transpose()*Ng)[0]

        KratosMultiphysics.Logger.Print("    - Substitution in the variational formulation")
        KratosSympy.SubstituteMatrixValue(rv_gauss, Ug, U_gauss)
        KratosSympy.SubstituteMatrixValue(rv_gauss, H, grad_U_gauss)
        KratosSympy.SubstituteMatrixValue(rv_gauss, V, w_gauss)

        KratosSympy.SubstituteScalarValue(rv_gauss, sc_params.alpha, alpha_sc_gauss)
        KratosSympy.SubstituteScalarValue(rv_gauss, sc_params.mu, mu_sc_gauss)
        KratosSympy.SubstituteScalarValue(rv_gauss, sc_params.beta, beta_sc_gauss)
        KratosSympy.SubstituteScalarValue(rv_gauss, sc_params.lamb, lamb_sc_gauss)

        return rv_gauss

    def _BuildGradientsMatrixArray(self):
        dim = self.geometry.ndims
        n_gauss = self.geometry.ngauss
        block_size = self.geometry.blocksize

        grad_U = [defs.ZeroMatrix(block_size, dim) for _ in range(n_gauss)]

        for g in range(n_gauss):
            grad_U[g][   0, :] = defs.Vector('data.gradients[{}].density'.format(g),  dim).T
            grad_U[g][1:-1, :] = defs.Matrix('data.gradients[{}].momentum'.format(g), dim, dim).T
            grad_U[g][  -1, :] = defs.Vector('data.gradients[{}].total_energy'.format(g), dim).T

        return grad_U

    def _DefineImposedFluxesGaussPointMatrix(self):
        dim = self.geometry.ndims
        n_gauss = self.geometry.ngauss
        block_size = self.geometry.blocksize

        F_gauss = [defs.ZeroVector(block_size) for _ in range(n_gauss)]

        for g in range(n_gauss):
            F_gauss[g][   0, 0] = sympy.Symbol('data.fluxes[{}].density'.format(g))
            F_gauss[g][1:-1, 0] = defs.Vector('data.fluxes[{}].momentum'.format(g), dim)
            F_gauss[g][  -1, 0] = sympy.Symbol('data.fluxes[{}].total_energy'.format(g))

        return F_gauss


    def _ComputeLHSandRHS(self, rv_tot, U, w):
        KratosMultiphysics.Logger.Print("    Compute RHS")

        testfunc = defs.ZeroVector(self.geometry.ndofs)
        for i in range(0, self.geometry.nnodes):
            for j in range(0, self.geometry.blocksize):
                testfunc[i*self.geometry.blocksize + j] = w[i, j]
        rhs = KratosSympy.Compute_RHS(rv_tot.copy(), testfunc, self.simplify)

        if self.is_explicit:
            return(None, rhs)

        KratosMultiphysics.Logger.Print("    Compute LHS")
        dofs = defs.ZeroVector(self.geometry.ndofs)
        for i in range(0, self.geometry.nnodes):
            for j in range(0, self.geometry.blocksize):
                dofs[i*self.geometry.blocksize + j] = U[i, j]
        lhs = KratosSympy.Compute_LHS(rhs, testfunc, dofs, self.simplify)

        return (lhs, rhs)


    def _OutputLHSandRHS(self, lhs, rhs):
        # Reading and filling the template file
        KratosMultiphysics.Logger.Print("    Substituting outstring from {}".format(self.template_filename))

        if self.geometry.symbolic_integration:
            rhs_name = "rRightHandSideBoundedVector"
            lhs_name = "rLeftHandSideBoundedVector"
        else:
            rhs_name = "rhs_gauss"
            lhs_name = "lhs_gauss"

        target = "//substitute_rhs_{}D_fluxes".format(self.geometry.ndims)
        self.CollectAndReplace(target, rhs, rhs_name)

        if self.is_explicit:
            return

        target = "//substitute_rhs_{}D_fluxes".format(self.geometry.ndims)
        self.CollectAndReplace(target, lhs, lhs_name)

    def _GenerateNeumann(self):
        n_nodes = self.geometry.nnodes
        block_size = self.geometry.blocksize

        # Construction of the variational equation
        flux = defs.Vector('F', block_size, real=True)
        V   = defs.Vector('V', block_size, real=True)              # Test function

        # Test functions defintion
        w = defs.Matrix('w', n_nodes, block_size, real=True)

        # Calculate the Gauss point residual
        rv = self._ComputeVariationalFormulation(flux, V)

        # Algebraic form calculation
        flux_gauss = self._DefineImposedFluxesGaussPointMatrix()

        # Substitution of the discretized values at the gauss points
        # Loop and accumulate the residual in each Gauss point
        rv_tot = defs.ZeroMatrix(1, 1)
        KratosMultiphysics.Logger.Print(" - Substitution of the discretized values at the gauss points")
        for i_gauss in self.geometry.SymbolicIntegrationPoints():
            rv_gauss = rv.copy()
            rv_tot += self._ComputeResidualAtGaussPointNeumann(flux, flux_gauss, i_gauss, rv_gauss, V, w)

        # Compute functionals and substitute text in file
        (lhs, rhs) = self._ComputeLHSandRHS(rv_tot, None, w)
        self._OutputLHSandRHS(lhs, rhs)

    def _GenerateFree(self):
        dim = self.geometry.ndims
        n_nodes = self.geometry.nnodes
        block_size = self.geometry.blocksize

        params = FormulationParameters(self.geometry, self.write_language)

        # Unknowns
        U = defs.Matrix('data.U', n_nodes, block_size, real=True)
        grad_U = self._BuildGradientsMatrixArray()

        # Nodal artificial magnitudes
        sc_nodes = ShockCapturingNodalParameters(self.geometry)
        sc_params = ShockCapturingParameters() if self.shock_capturing else None

        # Construction of the variational equation
        Ug  = defs.Vector('Ug', block_size, real=True)             # Dofs vector
        H   = defs.Matrix('H', block_size, dim, real=True)         # Gradient of Ug
        V   = defs.Vector('V', block_size, real=True)              # Test function
        n   = defs.Vector('data.unit_normal', dim, real=True)      # Normal vector

        # Test functions defintion
        w = defs.Matrix('w', n_nodes, block_size, real=True)

        # Imposed fluxes
        flux = self._ComputeFluxes(n, Ug, H, params, sc_params)

        # Calculate the Gauss point residual
        rv = self._ComputeVariationalFormulation(flux, V)

        # Algebraic form calculation

        # Substitution of the discretized values at the gauss points
        # Loop and accumulate the residual in each Gauss point
        rv_tot = defs.ZeroMatrix(1, 1)
        KratosMultiphysics.Logger.Print(" - Substitution of the discretized values at the gauss points")
        for i_gauss in self.geometry.SymbolicIntegrationPoints():
            rv_gauss = rv.copy()
            rv_tot += self._ComputeResidualAtGaussPointFree(grad_U, H, i_gauss, rv_gauss, sc_nodes, sc_params, U, Ug, V, w)

        # Compute functionals and substitute text in file
        (lhs, rhs) = self._ComputeLHSandRHS(rv_tot, U, w)
        self._OutputLHSandRHS(lhs, rhs)

    def Generate(self):
        self.WriteWarningMessage()

        KratosMultiphysics.Logger.Print("Computing geometry: {}".format(self.geometry))

        if self.is_neumann:
            self._GenerateNeumann()
        else:
            self._GenerateFree()
