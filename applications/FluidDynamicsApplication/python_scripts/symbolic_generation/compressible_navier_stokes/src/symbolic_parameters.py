import sympy
import KratosMultiphysics.sympy_fe_utilities as KratosSympy
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src.defines \
    import CompressibleNavierStokesDefines as defs

class FormulationParameters:
    "Dictionary of the constant parameters used in the Variational Formulation."

    def __init__(self, geometry, output_lang):
        self.mu      = sympy.Symbol('data.mu', positive = True)       # Dynamic viscosity
        self.h       = sympy.Symbol('data.h', positive = True)        # Element size
        self.c_v     = sympy.Symbol('data.c_v', positive = True)      # Specific Heat at Constant volume
        self.gamma   = sympy.Symbol('data.gamma',positive = True)     # Gamma (Cp/Cv)
        self.stab_c1 = sympy.Symbol('stab_c1', positive = True)  # Algorithm constant
        self.stab_c2 = sympy.Symbol('stab_c2', positive = True)  # Algorithm constant
        self.stab_c3 = sympy.Symbol('stab_c3', positive = True)  # Algorithm constant
        self.dim = geometry.ndims

        lamb = 'data.lambda_' if output_lang == 'python' else 'data.lambda'
        self.lamb    = sympy.Symbol(lamb, positive = True)   # Thermal Conductivity of the fluid

class ShockCapturingParameters:
    def __init__(self):
        self.alpha = sympy.Symbol('data.alpha_sc', positive = True) # Artificial density diffusivity for shock capturing
        self.mu    = sympy.Symbol('data.mu_sc', positive = True)    # Artificial dynamic viscosity for shock capturing
        self.beta  = sympy.Symbol('data.beta_sc', positive = True)  # Artificial bulk viscosity for shock capturing
        self.lamb  = sympy.Symbol('data.lamb_sc', positive = True)  # Artificial thermal conductivity for shock capturing

class ShockCapturingNodalParameters:
    def __init__(self, geometry):
        self.alpha = defs.Vector('data.alpha_sc_nodes', geometry.nnodes, positive=True) # Nodal artificial mass diffusivity
        self.mu    = defs.Vector('data.mu_sc_nodes', geometry.nnodes, positive=True)    # Nodal artificial dynamic viscosity
        self.beta  = defs.Vector('data.beta_sc_nodes', geometry.nnodes, positive=True)  # Nodal artificial bulk viscosity
        self.lamb  = defs.Vector('data.lamb_sc_nodes', geometry.nnodes, positive=True)  # Nodal artificial bulk viscosity

class PrimitiveMagnitudes:
    """
    Primitive variables and their gradients.
    Gradients defined as:
    ```
    grad_f[i, j] := df_j/dx_i.
    ```
    """

    def __init__(self, Ug, grad_Ug, params):
        (self.P, self.V, self.T) = self._PrimitivesFromConservatives(Ug, params)
        (self.grad_P, self.grad_V, self.grad_T) = self._PrimitiveGradientsFromConservatives(Ug, grad_Ug, params)

    @classmethod
    def _PrimitivesFromConservatives(cls, U, params):
        rho = U[0]
        mom = sympy.Matrix(U[1:-1])
        e_tot = U[-1]

        V = cls._velocity(rho, mom)
        T = cls._temperature(rho, V, e_tot, params)
        P = cls._pressure(rho, T, params)

        V.simplify()
        T.simplify()
        P.simplify()

        return (P, V, T)

    @classmethod
    def _PrimitiveGradientsFromConservatives(cls, U, grad_U, params):
        rho = U[0]
        mom = sympy.Matrix(U[1:-1])
        e_tot = U[-1]

        V = cls._velocity(rho, mom)
        T = cls._temperature(rho, V, e_tot, params)

        grad_rho = sympy.Matrix(grad_U[0, :])                   # d rho / dx_i
        grad_mom = sympy.Matrix(grad_U[1:-1, :])                # d V_j / dx_i
        grad_e_tot = sympy.Matrix(grad_U[-1, :])                # d e_tot / dx_i

        grad_V = cls._velocity_gradient(rho, mom, grad_rho, grad_mom)
        grad_T = cls._temperature_gradient(rho, V, T, grad_rho, grad_e_tot, grad_V, params)
        grad_P = cls._pressure_gradient(rho, grad_rho, T, grad_T, params)

        grad_V.simplify()
        grad_T.simplify()
        grad_P.simplify()

        return (grad_P.T, grad_V.T, grad_T.T)

    @classmethod
    def _velocity(cls, rho, mom):
        return mom/rho

    @classmethod
    def _temperature(cls, rho, vel, e_tot, params):
        e_kinetic = 0.5 * rho * sum([v**2 for v in vel])
        return (e_tot - e_kinetic) / (rho * params.c_v)

    @classmethod
    def _pressure(cls, rho, T, params):
        R = (params.gamma - 1) * params.c_v
        return rho * R * T

    @classmethod
    def _velocity_gradient(cls, rho, mom, grad_rho, grad_mom):
        """
        Velocity gradient. Gradients defined as:

        grad_f := df_j/dx_i

        """
        return (grad_mom*rho - mom*grad_rho) / rho**2


    @classmethod
    def _temperature_gradient(cls, rho, vel, T, grad_rho, grad_e_tot, grad_vel, params):
        """
        Temperature gradient.  Gradients defined as:

        grad_f := df_j/dx_i

        Explanation
        ---
        Temperature is defined as:
        ```
            T = (e_t - e_k) / (rho * c_v)
            [with e_k = 1/2 * rho * V²]
        ```
        where `e_k` is the kinetic energy.

        The gradient is:
        ```
            grad(T) = grad(e_t - e_k) / (rho * c_v)   -   (e_t - e_k)/(rho²*c_v) * grad(rho)
        ```
        Rearranging yields:
        ```
            grad(T) = (grad(e_t) - grad(e_k)) / (rho * c_v)   -   T/rho * grad(rho)
            [with grad(e_k) = rho*transp(V)*grad(V) + 1/2 * V² * grad(rho)]
        ```

        """
        V2 = vel.transpose() * vel
        grad_e_k = rho*vel.transpose()*grad_vel + 0.5 * V2 * grad_rho
        return (grad_e_tot - grad_e_k)/(rho*params.c_v) - T/rho * grad_rho

    @classmethod
    def _pressure_gradient(cls, rho, grad_rho, T, grad_T, params):
        """
        Pressure gradient.  Gradients defined as:

        grad_f := df_j/dx_i

        """
        R = (params.gamma - 1) * params.c_v
        return R * (grad_rho*T + rho*grad_T)
