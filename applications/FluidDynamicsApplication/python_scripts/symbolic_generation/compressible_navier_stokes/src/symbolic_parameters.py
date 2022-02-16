import sympy
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

    There are two interpolation methods.

    For a given `f:U->V`, where `U` is the
    set of conservative variables and `V` is the set of primitive ones,

    1. Nodal interpolation computes V at the nodes and then interpolates it:
    ```
    V(x) = Σ_i N_i(x)·f(U(x_i))     for i=1,2...nnodes
    ```

    2. Gaussian interpolation interpolates U and then computes V at the gauss points:
    ```
    V(x) = f(Σ_i N_i(x)·U(x_i))     for i=1,2...nnodes
    ```
    """

    def __init__(self, geometry):
        self.P = sympy.Symbol('pressure', real=True)
        self.V = defs.Vector('velocity', geometry.ndims, real=True)
        self.T = sympy.Symbol('temperature', real=True)
        self.grad_P = defs.Vector('grad_pressure', geometry.ndims, real=True)
        self.grad_V = defs.Matrix('grad_velocity', geometry.ndims, geometry.ndims, real=True)
        self.grad_T = defs.Vector('grad_temperature', geometry.ndims, real=True)

        self.nnodes = geometry.nnodes
        self.ndims = geometry.ndims

    def AsVector(self):
        V = [[self.P]]
        for v in self.V:
            V.append([v])
        V.append([self.T])
        return sympy.Matrix(V)
