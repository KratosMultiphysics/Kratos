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