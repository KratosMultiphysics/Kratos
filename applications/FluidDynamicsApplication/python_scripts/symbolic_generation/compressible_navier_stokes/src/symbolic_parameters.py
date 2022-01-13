## Dictionary of the constant parameters used in the Variational Formulation
import sympy
from dataclasses import dataclass
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src.defines \
    import DefineVector

@dataclass
class FormulationParameters:
    def __init__(self, geometry):
        self.mu      = sympy.Symbol('data.mu', positive = True)       # Dynamic viscosity
        self.h       = sympy.Symbol('data.h', positive = True)        # Element size
        self.lambda_ = sympy.Symbol('data.lambda', positive = True)   # Thermal Conductivity of the fluid
        self.c_v     = sympy.Symbol('data.c_v', positive = True)      # Specific Heat at Constant volume
        self.gamma   = sympy.Symbol('data.gamma',positive = True)     # Gamma (Cp/Cv)
        self.stab_c1 = sympy.Symbol('stab_c1', positive = True)  # Algorithm constant
        self.stab_c2 = sympy.Symbol('stab_c2', positive = True)  # Algorithm constant
        self.stab_c3 = sympy.Symbol('stab_c3', positive = True)  # Algorithm constant
        self.dim = geometry.ndims

@dataclass
class ShockCapturingParameters:
    def __init__(self):
        self.alpha   = sympy.Symbol('data.alpha_sc', positive = True) # Artificial density diffusivity for shock capturing
        self.mu      = sympy.Symbol('data.mu_sc', positive = True) # Artificial dynamic viscosity for shock capturing
        self.beta    = sympy.Symbol('data.beta_sc', positive = True) # Artificial bulk viscosity for shock capturing
        self.lambda_ = sympy.Symbol('data.lamb_sc', positive = True) # Artificial thermal conductivity for shock capturing

@dataclass
class ShockCapturingNodalParameters:
    def __init__(self, geometry):
        self.alpha   = DefineVector('data.alpha_sc_nodes', geometry.nnodes) # Nodal artificial mass diffusivity
        self.mu      = DefineVector('data.mu_sc_nodes', geometry.nnodes)    # Nodal artificial dynamic viscosity
        self.beta    = DefineVector('data.beta_sc_nodes', geometry.nnodes)  # Nodal artificial bulk viscosity
        self.lambda_ = DefineVector('data.lamb_sc_nodes', geometry.nnodes)  # Nodal artificial bulk viscosity