## Dictionary of the constant parameters used in the Variational Formulation
import sympy
import KratosMultiphysics.sympy_fe_utilities as KratosSympy
from dataclasses import dataclass

@dataclass
class FormulationParameters:
    mu: sympy.Symbol
    h : sympy.Symbol
    lambda_ : sympy.Symbol
    c_v : sympy.Symbol
    gamma : sympy.Symbol
    stab_c1 : sympy.Symbol
    stab_c2 : sympy.Symbol
    stab_c3 : sympy.Symbol
    dim: int

    def __init__(self, geometry):
        self.mu      = sympy.Symbol('mu', positive = True),       # Dynamic viscosity
        self.h       = sympy.Symbol('h', positive = True),        # Element size
        self.lambda_ = sympy.Symbol('lambda', positive = True),   # Thermal Conductivity of the fluid
        self.c_v     = sympy.Symbol('c_v', positive = True),      # Specific Heat at Constant volume
        self.gamma   = sympy.Symbol('gamma',positive = True),     # Gamma (Cp/Cv)
        self.stab_c1 = sympy.Symbol('stab_c1', positive = True),  # Algorithm constant
        self.stab_c2 = sympy.Symbol('stab_c2', positive = True),  # Algorithm constant
        self.stab_c3 = sympy.Symbol('stab_c3', positive = True),  # Algorithm constant
        self.dim = geometry.ndims

@dataclass
class ShockCapturingParameters:
    alpha: sympy.Symbol
    mu: sympy.Symbol
    beta : sympy.Symbol
    lambda_ : sympy.Symbol

    def __init__(self):
        self.alpha   = sympy.Symbol('alpha_sc', positive = True) # Artificial density diffusivity for shock capturing
        self.mu      = sympy.Symbol('mu_sc', positive = True) # Artificial dynamic viscosity for shock capturing
        self.beta    = sympy.Symbol('beta_sc', positive = True) # Artificial bulk viscosity for shock capturing
        self.lambda_ = sympy.Symbol('lamb_sc', positive = True) # Artificial thermal conductivity for shock capturing

class ShockCapturingNodalParameters:
    alpha: sympy.Symbol
    mu: sympy.Symbol
    beta : sympy.Symbol
    lambda_ : sympy.Symbol

    def __init__(self, geometry):
        self.alpha   = KratosSympy.DefineVector('alpha_sc_nodes', geometry.nnodes),# Nodal artificial mass diffusivity
        self.mu      = KratosSympy.DefineVector('mu_sc_nodes', geometry.nnodes),   # Nodal artificial dynamic viscosity
        self.beta    = KratosSympy.DefineVector('beta_sc_nodes', geometry.nnodes), # Nodal artificial bulk viscosity
        self.lambda_ = KratosSympy.DefineVector('lamb_sc_nodes', geometry.nnodes)  # Nodal artificial bulk viscosity