from KratosMultiphysics.sympy_fe_utilities import *
from sympy import *

"""
Auxiliary file used as a library for automatic
differentiation files
"""

# 2D calculations
def SetUp2DProblem():
    Strain0 = Symbol("r_strain[0]")
    Strain1 = Symbol("r_strain[1]")
    Strain2 = Symbol("r_strain[2]")

    Ct = DefineMatrix('r_Ct', 3, 3)

    # Stress (effective and integrated and deviatoric)
    Seff = DefineVector('Seff',3)
    Stress = DefineVector('Stress', 3)
    Deviator = DefineVector('Deviator', 3)

    # material parameters
    Young = Symbol("Young")
    nu = Symbol("nu")
    threshold = Symbol("threshold")
    Gf = Symbol("Gf")
    characteristic_length = Symbol("characteristic_length")

    return Strain0, Strain1, Strain2, Seff, Stress, Deviator, Ct, Young, nu, threshold, Gf, characteristic_length


def ComputePredictorStressVector2D(Young, nu, Strain0, Strain1, Strain2):
    c0 = Young / ((1.00 + nu)*(1 - 2 * nu));
    c1 = (1.00 - nu)*c0;
    c2 = c0 * nu;
    c3 = (0.5 - nu)*c0;

    S0 = c1 * Strain0 + c2 *Strain1;
    S1 = c2 * Strain0 + c1 * Strain1;
    S2 = c3 * Strain2;

    return S0, S1, S2

def ComputeJ2Invariant2D(Deviator, pmean):
    return 0.5 * (Deviator[0]**2.0 + Deviator[1]**2.0 + pmean**2.0) + Deviator[2]**2.0

def ComputePmean2D(StressVector):
    return (StressVector[0] + StressVector[1]) / 3.0

def ComputeDamageLinear(UniaxialStress, Threshold, A):
    return (1.0 - Threshold / UniaxialStress) / (1.0 + A)

def ComputeAParameterLinear(Threshold, Gf, L, Young):
    return -Threshold**2 / (2.0 * Young * Gf / L)

