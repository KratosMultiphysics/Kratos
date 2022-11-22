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

def ComputeDamage(UniaxialStress, Threshold, Gf, L, Young, Softening):
    if (Softening ==  "Linear"):
        A = ComputeAParameterLinear(Threshold, Gf, L, Young)
        return ComputeDamageLinear(UniaxialStress, Threshold, A)
    else: # Exponential
        A = ComputeAParameterExponential(Threshold, Gf, L, Young)
        return ComputeDamageExponential(UniaxialStress, Threshold, A)

def ComputeDamageLinear(UniaxialStress, Threshold, A):
    return (1.0 - Threshold / UniaxialStress) / (1.0 + A)

def ComputeDamageExponential(UniaxialStress, Threshold, A):
    return 1.0 - (Threshold / UniaxialStress) * exp(A * (1.0 - UniaxialStress / Threshold))

def ComputeAParameterLinear(Threshold, Gf, L, Young):
    return -Threshold**2 / (2.0 * Young * Gf / L)

def ComputeAParameterExponential(Threshold, Gf, L, Young):
    return 1.0 / (Gf * Young / (L * Threshold**2) - 0.5)




"""
Here we define the problem in terms of dimension, yield surface and linear/exponential softening
"""

# INPUT DATA
dimension = 2 # 3
yield_surface = "VonMises" # DruckerPrager, ModifiedMohrCoulomb, Rankine
softening = "Linear" # Exponential

# Common variables
mode = "c"

if (dimension == 2):
    Strain0, Strain1, Strain2, Seff, Stress, Deviator, Ct, Young, nu, threshold, Gf, characteristic_length = SetUp2DProblem()
    Seff[0], Seff[1], Seff[2] = ComputePredictorStressVector2D(Young, nu, Strain0, Strain1, Strain2)

    pmean = ComputePmean2D(Seff)
    Deviator[0] = Seff[0] - pmean
    Deviator[1] = Seff[1] - pmean
    Deviator[2] = Seff[2]
    J2 = 0.0; uniaxial_stress = 0.0; A = 0.0; damage = 0.0;
    if (yield_surface == "VonMises"):
        J2 =  ComputeJ2Invariant2D(Deviator, pmean)
        VonMisesStress = sqrt(3.0*J2)
        damage = ComputeDamage(VonMisesStress, threshold, Gf, characteristic_length, Young, softening)

    Stress[0] = (1.0 - damage)*Seff[0]
    Stress[1] = (1.0 - damage)*Seff[1]
    Stress[2] = (1.0 - damage)*Seff[2]

    # Ct = dS/dE
    Ct[0,0] = ((Stress[0]).diff(Strain0))
    Ct[0,1] = ((Stress[0]).diff(Strain1))
    Ct[0,2] = ((Stress[0]).diff(Strain2))
    Ct[1,0] = ((Stress[1]).diff(Strain0))
    Ct[1,1] = ((Stress[1]).diff(Strain1))
    Ct[1,2] = ((Stress[1]).diff(Strain2))
    Ct[2,0] = ((Stress[2]).diff(Strain0))
    Ct[2,1] = ((Stress[2]).diff(Strain1))
    Ct[2,2] = ((Stress[2]).diff(Strain2))

    out = OutputMatrix_CollectingFactors(Ct, "r_Ct", mode)
    print(out)
