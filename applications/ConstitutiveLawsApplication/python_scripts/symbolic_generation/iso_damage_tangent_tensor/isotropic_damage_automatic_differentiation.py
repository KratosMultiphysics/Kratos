from KratosMultiphysics.sympy_fe_utilities import *
from sympy import *
import math

"""
Auxiliary file used as a library for automatic
differentiation files
"""

# 2D calculations
def SetUp2DProblem():
    Strain0 = Symbol("r_strain[0]")
    Strain1 = Symbol("r_strain[1]")
    Strain2 = Symbol("r_strain[2]")

    Ct = sympy.zeros(3, 3)

    # Stress (effective and integrated and deviatoric)
    Seff = sympy.zeros(3)
    Stress = sympy.zeros(3)
    Deviator = sympy.zeros(3)

    # material parameters
    Young = Symbol("Young")
    nu = Symbol("nu")
    threshold = Symbol("threshold")
    Gf = Symbol("Gf")
    characteristic_length = Symbol("characteristic_length")

    return Strain0, Strain1, Strain2, Seff, Stress, Deviator, Ct, Young, nu, threshold, Gf, characteristic_length

def SetUp3DProblem():
    Strain0 = Symbol("r_strain[0]")
    Strain1 = Symbol("r_strain[1]")
    Strain2 = Symbol("r_strain[2]")
    Strain3 = Symbol("r_strain[3]")
    Strain4 = Symbol("r_strain[4]")
    Strain5 = Symbol("r_strain[5]")

    Ct = sympy.zeros(6, 6)

    # Stress (effective and integrated and deviatoric)
    Seff = sympy.zeros(6)
    Stress = sympy.zeros(6)
    Deviator = sympy.zeros(6)

    # material parameters
    Young = Symbol("Young")
    nu = Symbol("nu")
    threshold = Symbol("threshold")
    Gf = Symbol("Gf")
    characteristic_length = Symbol("characteristic_length")

    return Strain0, Strain1, Strain2, Strain3, Strain4, Strain5, Seff, Stress, Deviator, Ct, Young, nu, threshold, Gf, characteristic_length


def ComputePredictorStressVector2D(Young, nu, Strain0, Strain1, Strain2):
    c0 = Young / ((1.0 + nu)*(1.0 - 2.0 * nu))
    c1 = (1.0 - nu)*c0
    c2 = c0 * nu
    c3 = (0.5 - nu)*c0

    S0 = c1 * Strain0 + c2 *Strain1
    S1 = c2 * Strain0 + c1 * Strain1
    S2 = c3 * Strain2

    return S0, S1, S2
def ComputePredictorStressVector3D(Young, nu, Strain0, Strain1, Strain2, Strain3, Strain4, Strain5):
    c1 = Young / ((1.0 + nu) * (1.0 - 2.0 * nu))
    c2 = c1 * (1.0 - nu)
    c3 = c1 * nu
    c4 = c1 * 0.5 * (1.0 - 2.0 * nu)

    S0 = c2 * Strain0 + c3 * Strain1 + c3 * Strain2
    S1 = c3 * Strain0 + c2 * Strain1 + c3 * Strain2
    S2 = c3 * Strain0 + c3 * Strain1 + c2 * Strain2
    S3 = c4 * Strain3
    S4 = c4 * Strain4
    S5 = c4 * Strain5

    return S0, S1, S2, S3, S4, S5

def ComputeJ2Invariant2D(Deviator, pmean):
    return 0.5 * (Deviator[0]**2.0 + Deviator[1]**2.0 + pmean**2.0) + Deviator[2]**2.0
def ComputeJ2Invariant3D(Deviator, pmean):
    return 0.5*(Deviator[0]**2+Deviator[1]**2+Deviator[2]**2) + (Deviator[3]**2+Deviator[4]**2+Deviator[5]**2)

def ComputeJ3Invariant2D(Deviator):
    return Deviator[0] * Deviator[1] - (Deviator[2]**2)
def ComputeJ3Invariant3D(Deviator):
    return Deviator[0] * (Deviator[1] * Deviator[2] - Deviator[4] * Deviator[4]) + Deviator[3] * (-Deviator[3] * Deviator[2] + Deviator[5] * Deviator[4]) + Deviator[5] * (Deviator[3] * Deviator[4] - Deviator[5] * Deviator[1])

def ComputePmean2D(StressVector):
    return (StressVector[0] + StressVector[1]) / 3.0
def ComputePmean3D(StressVector):
    return (StressVector[0] + StressVector[1] + StressVector[2]) / 3.0

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
---------------------------------------------------------------------------------------------------
"""
dimension = 2 # 2, 3
yield_surface = "DruckerPrager" # DruckerPrager, ModifiedMohrCoulomb, Rankine, VonMises
softening = "Linear" # Exponential, Linear
"""
---------------------------------------------------------------------------------------------------
"""

# Common variables
mode = "c"

J2 = 0.0; uniaxial_stress = 0.0; A = 0.0; damage = 0.0;
if (dimension == 2):
    Strain0, Strain1, Strain2, Seff, Stress, Deviator, Ct, Young, nu, threshold, Gf, characteristic_length = SetUp2DProblem()

    Seff[0], Seff[1], Seff[2] = ComputePredictorStressVector2D(Young, nu, Strain0, Strain1, Strain2)

    pmean = ComputePmean2D(Seff)

    Deviator[0] = Seff[0] - pmean
    Deviator[1] = Seff[1] - pmean
    Deviator[2] = Seff[2]

    if (yield_surface == "VonMises"):
        J2 =  ComputeJ2Invariant2D(Deviator, pmean)
        VonMisesStress = sqrt(3.0*J2)
        damage = ComputeDamage(VonMisesStress, threshold, Gf, characteristic_length, Young, softening)
    elif (yield_surface == "Rankine"):
        rankine_yield = 0.5 * (Seff[0] + Seff[1]) + sqrt((0.5 * (Seff[0] - Seff[1]))**2  + (Seff[2])**2)
        damage = ComputeDamage(rankine_yield, threshold, Gf, characteristic_length, Young, softening)
    elif (yield_surface == "DruckerPrager"):
        sin_phi = Symbol("sin_phi")
        J2 =  ComputeJ2Invariant2D(Deviator, pmean)
        root_3 = math.sqrt(3.0)
        CFL = -root_3 * (3.0 - sin_phi) / (3.0 * sin_phi - 3.0)
        TEN0 = 6.0 * pmean * sin_phi / (root_3 * (3.0 - sin_phi)) + sqrt(J2)
        DruckerPragerStress = CFL*TEN0
        damage_threshold = abs(threshold * (3.0 + sin_phi) / (3.0 * sin_phi - 3.0))
        if (softening == "Exponential"):
            A = 1.0 / (Gf * Young / (characteristic_length * threshold**2) - 0.5)
            damage_threshold = abs(threshold * (3.0 + sin_phi) / (3.0 * sin_phi - 3.0))
            damage = 1.0 - (damage_threshold / DruckerPragerStress) * exp(A * (1.0 - DruckerPragerStress / damage_threshold))
        else:
            A = -threshold**2 / (2.0 * Young * Gf / characteristic_length)
            damage = (1.0 - damage_threshold / DruckerPragerStress) / (1.0 + A)
    elif (yield_surface == "ModifiedMohrCoulomb"):
        phi = Symbol("phi")
        threshold_tension = Symbol("threshold_tension") # f_t
        threshold_compression = Symbol("threshold_compression") # f_t
        J2 =  ComputeJ2Invariant2D(Deviator, pmean)
        R = abs(threshold_compression / threshold_tension)
        Rmohr = (tan((math.pi / 4.0) + phi / 2.0))**2
        alpha_r = R / Rmohr
        sin_phi = sin(phi);
        J3 = ComputeJ3Invariant2D(Deviator)
        K1 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) * sin_phi;
        K2 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) / sin_phi;
        K3 = 0.5 * (1.0 + alpha_r) * sin_phi - 0.5 * (1.0 - alpha_r);
        sint3 = (-3.0 * sqrt(3.0) * J3) / (2.0 * J2 * sqrt(J2))
        LodeAngle = asin(sint3) / 3.0;
        ModifiedMohrCoulombStress = (2.0 * tan(math.pi * 0.25 + phi * 0.5) / cos(phi)) * ((3.0*pmean * K3 / 3.0) +sqrt(J2) * (K1 * cos(LodeAngle) - K2 * sin(LodeAngle) * sin_phi / sqrt(3.0)))
        if (softening == "Exponential"):
            A = 1.0 / (Gf * Young * R * R / (characteristic_length * threshold_compression**2) - 0.5)
            damage = 1.0 - (threshold_compression / ModifiedMohrCoulombStress) * exp(A * (1.0 - ModifiedMohrCoulombStress / threshold_compression))
        else:
            A = -threshold_compression**2 / (2.0 * Young * Gf *R*R / characteristic_length)
            damage = (1.0 - threshold_compression / ModifiedMohrCoulombStress) / (1.0 + A)
    else:
        raise Exception("yield_surface not correctly defined, options available: VonMises, DruckerPrager, Rankine and ModifiedMohrCoulomb")

    Stress = (1.0 - damage)*Seff
    Strain = [Strain0, Strain1, Strain2]

    # Ct = dS/dE
    for i in range(3):
        for j in range(3):
            Ct[i,j] = Stress[i].diff(Strain[j])

    out = OutputMatrix_CollectingFactors(Ct, "r_Ct", mode)
    print(out)

else: # 3D
    Strain0, Strain1, Strain2, Strain3, Strain4, Strain5, Seff, Stress, Deviator, Ct, Young, nu, threshold, Gf, characteristic_length = SetUp3DProblem()
    Seff[0], Seff[1], Seff[2], Seff[3], Seff[4], Seff[5] = ComputePredictorStressVector3D(Young, nu, Strain0, Strain1, Strain2, Strain3, Strain4, Strain5)
    pmean = ComputePmean3D(Seff)

    Deviator[0] = Seff[0] - pmean
    Deviator[1] = Seff[1] - pmean
    Deviator[2] = Seff[2] - pmean
    Deviator[3] = Seff[3]
    Deviator[4] = Seff[4]
    Deviator[5] = Seff[5]

    if (yield_surface == "VonMises"):
        J2 =  ComputeJ2Invariant3D(Deviator, pmean)
        VonMisesStress = sqrt(3.0*J2)
        damage = ComputeDamage(VonMisesStress, threshold, Gf, characteristic_length, Young, softening)
    elif (yield_surface == "Rankine"):
        norm_stress_vector = DefineVector('norm_stress_vector', 6)
        PrincipalStresses = DefineVector('PrincipalStresses', 3)
        norm_frobenius = sqrt(Seff[0]**2 + Seff[1]**2 + Seff[2]**2+Seff[3]**2 + Seff[4]**2 + Seff[5]**2)
        norm_stress_vector[0] = Seff[0] / norm_frobenius
        norm_stress_vector[1] = Seff[1] / norm_frobenius
        norm_stress_vector[2] = Seff[2] / norm_frobenius
        norm_stress_vector[3] = Seff[3] / norm_frobenius
        norm_stress_vector[4] = Seff[4] / norm_frobenius
        norm_stress_vector[5] = Seff[5] / norm_frobenius

        I1 = norm_stress_vector[0] + norm_stress_vector[1] + norm_stress_vector[2]
        I2 = (norm_stress_vector[0] + norm_stress_vector[2]) * norm_stress_vector[1] + norm_stress_vector[0] * norm_stress_vector[2] -norm_stress_vector[3] * norm_stress_vector[3] - norm_stress_vector[4] * norm_stress_vector[4] - norm_stress_vector[5] * norm_stress_vector[5];
        I3 = (norm_stress_vector[1] * norm_stress_vector[2] - norm_stress_vector[4] * norm_stress_vector[4]) * norm_stress_vector[0] - norm_stress_vector[1] * norm_stress_vector[5] * norm_stress_vector[5] - norm_stress_vector[2] * norm_stress_vector[3] * norm_stress_vector[3] + 2.0 * norm_stress_vector[3] * norm_stress_vector[4] * norm_stress_vector[5]
        II1 = I1**2
        R = (2.0 * II1 * I1 - 9.0 * I2 * I1 + 27.0 * I3) / 54.0;
        Q = (3.0 * I2 - II1) / 9.0;
        cos_phi = R / (sqrt(-Q**3));
        phi = acos(cos_phi);
        phi_3 = phi / 3.0;
        aux1 = 2.0 * sqrt(-Q);
        aux2 = I1 / 3.0;
        deg_120 = 2.0 / 3.0 * math.pi;
        for i in range(3):
            PrincipalStresses[i] = norm_frobenius * (aux2 + aux1 * cos(phi_3 + deg_120 * i));

        rankine_yield = PrincipalStresses[0]
        damage = ComputeDamage(rankine_yield, threshold, Gf, characteristic_length, Young, softening)
    elif (yield_surface == "DruckerPrager"):
        sin_phi = Symbol("sin_phi")
        J2 =  ComputeJ2Invariant3D(Deviator, pmean)
        root_3 = math.sqrt(3.0)
        CFL = -root_3 * (3.0 - sin_phi) / (3.0 * sin_phi - 3.0)
        TEN0 = 6.0 * pmean * sin_phi / (root_3 * (3.0 - sin_phi)) + sqrt(J2)
        DruckerPragerStress = CFL*TEN0
        damage_threshold = abs(threshold * (3.0 + sin_phi) / (3.0 * sin_phi - 3.0))
        if (softening == "Exponential"):
            A = 1.0 / (Gf * Young / (characteristic_length * threshold**2) - 0.5)
            damage_threshold = abs(threshold * (3.0 + sin_phi) / (3.0 * sin_phi - 3.0))
            damage = 1.0 - (damage_threshold / DruckerPragerStress) * exp(A * (1.0 - DruckerPragerStress / damage_threshold))
        else:
            A = -threshold**2 / (2.0 * Young * Gf / characteristic_length)
            damage = (1.0 - damage_threshold / DruckerPragerStress) / (1.0 + A)
    elif (yield_surface == "ModifiedMohrCoulomb"):
        phi = Symbol("phi")
        threshold_tension = Symbol("threshold_tension") # f_t
        threshold_compression = Symbol("threshold_compression") # f_t
        J2 =  ComputeJ2Invariant3D(Deviator, pmean)
        R = abs(threshold_compression / threshold_tension)
        Rmohr = (tan((math.pi / 4.0) + phi / 2.0))**2
        alpha_r = R / Rmohr
        sin_phi = sin(phi);
        J3 = ComputeJ3Invariant3D(Deviator)
        K1 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) * sin_phi;
        K2 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) / sin_phi;
        K3 = 0.5 * (1.0 + alpha_r) * sin_phi - 0.5 * (1.0 - alpha_r);
        sint3 = (-3.0 * sqrt(3.0) * J3) / (2.0 * J2 * sqrt(J2))
        LodeAngle = asin(sint3) / 3.0;
        ModifiedMohrCoulombStress = (2.0 * tan(math.pi * 0.25 + phi * 0.5) / cos(phi)) * ((3.0*pmean * K3 / 3.0) +sqrt(J2) * (K1 * cos(LodeAngle) - K2 * sin(LodeAngle) * sin_phi / sqrt(3.0)))
        if (softening == "Exponential"):
            A = 1.0 / (Gf * Young * R * R / (characteristic_length * threshold_compression**2) - 0.5)
            damage = 1.0 - (threshold_compression / ModifiedMohrCoulombStress) * exp(A * (1.0 - ModifiedMohrCoulombStress / threshold_compression))
        else:
            A = -threshold_compression**2 / (2.0 * Young * Gf *R*R / characteristic_length)
            damage = (1.0 - threshold_compression / ModifiedMohrCoulombStress) / (1.0 + A)
    else:
        raise Exception("yield_surface not correctly defined, options available: VonMises, DruckerPrager, Rankine and ModifiedMohrCoulomb")

    Stress = (1.0 - damage)*Seff
    Strain = [Strain0, Strain1, Strain2, Strain3, Strain4, Strain5]

    # Ct = dS/dE
    for i in range(6):
        for j in range(6):
            Ct[i,j] = Stress[i].diff(Strain[j])

    out = OutputMatrix_CollectingFactors(Ct, "r_Ct", mode)
    print(out)