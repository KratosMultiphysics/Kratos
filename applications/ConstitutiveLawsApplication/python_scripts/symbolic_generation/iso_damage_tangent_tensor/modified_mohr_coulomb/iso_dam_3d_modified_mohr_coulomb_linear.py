from KratosMultiphysics.sympy_fe_utilities import *
from sympy import *
import math
# import KratosMultiphysics as KM


mode = "c"

# Strain
Strain0 = Symbol("r_strain[0]")
Strain1 = Symbol("r_strain[1]")
Strain2 = Symbol("r_strain[2]")
Strain3 = Symbol("r_strain[3]")
Strain4 = Symbol("r_strain[4]")
Strain5 = Symbol("r_strain[5]")

# Tangent constitutive tensor
Ct = DefineMatrix('r_Ct', 6, 6)

# Stress (effective and integrated and deviatoric)
Seff = DefineVector('Seff', 6)
Stress = DefineVector('Stress', 6)
Deviator = DefineVector('Deviator', 6)


# material parameters
Young = Symbol("Young")
nu = Symbol("nu")
threshold_tension = Symbol("threshold_tension") # f_t
threshold_compression = Symbol("threshold_compression") # f_t
Gf = Symbol("Gf")
phi = Symbol("phi") # in rad
characteristic_length = Symbol("characteristic_length")


c1 = Young / ((1.0 + nu) * (1.0 - 2.0 * nu))
c2 = c1 * (1.0 - nu)
c3 = c1 * nu
c4 = c1 * 0.5 * (1.0 - 2.0 * nu)

Seff[0] = c2 * Strain0 + c3 * Strain1 + c3 * Strain2
Seff[1] = c3 * Strain0 + c2 * Strain1 + c3 * Strain2
Seff[2] = c3 * Strain0 + c3 * Strain1 + c2 * Strain2
Seff[3] = c4 * Strain3
Seff[4] = c4 * Strain4
Seff[5] = c4 * Strain5

I1 = (Seff[0] + Seff[1] + Seff[2])
pmean = I1 / 3.0
Deviator[0] = Seff[0] - pmean
Deviator[1] = Seff[1] - pmean
Deviator[2] = Seff[2] - pmean
Deviator[3] = Seff[3]
Deviator[4] = Seff[4]
Deviator[5] = Seff[5]


J2 = 0.5*(Deviator[0]**2+Deviator[1]**2+Deviator[2]**2) + (Deviator[3]**2+Deviator[4]**2+Deviator[5]**2)


J3 = Deviator[0] * (Deviator[1] * Deviator[2] - Deviator[4] * Deviator[4]) + Deviator[3] * (-Deviator[3] * Deviator[2] + Deviator[5] * Deviator[4]) + Deviator[5] * (Deviator[3] * Deviator[4] - Deviator[5] * Deviator[1])

R = abs(threshold_compression / threshold_tension)
Rmohr = (tan((math.pi / 4.0) + phi / 2.0))**2

alpha_r = R / Rmohr
sin_phi = sin(phi);

K1 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) * sin_phi;
K2 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) / sin_phi;
K3 = 0.5 * (1.0 + alpha_r) * sin_phi - 0.5 * (1.0 - alpha_r);

sint3 = (-3.0 * sqrt(3.0) * J3) / (2.0 * J2 * sqrt(J2))

LodeAngle = asin(sint3) / 3.0;

ModifiedMohrCoulombStress = (2.0 * tan(math.pi * 0.25 + phi * 0.5) / cos(phi)) * ((I1 * K3 / 3.0) +sqrt(J2) * (K1 * cos(LodeAngle) - K2 * sin(LodeAngle) * sin_phi / sqrt(3.0)))

A = -threshold_compression**2 / (2.0 * Young * Gf *R*R / characteristic_length);


# only for linear softening!
damage = (1.0 - threshold_compression / ModifiedMohrCoulombStress) / (1.0 + A);

# # Integrated stress
Stress[0] = (1.0 - damage)*Seff[0]
Stress[1] = (1.0 - damage)*Seff[1]
Stress[2] = (1.0 - damage)*Seff[2]
Stress[3] = (1.0 - damage)*Seff[3]
Stress[4] = (1.0 - damage)*Seff[4]
Stress[5] = (1.0 - damage)*Seff[5]

# # Ct = dS/dE
Ct[0,0] = ((Stress[0]).diff(Strain0))
Ct[0,1] = ((Stress[0]).diff(Strain1))
Ct[0,2] = ((Stress[0]).diff(Strain2))
Ct[0,3] = ((Stress[0]).diff(Strain3))
Ct[0,4] = ((Stress[0]).diff(Strain4))
Ct[0,5] = ((Stress[0]).diff(Strain5))

Ct[1,0] = ((Stress[1]).diff(Strain0))
Ct[1,1] = ((Stress[1]).diff(Strain1))
Ct[1,2] = ((Stress[1]).diff(Strain2))
Ct[1,3] = ((Stress[1]).diff(Strain3))
Ct[1,4] = ((Stress[1]).diff(Strain4))
Ct[1,5] = ((Stress[1]).diff(Strain5))

Ct[2,0] = ((Stress[2]).diff(Strain0))
Ct[2,1] = ((Stress[2]).diff(Strain1))
Ct[2,2] = ((Stress[2]).diff(Strain2))
Ct[2,3] = ((Stress[2]).diff(Strain3))
Ct[2,4] = ((Stress[2]).diff(Strain4))
Ct[2,5] = ((Stress[2]).diff(Strain5))

Ct[3,0] = ((Stress[3]).diff(Strain0))
Ct[3,1] = ((Stress[3]).diff(Strain1))
Ct[3,2] = ((Stress[3]).diff(Strain2))
Ct[3,3] = ((Stress[3]).diff(Strain3))
Ct[3,4] = ((Stress[3]).diff(Strain4))
Ct[3,5] = ((Stress[3]).diff(Strain5))

Ct[4,0] = ((Stress[4]).diff(Strain0))
Ct[4,1] = ((Stress[4]).diff(Strain1))
Ct[4,2] = ((Stress[4]).diff(Strain2))
Ct[4,3] = ((Stress[4]).diff(Strain3))
Ct[4,4] = ((Stress[4]).diff(Strain4))
Ct[4,5] = ((Stress[4]).diff(Strain5))

Ct[5,0] = ((Stress[5]).diff(Strain0))
Ct[5,1] = ((Stress[5]).diff(Strain1))
Ct[5,2] = ((Stress[5]).diff(Strain2))
Ct[5,3] = ((Stress[5]).diff(Strain3))
Ct[5,4] = ((Stress[5]).diff(Strain4))
Ct[5,5] = ((Stress[5]).diff(Strain5))

out = OutputMatrix_CollectingFactors(Ct, "r_Ct", mode)
print(out)

# out = OutputVector_CollectingFactors(Stress, "stress", mode)
# print(out)
