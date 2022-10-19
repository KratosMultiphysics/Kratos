from KratosMultiphysics.sympy_fe_utilities import *
from sympy import *
import math
# import KratosMultiphysics as KM


mode = "c"

# Strain
Strain0 = Symbol("r_strain[0]")
Strain1 = Symbol("r_strain[1]")
Strain2 = Symbol("r_strain[2]")


# Tangent constitutive tensor
Ct = DefineMatrix('r_Ct', 3, 3)

# Stress (effective and integrated and deviatoric)
Seff = DefineVector('Seff',3)
Stress = DefineVector('Stress', 3)
Deviator = DefineVector('Deviator', 3)


# material parameters
Young = Symbol("Young")
nu = Symbol("nu")
threshold_tension = Symbol("threshold_tension") # f_t
threshold_compression = Symbol("threshold_compression") # f_t
Gf = Symbol("Gf")
phi = Symbol("phi") # in rad
characteristic_length = Symbol("characteristic_length")


#plane strain
c0 = Young / ((1.00 + nu)*(1 - 2 * nu));
c1 = (1.00 - nu)*c0;
c2 = c0 * nu;
c3 = (0.5 - nu)*c0;

Seff[0] = c1 * Strain0 + c2 *Strain1;
Seff[1] = c2 * Strain0 + c1 * Strain1;
Seff[2] = c3 * Strain2;


I1 = (Seff[0] + Seff[1])
pmean = I1 / 3.0
Deviator[0] = Seff[0] - pmean
Deviator[1] = Seff[1] - pmean
Deviator[2] = Seff[2]


J2 =  0.5 * (Deviator[0]**2.0 + Deviator[1]**2.0 + pmean**2.0) + Deviator[2]**2.0;

R = abs(threshold_compression / threshold_tension)
Rmohr = (tan((math.pi / 4.0) + phi / 2.0))**2

alpha_r = R / Rmohr
sin_phi = sin(phi);

J3 = Deviator[0] * Deviator[1] - (Deviator[2]**2)


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


# # Ct = dS/dE
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

# out = OutputVector_CollectingFactors(Stress, "stress", mode)
# print(out)
