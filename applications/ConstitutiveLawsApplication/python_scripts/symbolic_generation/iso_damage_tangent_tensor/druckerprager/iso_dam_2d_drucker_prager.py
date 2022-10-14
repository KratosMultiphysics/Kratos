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
threshold = Symbol("threshold") # f_t
Gf = Symbol("Gf")
sin_phi = Symbol("sin_phi")
characteristic_length = Symbol("characteristic_length")


#plane strain
c0 = Young / ((1.00 + nu)*(1 - 2 * nu));
c1 = (1.00 - nu)*c0;
c2 = c0 * nu;
c3 = (0.5 - nu)*c0;

Seff[0] = c1 * Strain0 + c2 *Strain1;
Seff[1] = c2 * Strain0 + c1 * Strain1;
Seff[2] = c3 * Strain2;



pmean = (Seff[0] + Seff[1]) / 3.0
Deviator[0] = Seff[0] - pmean
Deviator[1] = Seff[1] - pmean
Deviator[2] = Seff[2]


J2 =  0.5 * (Deviator[0]**2.0 + Deviator[1]**2.0 + pmean**2.0) + Deviator[2]**2.0;

root_3 = math.sqrt(3.0)

CFL = -root_3 * (3.0 - sin_phi) / (3.0 * sin_phi - 3.0)
TEN0 = 6.0 * pmean * sin_phi / (root_3 * (3.0 - sin_phi)) + sqrt(J2)
DruckerPragerStress = CFL*TEN0


A = 1.0 / (Gf * Young / (characteristic_length * threshold**2) - 0.5)

# for drucker prager
damage_threshold = abs(threshold * (3.0 + sin_phi) / (3.0 * sin_phi - 3.0))


damage = 1.0 - (damage_threshold / DruckerPragerStress) * exp(A * (1.0 - DruckerPragerStress / damage_threshold))


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
