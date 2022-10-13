from KratosMultiphysics.sympy_fe_utilities import *
from sympy import *
import math

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
PrincipalStresses = DefineVector('PrincipalStresses', 3)
norm_stress_vector = DefineVector('norm_stress_vector', 6)

# material parameters
Young = Symbol("Young")
nu = Symbol("nu")
threshold = Symbol("threshold")
Gf = Symbol("Gf")
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


pmean = (Seff[0] + Seff[1] + Seff[2]) / 3.0
Deviator[0] = Seff[0] - pmean
Deviator[1] = Seff[1] - pmean
Deviator[2] = Seff[2] - pmean
Deviator[3] = Seff[3]
Deviator[4] = Seff[4]
Deviator[5] = Seff[5]

#
# Rankine yield
#
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


# Assuming Von Mises stress!
A = -threshold**2 / (2.0 * Young * Gf / characteristic_length);


# only for linear softening!
damage = (1.0 - threshold / rankine_yield) / (1.0 + A);

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
# ...

out = OutputMatrix_CollectingFactors(Ct, "r_Ct", mode)
print(out)


# out = OutputVector_CollectingFactors(Stress, "stress", mode)
# print(out)