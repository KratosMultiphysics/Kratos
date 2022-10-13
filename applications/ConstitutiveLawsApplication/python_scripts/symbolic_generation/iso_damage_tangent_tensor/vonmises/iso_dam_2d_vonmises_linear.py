from KratosMultiphysics.sympy_fe_utilities import *
from sympy import *


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
threshold = Symbol("threshold")
Gf = Symbol("Gf")
characteristic_length = Symbol("characteristic_length")


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


# ONLY FOR VON MISES!!
J2 =  0.5 * (Deviator[0]**2.0 + Deviator[1]**2.0 + pmean**2.0) + Deviator[2]**2.0;
VonMisesStress = sqrt(3.0*J2)


# Assuming Von Mises stress!
A = -threshold**2 / (2.0 * Young * Gf / characteristic_length);


# only for linear softening!
damage = (1.0 - threshold / VonMisesStress) / (1.0 + A);

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