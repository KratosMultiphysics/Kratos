from KratosMultiphysics.sympy_fe_utilities import *
from sympy import *
import KratosMultiphysics.ConstitutiveLawsApplication.symbolic_generation.iso_damage_tangent_tensor.auxiliary_methods_automatic_differentiation as AuxLib

mode = "c"

Strain0, Strain1, Strain2, Seff, Stress, Deviator, Ct, Young, nu, threshold, Gf, characteristic_length = AuxLib.SetUp2DProblem()


Seff[0], Seff[1], Seff[2] = AuxLib.ComputePredictorStressVector2D(Young, nu, Strain0, Strain1, Strain2)


pmean = AuxLib.ComputePmean2D(Seff)
Deviator[0] = Seff[0] - pmean
Deviator[1] = Seff[1] - pmean
Deviator[2] = Seff[2]

# ONLY FOR VON MISES!!
J2 =  AuxLib.ComputeJ2Invariant2D(Deviator, pmean)
VonMisesStress = sqrt(3.0*J2)

# Assuming Von Mises stress!
A = AuxLib.ComputeAParameterLinear(threshold, Gf, characteristic_length, Young)

# only for linear softening!
damage = AuxLib.ComputeDamageLinear(VonMisesStress, threshold, A)

# # Integrated stress
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