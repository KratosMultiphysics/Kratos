from KratosMultiphysics.sympy_fe_utilities import *
from sympy import *
# import KratosMultiphysics as KM


mode = "c"

# Strain
# Strain = DefineVector('Strain', 6)
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

# ONLY FOR VON MISES!!
J2 = 0.5*(Deviator[0]**2+Deviator[1]**2+Deviator[2]**2) + (Deviator[3]**2+Deviator[4]**2+Deviator[5]**2)
VonMisesStress = sqrt(3.0*J2)



# Assuming Von Mises stress!
A = 1.0 / (Gf * Young / (characteristic_length * threshold**2) - 0.5)




# only for exponential softening!
damage = 1.0 - (threshold / VonMisesStress) * exp(A * (1.0 - VonMisesStress / threshold))

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


# out = OutputVector_CollectingFactors(Seff, "stress", mode)
# print(out)

# Output:
'''
const double cr_Ct0 = nu - 1.0;
const double cr_Ct1 = -nu;
const double cr_Ct2 = std::pow(cr_Ct1 + 0.5, -2);
const double cr_Ct3 = nu*r_strain[0];
const double cr_Ct4 = cr_Ct3;
const double cr_Ct5 = nu*r_strain[1];
const double cr_Ct6 = 0.5*cr_Ct5;
const double cr_Ct7 = -cr_Ct6;
const double cr_Ct8 = nu*r_strain[2];
const double cr_Ct9 = 0.5*cr_Ct8;
const double cr_Ct10 = -cr_Ct9;
const double cr_Ct11 = cr_Ct1 + 1.0;
const double cr_Ct12 = cr_Ct11*r_strain[0];
const double cr_Ct13 = cr_Ct11*r_strain[1];
const double cr_Ct14 = 0.5*cr_Ct13;
const double cr_Ct15 = cr_Ct11*r_strain[2];
const double cr_Ct16 = 0.5*cr_Ct15;
const double cr_Ct17 = cr_Ct10 - cr_Ct12 + cr_Ct14 + cr_Ct16 + cr_Ct4 + cr_Ct7;
const double cr_Ct18 = cr_Ct0*r_strain[1];
const double cr_Ct19 = cr_Ct5;
const double cr_Ct20 = cr_Ct0*r_strain[2];
const double cr_Ct21 = 0.5*cr_Ct20;
const double cr_Ct22 = 0.5*cr_Ct3;
const double cr_Ct23 = cr_Ct0*r_strain[0];
const double cr_Ct24 = 0.5*cr_Ct23;
const double cr_Ct25 = -cr_Ct22 - cr_Ct24;
const double cr_Ct26 = cr_Ct10 + cr_Ct18 + cr_Ct19 - cr_Ct21 + cr_Ct25;
const double cr_Ct27 = std::pow(nu - 0.5, -2);
const double cr_Ct28 = 0.22222222222222224*cr_Ct27;
const double cr_Ct29 = cr_Ct8;
const double cr_Ct30 = 0.5*cr_Ct18;
const double cr_Ct31 = cr_Ct20 + cr_Ct25 + cr_Ct29 - cr_Ct30 + cr_Ct7;
const double cr_Ct32 = std::pow(r_strain[3], 2);
const double cr_Ct33 = std::pow(r_strain[4], 2);
const double cr_Ct34 = std::pow(r_strain[5], 2);
const double cr_Ct35 = cr_Ct32 + cr_Ct33 + cr_Ct34;
const double cr_Ct36 = 0.22222222222222224*std::pow(cr_Ct17, 2)*cr_Ct2 + std::pow(cr_Ct26, 2)*cr_Ct28 + cr_Ct28*std::pow(cr_Ct31, 2) + cr_Ct35;
const double cr_Ct37 = nu + 1.0;
const double cr_Ct38 = std::pow(Young, 2)/std::pow(cr_Ct37, 2);
const double cr_Ct39 = cr_Ct36*cr_Ct38;
const double cr_Ct40 = std::sqrt(cr_Ct39);
const double cr_Ct41 = threshold/cr_Ct40;
const double cr_Ct42 = 1.1547005383792517*cr_Ct41;
const double cr_Ct43 = 2.0*nu;
const double cr_Ct44 = cr_Ct43 - 1.0;
const double cr_Ct45 = cr_Ct27*cr_Ct44;
const double cr_Ct46 = 0.25*cr_Ct45;
const double cr_Ct47 = cr_Ct31*cr_Ct46;
const double cr_Ct48 = -cr_Ct47;
const double cr_Ct49 = cr_Ct26*cr_Ct46;
const double cr_Ct50 = -cr_Ct49;
const double cr_Ct51 = 4.0*nu;
const double cr_Ct52 = cr_Ct51 - 2.0;
const double cr_Ct53 = cr_Ct17*cr_Ct2*cr_Ct52;
const double cr_Ct54 = cr_Ct48 + cr_Ct50 + 0.25*cr_Ct53;
const double cr_Ct55 = -cr_Ct5;
const double cr_Ct56 = -cr_Ct8;
const double cr_Ct57 = cr_Ct23 + cr_Ct55 + cr_Ct56;
const double cr_Ct58 = 1.0/(Gf*Young/(characteristic_length*std::pow(threshold, 2)) - 0.5);
const double cr_Ct59 = 1.0/cr_Ct36;
const double cr_Ct60 = cr_Ct58*cr_Ct59;
const double cr_Ct61 = 0.44444444444444448*cr_Ct60;
const double cr_Ct62 = cr_Ct57*cr_Ct61;
const double cr_Ct63 = 2.0 - cr_Ct51;
const double cr_Ct64 = -cr_Ct4;
const double cr_Ct65 = -cr_Ct16 + cr_Ct9;
const double cr_Ct66 = -cr_Ct14 + cr_Ct6;
const double cr_Ct67 = cr_Ct12 + cr_Ct64 + cr_Ct65 + cr_Ct66;
const double cr_Ct68 = -cr_Ct19;
const double cr_Ct69 = -0.5*cr_Ct12 + cr_Ct22;
const double cr_Ct70 = cr_Ct13 + cr_Ct65 + cr_Ct68 + cr_Ct69;
const double cr_Ct71 = cr_Ct44*cr_Ct70;
const double cr_Ct72 = -cr_Ct29;
const double cr_Ct73 = cr_Ct15 + cr_Ct66 + cr_Ct69 + cr_Ct72;
const double cr_Ct74 = cr_Ct44*cr_Ct73;
const double cr_Ct75 = cr_Ct63*cr_Ct67 + cr_Ct71 + cr_Ct74;
const double cr_Ct76 = cr_Ct38/std::pow(cr_Ct39, 3.0/2.0);
const double cr_Ct77 = cr_Ct76*threshold;
const double cr_Ct78 = cr_Ct2*cr_Ct77;
const double cr_Ct79 = 0.12830005981991685*cr_Ct78;
const double cr_Ct80 = cr_Ct57*cr_Ct79;
const double cr_Ct81 = 1.0/(cr_Ct43 - 1.0);
const double cr_Ct82 = 0.8660254037844386/threshold;
const double cr_Ct83 = cr_Ct58;
const double cr_Ct84 = Young/cr_Ct37;
const double cr_Ct85 = cr_Ct84*std::exp(cr_Ct83*(-cr_Ct40*cr_Ct82 + 1.0));
const double cr_Ct86 = cr_Ct81*cr_Ct85;
const double cr_Ct87 = cr_Ct42*nu;
const double cr_Ct88 = cr_Ct44*cr_Ct67;
const double cr_Ct89 = cr_Ct63*cr_Ct70 + cr_Ct74 + cr_Ct88;
const double cr_Ct90 = cr_Ct81/(1.0 - cr_Ct43);
const double cr_Ct91 = cr_Ct52*cr_Ct90;
const double cr_Ct92 = cr_Ct26*cr_Ct91;
const double cr_Ct93 = cr_Ct17*cr_Ct44*cr_Ct90;
const double cr_Ct94 = cr_Ct48 - cr_Ct92 + cr_Ct93;
const double cr_Ct95 = cr_Ct63*cr_Ct73 + cr_Ct71 + cr_Ct88;
const double cr_Ct96 = cr_Ct31*cr_Ct91;
const double cr_Ct97 = cr_Ct50 + cr_Ct93 - cr_Ct96;
const double cr_Ct98 = 1.1547005383792517*threshold;
const double cr_Ct99 = cr_Ct86*(cr_Ct59*cr_Ct83 + cr_Ct76*cr_Ct98);
const double cr_Ct100 = cr_Ct57*cr_Ct99;
const double cr_Ct101 = -cr_Ct3;
const double cr_Ct102 = cr_Ct101 + cr_Ct18 + cr_Ct56;
const double cr_Ct103 = cr_Ct102*cr_Ct61;
const double cr_Ct104 = cr_Ct102*cr_Ct79;
const double cr_Ct105 = cr_Ct21 + cr_Ct9;
const double cr_Ct106 = cr_Ct30 + cr_Ct6;
const double cr_Ct107 = -cr_Ct18;
const double cr_Ct108 = cr_Ct22 + cr_Ct24;
const double cr_Ct109 = -cr_Ct20;
const double cr_Ct110 = cr_Ct28*std::pow(cr_Ct105 + cr_Ct106 - cr_Ct23 + cr_Ct64, 2) + cr_Ct28*std::pow(cr_Ct105 + cr_Ct107 + cr_Ct108 + cr_Ct68, 2) + cr_Ct28*std::pow(cr_Ct106 + cr_Ct108 + cr_Ct109 + cr_Ct72, 2) + cr_Ct35;
const double cr_Ct111 = cr_Ct110*cr_Ct38;
const double cr_Ct112 = std::sqrt(cr_Ct111);
const double cr_Ct113 = cr_Ct0*cr_Ct98/cr_Ct112;
const double cr_Ct114 = cr_Ct107 + cr_Ct3 + cr_Ct8;
const double cr_Ct115 = -cr_Ct93;
const double cr_Ct116 = 0.44444444444444448*cr_Ct58/cr_Ct110;
const double cr_Ct117 = 0.5132002392796674*cr_Ct38*threshold/std::pow(cr_Ct111, 3.0/2.0);
const double cr_Ct118 = cr_Ct81*cr_Ct84*std::exp(cr_Ct83*(-cr_Ct112*cr_Ct82 + 1.0));
const double cr_Ct119 = cr_Ct102*cr_Ct99;
const double cr_Ct120 = cr_Ct101 + cr_Ct20 + cr_Ct55;
const double cr_Ct121 = cr_Ct120*cr_Ct61;
const double cr_Ct122 = cr_Ct120*cr_Ct79;
const double cr_Ct123 = cr_Ct109 + cr_Ct3 + cr_Ct5;
const double cr_Ct124 = cr_Ct120*cr_Ct99;
const double cr_Ct125 = 0.055555555555555559*cr_Ct45;
const double cr_Ct126 = -cr_Ct125*cr_Ct31;
const double cr_Ct127 = -cr_Ct125*cr_Ct26;
const double cr_Ct128 = 0.064150029909958425*cr_Ct78;
const double cr_Ct129 = cr_Ct128*cr_Ct75 + cr_Ct60*(cr_Ct126 + cr_Ct127 + 0.055555555555555559*cr_Ct53);
const double cr_Ct130 = cr_Ct85*r_strain[3];
const double cr_Ct131 = 0.22222222222222224*cr_Ct93;
const double cr_Ct132 = cr_Ct128*cr_Ct89 + cr_Ct60*(cr_Ct126 + cr_Ct131 - 0.22222222222222224*cr_Ct92);
const double cr_Ct133 = cr_Ct128*cr_Ct95 + cr_Ct60*(cr_Ct127 + cr_Ct131 - 0.22222222222222224*cr_Ct96);
const double cr_Ct134 = 0.57735026918962584*cr_Ct41;
const double cr_Ct135 = 0.5*cr_Ct60;
const double cr_Ct136 = 0.57735026918962584*cr_Ct77;
const double cr_Ct137 = cr_Ct135 + cr_Ct136;
const double cr_Ct138 = cr_Ct130*cr_Ct137;
const double cr_Ct139 = -cr_Ct138*r_strain[4];
const double cr_Ct140 = -cr_Ct138*r_strain[5];
const double cr_Ct141 = cr_Ct85*r_strain[4];
const double cr_Ct142 = -cr_Ct137*cr_Ct141*r_strain[5];
const double cr_Ct143 = cr_Ct85*r_strain[5];
r_Ct(0,0)=cr_Ct86*(cr_Ct0*cr_Ct42 - cr_Ct54*cr_Ct62 - cr_Ct75*cr_Ct80);
r_Ct(0,1)=-cr_Ct86*(cr_Ct62*cr_Ct94 + cr_Ct80*cr_Ct89 + cr_Ct87);
r_Ct(0,2)=-cr_Ct86*(cr_Ct62*cr_Ct97 + cr_Ct80*cr_Ct95 + cr_Ct87);
r_Ct(0,3)=-cr_Ct100*r_strain[3];
r_Ct(0,4)=-cr_Ct100*r_strain[4];
r_Ct(0,5)=-cr_Ct100*r_strain[5];
r_Ct(1,0)=-cr_Ct86*(cr_Ct103*cr_Ct54 + cr_Ct104*cr_Ct75 + cr_Ct87);
r_Ct(1,1)=cr_Ct118*(cr_Ct113 - cr_Ct114*cr_Ct116*(cr_Ct115 + cr_Ct47 + cr_Ct92) + cr_Ct114*cr_Ct117*cr_Ct94);
r_Ct(1,2)=-cr_Ct86*(cr_Ct103*cr_Ct97 + cr_Ct104*cr_Ct95 + cr_Ct87);
r_Ct(1,3)=-cr_Ct119*r_strain[3];
r_Ct(1,4)=-cr_Ct119*r_strain[4];
r_Ct(1,5)=-cr_Ct119*r_strain[5];
r_Ct(2,0)=-cr_Ct86*(cr_Ct121*cr_Ct54 + cr_Ct122*cr_Ct75 + cr_Ct87);
r_Ct(2,1)=-cr_Ct86*(cr_Ct121*cr_Ct94 + cr_Ct122*cr_Ct89 + cr_Ct87);
r_Ct(2,2)=cr_Ct118*(cr_Ct113 - cr_Ct116*cr_Ct123*(cr_Ct115 + cr_Ct49 + cr_Ct96) + cr_Ct117*cr_Ct123*cr_Ct97);
r_Ct(2,3)=-cr_Ct124*r_strain[3];
r_Ct(2,4)=-cr_Ct124*r_strain[4];
r_Ct(2,5)=-cr_Ct124*r_strain[5];
r_Ct(3,0)=-cr_Ct129*cr_Ct130;
r_Ct(3,1)=-cr_Ct130*cr_Ct132;
r_Ct(3,2)=-cr_Ct130*cr_Ct133;
r_Ct(3,3)=cr_Ct85*(cr_Ct134 - cr_Ct135*cr_Ct32 - cr_Ct136*cr_Ct32);
r_Ct(3,4)=cr_Ct139;
r_Ct(3,5)=cr_Ct140;
r_Ct(4,0)=-cr_Ct129*cr_Ct141;
r_Ct(4,1)=-cr_Ct132*cr_Ct141;
r_Ct(4,2)=-cr_Ct133*cr_Ct141;
r_Ct(4,3)=cr_Ct139;
r_Ct(4,4)=cr_Ct85*(cr_Ct134 - cr_Ct135*cr_Ct33 - cr_Ct136*cr_Ct33);
r_Ct(4,5)=cr_Ct142;
r_Ct(5,0)=-cr_Ct129*cr_Ct143;
r_Ct(5,1)=-cr_Ct132*cr_Ct143;
r_Ct(5,2)=-cr_Ct133*cr_Ct143;
r_Ct(5,3)=cr_Ct140;
r_Ct(5,4)=cr_Ct142;
r_Ct(5,5)=cr_Ct85*(cr_Ct134 - cr_Ct135*cr_Ct34 - cr_Ct136*cr_Ct34);

'''