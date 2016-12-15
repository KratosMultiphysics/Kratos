// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_CONTACT2D2N2NDLM)
#define KRATOS_CONTACT2D2N2NDLM

// System includes

// External includes

// Project includes
#include "custom_conditions/mortar_contact_condition.h"
#include "custom_utilities/contact_utilities.h"
#include "structural_mechanics_application.h"
#include "structural_mechanics_application_variables.h"
#include <boost/math/special_functions/sign.hpp>

namespace Kratos 
{
    ///@name Kratos Globals
    ///@{

    ///@}
    ///@name Type Definitions
    ///@{
        
    typedef Point<3>                                  PointType;
    typedef Node<3>                                    NodeType;
    typedef Geometry<NodeType>                     GeometryType;
        
    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name  Functions
    ///@{

    ///@}
    ///@name Kratos Classes
    ///@{
        
class Contact2D2N2NDLM
{
public:
    
    static inline bounded_matrix<double,16,16> ComputeGaussPointActiveLHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    bounded_matrix<double,16,16> lhs;
    
    const Matrix normalmaster     = rContactData.Normal_m;
    const Matrix normalslave      = rContactData.Normal_s;
    const Matrix tan1slave        = rContactData.Tangent_xi_s;
    const Matrix lm               = rContactData.LagrangeMultipliers;
    const Matrix dlm              = rContactData.DoubleLagrangeMultipliers;
//     const double Dt               = rContactData.Dt;
    const double epsilon          = rContactData.epsilon;
//     const double epsilon_tangent  = rContactData.epsilon_tangent;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const Matrix v1 = rContactData.v1;
    const Matrix v2 = rContactData.v2;

    const std::vector<double> DeltaJs         = rContactData.DeltaJ_s;
    const std::vector<Matrix> DeltaNormals    = rContactData.Delta_Normal_s;
    const std::vector<Matrix> Deltatangentxis = rContactData.Delta_Tangent_xi_s;
    const std::vector<Matrix> DeltaNormalm    = rContactData.Delta_Normal_m;
    const std::vector<double> DeltaGap        = rContactData.DeltaGap;
    const std::vector<Vector> DeltaPhi        = rContactData.DeltaPhi;
    const std::vector<Vector> DeltaN2         = rContactData.DeltaN2;
    
    const double Dnormalslave11u111 =     DeltaNormals[3](1,1);
    const double Dnormalslave11u110 =     DeltaNormals[2](1,1);
    const double Dnormalslave11u101 =     DeltaNormals[1](1,1);
    const double Dnormalslave11u100 =     DeltaNormals[0](1,1);
    const double Dnormalslave10u111 =     DeltaNormals[3](1,0);
    const double Dnormalslave10u110 =     DeltaNormals[2](1,0);
    const double Dnormalslave10u101 =     DeltaNormals[1](1,0);
    const double Dnormalslave10u100 =     DeltaNormals[0](1,0);
    const double Dnormalslave01u111 =     DeltaNormals[3](0,1);
    const double Dnormalslave01u110 =     DeltaNormals[2](0,1);
    const double Dnormalslave01u101 =     DeltaNormals[1](0,1);
    const double Dnormalslave01u100 =     DeltaNormals[0](0,1);
    const double Dnormalslave00u111 =     DeltaNormals[3](0,0);
    const double Dnormalslave00u110 =     DeltaNormals[2](0,0);
    const double Dnormalslave00u101 =     DeltaNormals[1](0,0);
    const double Dnormalslave00u100 =     DeltaNormals[0](0,0);
    const double DPhi1u111 =     DeltaPhi[3][1];
    const double DPhi1u110 =     DeltaPhi[2][1];
    const double DPhi1u101 =     DeltaPhi[1][1];
    const double DPhi1u100 =     DeltaPhi[0][1];
    const double DPhi0u111 =     DeltaPhi[3][0];
    const double DPhi0u110 =     DeltaPhi[2][0];
    const double DPhi0u101 =     DeltaPhi[1][0];
    const double DPhi0u100 =     DeltaPhi[0][0];
    const double DN21u211 =     DeltaN2[7][1];
    const double DN21u210 =     DeltaN2[6][1];
    const double DN21u201 =     DeltaN2[5][1];
    const double DN21u200 =     DeltaN2[4][1];
    const double DN21u111 =     DeltaN2[3][1];
    const double DN21u110 =     DeltaN2[2][1];
    const double DN21u101 =     DeltaN2[1][1];
    const double DN21u100 =     DeltaN2[0][1];
    const double DN20u211 =     DeltaN2[7][0];
    const double DN20u210 =     DeltaN2[6][0];
    const double DN20u201 =     DeltaN2[5][0];
    const double DN20u200 =     DeltaN2[4][0];
    const double DN20u111 =     DeltaN2[3][0];
    const double DN20u110 =     DeltaN2[2][0];
    const double DN20u101 =     DeltaN2[1][0];
    const double DN20u100 =     DeltaN2[0][0];
    const double Dgapgu211 =     DeltaGap[7];
    const double Dgapgu210 =     DeltaGap[6];
    const double Dgapgu201 =     DeltaGap[5];
    const double Dgapgu200 =     DeltaGap[4];
    const double Dgapgu111 =     DeltaGap[3];
    const double Dgapgu110 =     DeltaGap[2];
    const double Dgapgu101 =     DeltaGap[1];
    const double Dgapgu100 =     DeltaGap[0];
    const double DdetJu111 =     DeltaJs[3];
    const double DdetJu110 =     DeltaJs[2];
    const double DdetJu101 =     DeltaJs[1];
    const double DdetJu100 =     DeltaJs[0];
 
    const double clhs0 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs1 =     DN20u200; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs2 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs3 =     Phi[0]; // PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs4 =     clhs3*dlm(0,0);
    const double clhs5 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs6 =     clhs5*dlm(1,0);
    const double clhs7 =     clhs4 + clhs6;
    const double clhs8 =     clhs3*lm(0,0);
    const double clhs9 =     clhs5*lm(1,0);
    const double clhs10 =     clhs8 + clhs9;
    const double clhs11 =     clhs2*(clhs10 + clhs7);
    const double clhs12 =     DN20u201; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs13 =     DN20u210; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs14 =     DN20u211; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs15 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs16 =     clhs15*clhs7;
    const double clhs17 =     DN20u100; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs18 =     clhs2*clhs7;
    const double clhs19 =     clhs10*clhs15;
    const double clhs20 =     clhs10*clhs2;
    const double clhs21 =     clhs0*clhs2;
    const double clhs22 =     DPhi0u100; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs23 =     clhs22*dlm(0,0);
    const double clhs24 =     DPhi1u100; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs25 =     clhs24*dlm(1,0);
    const double clhs26 =     clhs23 + clhs25;
    const double clhs27 =     clhs22*lm(0,0);
    const double clhs28 =     clhs24*lm(1,0);
    const double clhs29 =     clhs27 + clhs28;
    const double clhs30 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs31 =     clhs30*clhs7;
    const double clhs32 =     DN20u101; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs33 =     clhs10*clhs30;
    const double clhs34 =     DPhi0u101; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs35 =     clhs34*dlm(0,0);
    const double clhs36 =     DPhi1u101; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs37 =     clhs36*dlm(1,0);
    const double clhs38 =     clhs35 + clhs37;
    const double clhs39 =     clhs34*lm(0,0);
    const double clhs40 =     clhs36*lm(1,0);
    const double clhs41 =     clhs39 + clhs40;
    const double clhs42 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs43 =     clhs42*clhs7;
    const double clhs44 =     DN20u110; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs45 =     clhs10*clhs42;
    const double clhs46 =     DPhi0u110; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs47 =     clhs46*dlm(0,0);
    const double clhs48 =     DPhi1u110; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs49 =     clhs48*dlm(1,0);
    const double clhs50 =     clhs47 + clhs49;
    const double clhs51 =     clhs46*lm(0,0);
    const double clhs52 =     clhs48*lm(1,0);
    const double clhs53 =     clhs51 + clhs52;
    const double clhs54 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs55 =     clhs54*clhs7;
    const double clhs56 =     DN20u111; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs57 =     clhs10*clhs54;
    const double clhs58 =     DPhi0u111; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs59 =     clhs58*dlm(0,0);
    const double clhs60 =     DPhi1u111; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs61 =     clhs60*dlm(1,0);
    const double clhs62 =     clhs59 + clhs61;
    const double clhs63 =     clhs58*lm(0,0);
    const double clhs64 =     clhs60*lm(1,0);
    const double clhs65 =     clhs63 + clhs64;
    const double clhs66 =     clhs21*clhs3;
    const double clhs67 =     clhs21*clhs5;
    const double clhs68 =     clhs3*dlm(0,1);
    const double clhs69 =     clhs5*dlm(1,1);
    const double clhs70 =     clhs68 + clhs69;
    const double clhs71 =     clhs3*lm(0,1);
    const double clhs72 =     clhs5*lm(1,1);
    const double clhs73 =     clhs71 + clhs72;
    const double clhs74 =     clhs2*(clhs70 + clhs73);
    const double clhs75 =     clhs15*clhs70;
    const double clhs76 =     clhs2*clhs70;
    const double clhs77 =     clhs15*clhs73;
    const double clhs78 =     clhs2*clhs73;
    const double clhs79 =     clhs22*dlm(0,1);
    const double clhs80 =     clhs24*dlm(1,1);
    const double clhs81 =     clhs79 + clhs80;
    const double clhs82 =     clhs22*lm(0,1);
    const double clhs83 =     clhs24*lm(1,1);
    const double clhs84 =     clhs82 + clhs83;
    const double clhs85 =     clhs30*clhs70;
    const double clhs86 =     clhs30*clhs73;
    const double clhs87 =     clhs34*dlm(0,1);
    const double clhs88 =     clhs36*dlm(1,1);
    const double clhs89 =     clhs87 + clhs88;
    const double clhs90 =     clhs34*lm(0,1);
    const double clhs91 =     clhs36*lm(1,1);
    const double clhs92 =     clhs90 + clhs91;
    const double clhs93 =     clhs42*clhs70;
    const double clhs94 =     clhs42*clhs73;
    const double clhs95 =     clhs46*dlm(0,1);
    const double clhs96 =     clhs48*dlm(1,1);
    const double clhs97 =     clhs95 + clhs96;
    const double clhs98 =     clhs46*lm(0,1);
    const double clhs99 =     clhs48*lm(1,1);
    const double clhs100 =     clhs98 + clhs99;
    const double clhs101 =     clhs54*clhs70;
    const double clhs102 =     clhs54*clhs73;
    const double clhs103 =     clhs58*dlm(0,1);
    const double clhs104 =     clhs60*dlm(1,1);
    const double clhs105 =     clhs103 + clhs104;
    const double clhs106 =     clhs58*lm(0,1);
    const double clhs107 =     clhs60*lm(1,1);
    const double clhs108 =     clhs106 + clhs107;
    const double clhs109 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs110 =     DN21u200; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs111 =     DN21u201; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs112 =     DN21u210; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs113 =     DN21u211; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs114 =     DN21u100; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs115 =     clhs109*clhs2;
    const double clhs116 =     DN21u101; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs117 =     DN21u110; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs118 =     DN21u111; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs119 =     clhs115*clhs3;
    const double clhs120 =     clhs115*clhs5;
    const double clhs121 =     clhs16 + clhs19 + clhs2*clhs26 + clhs2*clhs29;
    const double clhs122 =     clhs2*clhs38 + clhs2*clhs41 + clhs31 + clhs33;
    const double clhs123 =     clhs2*clhs50 + clhs2*clhs53 + clhs43 + clhs45;
    const double clhs124 =     clhs2*clhs62 + clhs2*clhs65 + clhs55 + clhs57;
    const double clhs125 =     N1[0]*clhs2;
    const double clhs126 =     -clhs125*clhs3;
    const double clhs127 =     -clhs125*clhs5;
    const double clhs128 =     clhs2*clhs81 + clhs2*clhs84 + clhs75 + clhs77;
    const double clhs129 =     clhs2*clhs89 + clhs2*clhs92 + clhs85 + clhs86;
    const double clhs130 =     clhs100*clhs2 + clhs2*clhs97 + clhs93 + clhs94;
    const double clhs131 =     clhs101 + clhs102 + clhs105*clhs2 + clhs108*clhs2;
    const double clhs132 =     N1[1]*clhs2;
    const double clhs133 =     -clhs132*clhs3;
    const double clhs134 =     -clhs132*clhs5;
    const double clhs135 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs136 =     Dgapgu200; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs137 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs138 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs139 =     clhs2*clhs3*(GPnormal[0]*clhs137 + GPnormal[1]*clhs138);
    const double clhs140 =     -clhs136*clhs139;
    const double clhs141 =     Dgapgu201; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs142 =     -clhs139*clhs141;
    const double clhs143 =     Dgapgu210; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs144 =     -clhs139*clhs143;
    const double clhs145 =     Dgapgu211; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs146 =     -clhs139*clhs145;
    const double clhs147 =     clhs135*clhs2*clhs3;
    const double clhs148 =     Dgapgu100; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs149 =     clhs137*clhs2*clhs3;
    const double clhs150 =     clhs135*clhs137*clhs3;
    const double clhs151 =     clhs135*clhs137*clhs2;
    const double clhs152 =     clhs147*Dnormalslave00u100 + clhs148*clhs149 + clhs15*clhs150 + clhs151*clhs22; // CLHS147*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS148*CLHS149 + CLHS15*CLHS150 + CLHS151*CLHS22
    const double clhs153 =     clhs10 - clhs4 - clhs6;
    const double clhs154 =     clhs153*clhs3*epsilon;
    const double clhs155 =     clhs153*clhs2*epsilon;
    const double clhs156 =     clhs2*epsilon*(-clhs23 - clhs25 + clhs29);
    const double clhs157 =     clhs15*clhs154 + clhs152 + clhs155*clhs22 + clhs156*clhs3;
    const double clhs158 =     clhs138*clhs2*clhs3;
    const double clhs159 =     clhs135*clhs138*clhs3;
    const double clhs160 =     clhs135*clhs138*clhs2;
    const double clhs161 =     clhs147*Dnormalslave01u100 + clhs148*clhs158 + clhs15*clhs159 + clhs160*clhs22; // CLHS147*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS148*CLHS158 + CLHS15*CLHS159 + CLHS160*CLHS22
    const double clhs162 =     -clhs68 - clhs69 + clhs73;
    const double clhs163 =     clhs162*clhs3*epsilon;
    const double clhs164 =     clhs162*clhs2*epsilon;
    const double clhs165 =     clhs2*epsilon*(-clhs79 - clhs80 + clhs84);
    const double clhs166 =     clhs15*clhs163 + clhs161 + clhs164*clhs22 + clhs165*clhs3;
    const double clhs167 =     Dgapgu101; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs168 =     clhs147*Dnormalslave00u101 + clhs149*clhs167 + clhs150*clhs30 + clhs151*clhs34; // CLHS147*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS149*CLHS167 + CLHS150*CLHS30 + CLHS151*CLHS34
    const double clhs169 =     clhs2*epsilon*(-clhs35 - clhs37 + clhs41);
    const double clhs170 =     clhs154*clhs30 + clhs155*clhs34 + clhs168 + clhs169*clhs3;
    const double clhs171 =     clhs147*Dnormalslave01u101 + clhs158*clhs167 + clhs159*clhs30 + clhs160*clhs34; // CLHS147*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS158*CLHS167 + CLHS159*CLHS30 + CLHS160*CLHS34
    const double clhs172 =     clhs2*epsilon*(-clhs87 - clhs88 + clhs92);
    const double clhs173 =     clhs163*clhs30 + clhs164*clhs34 + clhs171 + clhs172*clhs3;
    const double clhs174 =     Dgapgu110; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs175 =     clhs147*Dnormalslave00u110 + clhs149*clhs174 + clhs150*clhs42 + clhs151*clhs46; // CLHS147*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS149*CLHS174 + CLHS150*CLHS42 + CLHS151*CLHS46
    const double clhs176 =     clhs2*epsilon*(-clhs47 - clhs49 + clhs53);
    const double clhs177 =     clhs154*clhs42 + clhs155*clhs46 + clhs175 + clhs176*clhs3;
    const double clhs178 =     clhs147*Dnormalslave01u110 + clhs158*clhs174 + clhs159*clhs42 + clhs160*clhs46; // CLHS147*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS158*CLHS174 + CLHS159*CLHS42 + CLHS160*CLHS46
    const double clhs179 =     clhs2*epsilon*(clhs100 - clhs95 - clhs96);
    const double clhs180 =     clhs163*clhs42 + clhs164*clhs46 + clhs178 + clhs179*clhs3;
    const double clhs181 =     Dgapgu111; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs182 =     clhs147*Dnormalslave00u111 + clhs149*clhs181 + clhs150*clhs54 + clhs151*clhs58; // CLHS147*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS149*CLHS181 + CLHS150*CLHS54 + CLHS151*CLHS58
    const double clhs183 =     clhs2*epsilon*(-clhs59 - clhs61 + clhs65);
    const double clhs184 =     clhs154*clhs54 + clhs155*clhs58 + clhs182 + clhs183*clhs3;
    const double clhs185 =     clhs147*Dnormalslave01u111 + clhs158*clhs181 + clhs159*clhs54 + clhs160*clhs58; // CLHS147*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS158*CLHS181 + CLHS159*CLHS54 + CLHS160*CLHS58
    const double clhs186 =     clhs2*epsilon*(-clhs103 - clhs104 + clhs108);
    const double clhs187 =     clhs163*clhs54 + clhs164*clhs58 + clhs185 + clhs186*clhs3;
    const double clhs188 =     std::pow(clhs3, 2);
    const double clhs189 =     GPnormal[0]*clhs2*epsilon;
    const double clhs190 =     clhs188*clhs189;
    const double clhs191 =     -clhs190;
    const double clhs192 =     GPnormal[1]*clhs2*epsilon;
    const double clhs193 =     clhs188*clhs192;
    const double clhs194 =     -clhs193;
    const double clhs195 =     clhs3*clhs5;
    const double clhs196 =     clhs189*clhs195;
    const double clhs197 =     -clhs196;
    const double clhs198 =     clhs192*clhs195;
    const double clhs199 =     -clhs198;
    const double clhs200 =     clhs2*clhs3*(GPtangent1[0]*clhs137 + GPtangent1[1]*clhs138);
    const double clhs201 =     -clhs136*clhs200;
    const double clhs202 =     -clhs141*clhs200;
    const double clhs203 =     -clhs143*clhs200;
    const double clhs204 =     -clhs145*clhs200;
    const double clhs205 =     GPtangent1[0]*clhs2*epsilon;
    const double clhs206 =     clhs188*clhs205;
    const double clhs207 =     -clhs206;
    const double clhs208 =     GPtangent1[1]*clhs2*epsilon;
    const double clhs209 =     clhs188*clhs208;
    const double clhs210 =     -clhs209;
    const double clhs211 =     clhs195*clhs205;
    const double clhs212 =     -clhs211;
    const double clhs213 =     clhs195*clhs208;
    const double clhs214 =     -clhs213;
    const double clhs215 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs216 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs217 =     clhs2*clhs5*(GPnormal[0]*clhs215 + GPnormal[1]*clhs216);
    const double clhs218 =     -clhs136*clhs217;
    const double clhs219 =     -clhs141*clhs217;
    const double clhs220 =     -clhs143*clhs217;
    const double clhs221 =     -clhs145*clhs217;
    const double clhs222 =     clhs135*clhs2*clhs5;
    const double clhs223 =     clhs2*clhs215*clhs5;
    const double clhs224 =     clhs135*clhs215*clhs5;
    const double clhs225 =     clhs135*clhs2*clhs215;
    const double clhs226 =     clhs148*clhs223 + clhs15*clhs224 + clhs222*Dnormalslave10u100 + clhs225*clhs24; // CLHS148*CLHS223 + CLHS15*CLHS224 + CLHS222*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS225*CLHS24
    const double clhs227 =     clhs153*clhs5*epsilon;
    const double clhs228 =     clhs15*clhs227 + clhs155*clhs24 + clhs156*clhs5 + clhs226;
    const double clhs229 =     clhs2*clhs216*clhs5;
    const double clhs230 =     clhs135*clhs216*clhs5;
    const double clhs231 =     clhs135*clhs2*clhs216;
    const double clhs232 =     clhs148*clhs229 + clhs15*clhs230 + clhs222*Dnormalslave11u100 + clhs231*clhs24; // CLHS148*CLHS229 + CLHS15*CLHS230 + CLHS222*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS231*CLHS24
    const double clhs233 =     clhs162*clhs5*epsilon;
    const double clhs234 =     clhs15*clhs233 + clhs164*clhs24 + clhs165*clhs5 + clhs232;
    const double clhs235 =     clhs167*clhs223 + clhs222*Dnormalslave10u101 + clhs224*clhs30 + clhs225*clhs36; // CLHS167*CLHS223 + CLHS222*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS224*CLHS30 + CLHS225*CLHS36
    const double clhs236 =     clhs155*clhs36 + clhs169*clhs5 + clhs227*clhs30 + clhs235;
    const double clhs237 =     clhs167*clhs229 + clhs222*Dnormalslave11u101 + clhs230*clhs30 + clhs231*clhs36; // CLHS167*CLHS229 + CLHS222*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS230*CLHS30 + CLHS231*CLHS36
    const double clhs238 =     clhs164*clhs36 + clhs172*clhs5 + clhs233*clhs30 + clhs237;
    const double clhs239 =     clhs174*clhs223 + clhs222*Dnormalslave10u110 + clhs224*clhs42 + clhs225*clhs48; // CLHS174*CLHS223 + CLHS222*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS224*CLHS42 + CLHS225*CLHS48
    const double clhs240 =     clhs155*clhs48 + clhs176*clhs5 + clhs227*clhs42 + clhs239;
    const double clhs241 =     clhs174*clhs229 + clhs222*Dnormalslave11u110 + clhs230*clhs42 + clhs231*clhs48; // CLHS174*CLHS229 + CLHS222*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS230*CLHS42 + CLHS231*CLHS48
    const double clhs242 =     clhs164*clhs48 + clhs179*clhs5 + clhs233*clhs42 + clhs241;
    const double clhs243 =     clhs181*clhs223 + clhs222*Dnormalslave10u111 + clhs224*clhs54 + clhs225*clhs60; // CLHS181*CLHS223 + CLHS222*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS224*CLHS54 + CLHS225*CLHS60
    const double clhs244 =     clhs155*clhs60 + clhs183*clhs5 + clhs227*clhs54 + clhs243;
    const double clhs245 =     clhs181*clhs229 + clhs222*Dnormalslave11u111 + clhs230*clhs54 + clhs231*clhs60; // CLHS181*CLHS229 + CLHS222*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS230*CLHS54 + CLHS231*CLHS60
    const double clhs246 =     clhs164*clhs60 + clhs186*clhs5 + clhs233*clhs54 + clhs245;
    const double clhs247 =     std::pow(clhs5, 2);
    const double clhs248 =     clhs189*clhs247;
    const double clhs249 =     -clhs248;
    const double clhs250 =     clhs192*clhs247;
    const double clhs251 =     -clhs250;
    const double clhs252 =     clhs2*clhs5*(GPtangent1[0]*clhs215 + GPtangent1[1]*clhs216);
    const double clhs253 =     -clhs136*clhs252;
    const double clhs254 =     -clhs141*clhs252;
    const double clhs255 =     -clhs143*clhs252;
    const double clhs256 =     -clhs145*clhs252;
    const double clhs257 =     clhs205*clhs247;
    const double clhs258 =     -clhs257;
    const double clhs259 =     clhs208*clhs247;
    const double clhs260 =     -clhs259;
    const double clhs261 =     clhs7 - clhs8 - clhs9;
    const double clhs262 =     clhs261*clhs3*epsilon;
    const double clhs263 =     clhs2*clhs261*epsilon;
    const double clhs264 =     clhs2*epsilon*(clhs26 - clhs27 - clhs28);
    const double clhs265 =     clhs15*clhs262 + clhs152 + clhs22*clhs263 + clhs264*clhs3;
    const double clhs266 =     clhs70 - clhs71 - clhs72;
    const double clhs267 =     clhs266*clhs3*epsilon;
    const double clhs268 =     clhs2*clhs266*epsilon;
    const double clhs269 =     clhs2*epsilon*(clhs81 - clhs82 - clhs83);
    const double clhs270 =     clhs15*clhs267 + clhs161 + clhs22*clhs268 + clhs269*clhs3;
    const double clhs271 =     clhs2*epsilon*(clhs38 - clhs39 - clhs40);
    const double clhs272 =     clhs168 + clhs262*clhs30 + clhs263*clhs34 + clhs271*clhs3;
    const double clhs273 =     clhs2*epsilon*(clhs89 - clhs90 - clhs91);
    const double clhs274 =     clhs171 + clhs267*clhs30 + clhs268*clhs34 + clhs273*clhs3;
    const double clhs275 =     clhs2*epsilon*(clhs50 - clhs51 - clhs52);
    const double clhs276 =     clhs175 + clhs262*clhs42 + clhs263*clhs46 + clhs275*clhs3;
    const double clhs277 =     clhs2*epsilon*(clhs97 - clhs98 - clhs99);
    const double clhs278 =     clhs178 + clhs267*clhs42 + clhs268*clhs46 + clhs277*clhs3;
    const double clhs279 =     clhs2*epsilon*(clhs62 - clhs63 - clhs64);
    const double clhs280 =     clhs182 + clhs262*clhs54 + clhs263*clhs58 + clhs279*clhs3;
    const double clhs281 =     clhs2*epsilon*(clhs105 - clhs106 - clhs107);
    const double clhs282 =     clhs185 + clhs267*clhs54 + clhs268*clhs58 + clhs281*clhs3;
    const double clhs283 =     clhs261*clhs5*epsilon;
    const double clhs284 =     clhs15*clhs283 + clhs226 + clhs24*clhs263 + clhs264*clhs5;
    const double clhs285 =     clhs266*clhs5*epsilon;
    const double clhs286 =     clhs15*clhs285 + clhs232 + clhs24*clhs268 + clhs269*clhs5;
    const double clhs287 =     clhs235 + clhs263*clhs36 + clhs271*clhs5 + clhs283*clhs30;
    const double clhs288 =     clhs237 + clhs268*clhs36 + clhs273*clhs5 + clhs285*clhs30;
    const double clhs289 =     clhs239 + clhs263*clhs48 + clhs275*clhs5 + clhs283*clhs42;
    const double clhs290 =     clhs241 + clhs268*clhs48 + clhs277*clhs5 + clhs285*clhs42;
    const double clhs291 =     clhs243 + clhs263*clhs60 + clhs279*clhs5 + clhs283*clhs54;
    const double clhs292 =     clhs245 + clhs268*clhs60 + clhs281*clhs5 + clhs285*clhs54;

    lhs(0,0)=clhs1*clhs11;
    lhs(0,1)=clhs11*clhs12;
    lhs(0,2)=clhs11*clhs13;
    lhs(0,3)=clhs11*clhs14;
    lhs(0,4)=clhs0*clhs16 + clhs0*clhs19 + clhs17*clhs18 + clhs17*clhs20 + clhs21*clhs26 + clhs21*clhs29;
    lhs(0,5)=clhs0*clhs31 + clhs0*clhs33 + clhs18*clhs32 + clhs20*clhs32 + clhs21*clhs38 + clhs21*clhs41;
    lhs(0,6)=clhs0*clhs43 + clhs0*clhs45 + clhs18*clhs44 + clhs20*clhs44 + clhs21*clhs50 + clhs21*clhs53;
    lhs(0,7)=clhs0*clhs55 + clhs0*clhs57 + clhs18*clhs56 + clhs20*clhs56 + clhs21*clhs62 + clhs21*clhs65;
    lhs(0,8)=clhs66;
    lhs(0,9)=0;
    lhs(0,10)=clhs67;
    lhs(0,11)=0;
    lhs(0,12)=clhs66;
    lhs(0,13)=0;
    lhs(0,14)=clhs67;
    lhs(0,15)=0;
    lhs(1,0)=clhs1*clhs74;
    lhs(1,1)=clhs12*clhs74;
    lhs(1,2)=clhs13*clhs74;
    lhs(1,3)=clhs14*clhs74;
    lhs(1,4)=clhs0*clhs75 + clhs0*clhs77 + clhs17*clhs76 + clhs17*clhs78 + clhs21*clhs81 + clhs21*clhs84;
    lhs(1,5)=clhs0*clhs85 + clhs0*clhs86 + clhs21*clhs89 + clhs21*clhs92 + clhs32*clhs76 + clhs32*clhs78;
    lhs(1,6)=clhs0*clhs93 + clhs0*clhs94 + clhs100*clhs21 + clhs21*clhs97 + clhs44*clhs76 + clhs44*clhs78;
    lhs(1,7)=clhs0*clhs101 + clhs0*clhs102 + clhs105*clhs21 + clhs108*clhs21 + clhs56*clhs76 + clhs56*clhs78;
    lhs(1,8)=0;
    lhs(1,9)=clhs66;
    lhs(1,10)=0;
    lhs(1,11)=clhs67;
    lhs(1,12)=0;
    lhs(1,13)=clhs66;
    lhs(1,14)=0;
    lhs(1,15)=clhs67;
    lhs(2,0)=clhs11*clhs110;
    lhs(2,1)=clhs11*clhs111;
    lhs(2,2)=clhs11*clhs112;
    lhs(2,3)=clhs11*clhs113;
    lhs(2,4)=clhs109*clhs16 + clhs109*clhs19 + clhs114*clhs18 + clhs114*clhs20 + clhs115*clhs26 + clhs115*clhs29;
    lhs(2,5)=clhs109*clhs31 + clhs109*clhs33 + clhs115*clhs38 + clhs115*clhs41 + clhs116*clhs18 + clhs116*clhs20;
    lhs(2,6)=clhs109*clhs43 + clhs109*clhs45 + clhs115*clhs50 + clhs115*clhs53 + clhs117*clhs18 + clhs117*clhs20;
    lhs(2,7)=clhs109*clhs55 + clhs109*clhs57 + clhs115*clhs62 + clhs115*clhs65 + clhs118*clhs18 + clhs118*clhs20;
    lhs(2,8)=clhs119;
    lhs(2,9)=0;
    lhs(2,10)=clhs120;
    lhs(2,11)=0;
    lhs(2,12)=clhs119;
    lhs(2,13)=0;
    lhs(2,14)=clhs120;
    lhs(2,15)=0;
    lhs(3,0)=clhs110*clhs74;
    lhs(3,1)=clhs111*clhs74;
    lhs(3,2)=clhs112*clhs74;
    lhs(3,3)=clhs113*clhs74;
    lhs(3,4)=clhs109*clhs75 + clhs109*clhs77 + clhs114*clhs76 + clhs114*clhs78 + clhs115*clhs81 + clhs115*clhs84;
    lhs(3,5)=clhs109*clhs85 + clhs109*clhs86 + clhs115*clhs89 + clhs115*clhs92 + clhs116*clhs76 + clhs116*clhs78;
    lhs(3,6)=clhs100*clhs115 + clhs109*clhs93 + clhs109*clhs94 + clhs115*clhs97 + clhs117*clhs76 + clhs117*clhs78;
    lhs(3,7)=clhs101*clhs109 + clhs102*clhs109 + clhs105*clhs115 + clhs108*clhs115 + clhs118*clhs76 + clhs118*clhs78;
    lhs(3,8)=0;
    lhs(3,9)=clhs119;
    lhs(3,10)=0;
    lhs(3,11)=clhs120;
    lhs(3,12)=0;
    lhs(3,13)=clhs119;
    lhs(3,14)=0;
    lhs(3,15)=clhs120;
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=-N1[0]*clhs121;
    lhs(4,5)=-N1[0]*clhs122;
    lhs(4,6)=-N1[0]*clhs123;
    lhs(4,7)=-N1[0]*clhs124;
    lhs(4,8)=clhs126;
    lhs(4,9)=0;
    lhs(4,10)=clhs127;
    lhs(4,11)=0;
    lhs(4,12)=clhs126;
    lhs(4,13)=0;
    lhs(4,14)=clhs127;
    lhs(4,15)=0;
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=-N1[0]*clhs128;
    lhs(5,5)=-N1[0]*clhs129;
    lhs(5,6)=-N1[0]*clhs130;
    lhs(5,7)=-N1[0]*clhs131;
    lhs(5,8)=0;
    lhs(5,9)=clhs126;
    lhs(5,10)=0;
    lhs(5,11)=clhs127;
    lhs(5,12)=0;
    lhs(5,13)=clhs126;
    lhs(5,14)=0;
    lhs(5,15)=clhs127;
    lhs(6,0)=0;
    lhs(6,1)=0;
    lhs(6,2)=0;
    lhs(6,3)=0;
    lhs(6,4)=-N1[1]*clhs121;
    lhs(6,5)=-N1[1]*clhs122;
    lhs(6,6)=-N1[1]*clhs123;
    lhs(6,7)=-N1[1]*clhs124;
    lhs(6,8)=clhs133;
    lhs(6,9)=0;
    lhs(6,10)=clhs134;
    lhs(6,11)=0;
    lhs(6,12)=clhs133;
    lhs(6,13)=0;
    lhs(6,14)=clhs134;
    lhs(6,15)=0;
    lhs(7,0)=0;
    lhs(7,1)=0;
    lhs(7,2)=0;
    lhs(7,3)=0;
    lhs(7,4)=-N1[1]*clhs128;
    lhs(7,5)=-N1[1]*clhs129;
    lhs(7,6)=-N1[1]*clhs130;
    lhs(7,7)=-N1[1]*clhs131;
    lhs(7,8)=0;
    lhs(7,9)=clhs133;
    lhs(7,10)=0;
    lhs(7,11)=clhs134;
    lhs(7,12)=0;
    lhs(7,13)=clhs133;
    lhs(7,14)=0;
    lhs(7,15)=clhs134;
    lhs(8,0)=clhs140;
    lhs(8,1)=clhs142;
    lhs(8,2)=clhs144;
    lhs(8,3)=clhs146;
    lhs(8,4)=-GPnormal[0]*clhs157 - GPnormal[1]*clhs166;
    lhs(8,5)=-GPnormal[0]*clhs170 - GPnormal[1]*clhs173;
    lhs(8,6)=-GPnormal[0]*clhs177 - GPnormal[1]*clhs180;
    lhs(8,7)=-GPnormal[0]*clhs184 - GPnormal[1]*clhs187;
    lhs(8,8)=clhs191;
    lhs(8,9)=clhs194;
    lhs(8,10)=clhs197;
    lhs(8,11)=clhs199;
    lhs(8,12)=clhs190;
    lhs(8,13)=clhs193;
    lhs(8,14)=clhs196;
    lhs(8,15)=clhs198;
    lhs(9,0)=clhs201;
    lhs(9,1)=clhs202;
    lhs(9,2)=clhs203;
    lhs(9,3)=clhs204;
    lhs(9,4)=-GPtangent1[0]*clhs157 - GPtangent1[1]*clhs166;
    lhs(9,5)=-GPtangent1[0]*clhs170 - GPtangent1[1]*clhs173;
    lhs(9,6)=-GPtangent1[0]*clhs177 - GPtangent1[1]*clhs180;
    lhs(9,7)=-GPtangent1[0]*clhs184 - GPtangent1[1]*clhs187;
    lhs(9,8)=clhs207;
    lhs(9,9)=clhs210;
    lhs(9,10)=clhs212;
    lhs(9,11)=clhs214;
    lhs(9,12)=clhs206;
    lhs(9,13)=clhs209;
    lhs(9,14)=clhs211;
    lhs(9,15)=clhs213;
    lhs(10,0)=clhs218;
    lhs(10,1)=clhs219;
    lhs(10,2)=clhs220;
    lhs(10,3)=clhs221;
    lhs(10,4)=-GPnormal[0]*clhs228 - GPnormal[1]*clhs234;
    lhs(10,5)=-GPnormal[0]*clhs236 - GPnormal[1]*clhs238;
    lhs(10,6)=-GPnormal[0]*clhs240 - GPnormal[1]*clhs242;
    lhs(10,7)=-GPnormal[0]*clhs244 - GPnormal[1]*clhs246;
    lhs(10,8)=clhs197;
    lhs(10,9)=clhs199;
    lhs(10,10)=clhs249;
    lhs(10,11)=clhs251;
    lhs(10,12)=clhs196;
    lhs(10,13)=clhs198;
    lhs(10,14)=clhs248;
    lhs(10,15)=clhs250;
    lhs(11,0)=clhs253;
    lhs(11,1)=clhs254;
    lhs(11,2)=clhs255;
    lhs(11,3)=clhs256;
    lhs(11,4)=-GPtangent1[0]*clhs228 - GPtangent1[1]*clhs234;
    lhs(11,5)=-GPtangent1[0]*clhs236 - GPtangent1[1]*clhs238;
    lhs(11,6)=-GPtangent1[0]*clhs240 - GPtangent1[1]*clhs242;
    lhs(11,7)=-GPtangent1[0]*clhs244 - GPtangent1[1]*clhs246;
    lhs(11,8)=clhs212;
    lhs(11,9)=clhs214;
    lhs(11,10)=clhs258;
    lhs(11,11)=clhs260;
    lhs(11,12)=clhs211;
    lhs(11,13)=clhs213;
    lhs(11,14)=clhs257;
    lhs(11,15)=clhs259;
    lhs(12,0)=clhs140;
    lhs(12,1)=clhs142;
    lhs(12,2)=clhs144;
    lhs(12,3)=clhs146;
    lhs(12,4)=-GPnormal[0]*clhs265 - GPnormal[1]*clhs270;
    lhs(12,5)=-GPnormal[0]*clhs272 - GPnormal[1]*clhs274;
    lhs(12,6)=-GPnormal[0]*clhs276 - GPnormal[1]*clhs278;
    lhs(12,7)=-GPnormal[0]*clhs280 - GPnormal[1]*clhs282;
    lhs(12,8)=clhs190;
    lhs(12,9)=clhs193;
    lhs(12,10)=clhs196;
    lhs(12,11)=clhs198;
    lhs(12,12)=clhs191;
    lhs(12,13)=clhs194;
    lhs(12,14)=clhs197;
    lhs(12,15)=clhs199;
    lhs(13,0)=clhs201;
    lhs(13,1)=clhs202;
    lhs(13,2)=clhs203;
    lhs(13,3)=clhs204;
    lhs(13,4)=-GPtangent1[0]*clhs265 - GPtangent1[1]*clhs270;
    lhs(13,5)=-GPtangent1[0]*clhs272 - GPtangent1[1]*clhs274;
    lhs(13,6)=-GPtangent1[0]*clhs276 - GPtangent1[1]*clhs278;
    lhs(13,7)=-GPtangent1[0]*clhs280 - GPtangent1[1]*clhs282;
    lhs(13,8)=clhs206;
    lhs(13,9)=clhs209;
    lhs(13,10)=clhs211;
    lhs(13,11)=clhs213;
    lhs(13,12)=clhs207;
    lhs(13,13)=clhs210;
    lhs(13,14)=clhs212;
    lhs(13,15)=clhs214;
    lhs(14,0)=clhs218;
    lhs(14,1)=clhs219;
    lhs(14,2)=clhs220;
    lhs(14,3)=clhs221;
    lhs(14,4)=-GPnormal[0]*clhs284 - GPnormal[1]*clhs286;
    lhs(14,5)=-GPnormal[0]*clhs287 - GPnormal[1]*clhs288;
    lhs(14,6)=-GPnormal[0]*clhs289 - GPnormal[1]*clhs290;
    lhs(14,7)=-GPnormal[0]*clhs291 - GPnormal[1]*clhs292;
    lhs(14,8)=clhs196;
    lhs(14,9)=clhs198;
    lhs(14,10)=clhs248;
    lhs(14,11)=clhs250;
    lhs(14,12)=clhs197;
    lhs(14,13)=clhs199;
    lhs(14,14)=clhs249;
    lhs(14,15)=clhs251;
    lhs(15,0)=clhs253;
    lhs(15,1)=clhs254;
    lhs(15,2)=clhs255;
    lhs(15,3)=clhs256;
    lhs(15,4)=-GPtangent1[0]*clhs284 - GPtangent1[1]*clhs286;
    lhs(15,5)=-GPtangent1[0]*clhs287 - GPtangent1[1]*clhs288;
    lhs(15,6)=-GPtangent1[0]*clhs289 - GPtangent1[1]*clhs290;
    lhs(15,7)=-GPtangent1[0]*clhs291 - GPtangent1[1]*clhs292;
    lhs(15,8)=clhs211;
    lhs(15,9)=clhs213;
    lhs(15,10)=clhs257;
    lhs(15,11)=clhs259;
    lhs(15,12)=clhs212;
    lhs(15,13)=clhs214;
    lhs(15,14)=clhs258;
    lhs(15,15)=clhs260;

    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,16,16> ComputeGaussPointStickLHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
//         const Matrix DPhi, 
        const double detJ, 
        const double mu, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    bounded_matrix<double,16,16> lhs;
    
    const Matrix normalmaster     = rContactData.Normal_m;
    const Matrix normalslave      = rContactData.Normal_s;
    const Matrix tan1slave        = rContactData.Tangent_xi_s;
    const Matrix lm               = rContactData.LagrangeMultipliers;
    const Matrix dlm              = rContactData.DoubleLagrangeMultipliers;
    const double Dt               = rContactData.Dt;
//     const double epsilon          = rContactData.epsilon;
//     const double epsilon_tangent  = rContactData.epsilon_tangent;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const Matrix v1 = rContactData.v1;
    const Matrix v2 = rContactData.v2;
    
    const std::vector<double> DeltaJs         = rContactData.DeltaJ_s;
    const std::vector<Matrix> DeltaNormals    = rContactData.Delta_Normal_s;
    const std::vector<Matrix> Deltatangentxis = rContactData.Delta_Tangent_xi_s;
    const std::vector<double> DeltaGap        = rContactData.DeltaGap;
    const std::vector<Vector> DeltaPhi        = rContactData.DeltaPhi;
    const std::vector<Vector> DeltaN2         = rContactData.DeltaN2;
    
    const double Dtan1slave11u111 =     Deltatangentxis[3](1,1);
    const double Dtan1slave11u110 =     Deltatangentxis[2](1,1);
    const double Dtan1slave11u101 =     Deltatangentxis[1](1,1);
    const double Dtan1slave11u100 =     Deltatangentxis[0](1,1);
    const double Dtan1slave10u111 =     Deltatangentxis[3](1,0);
    const double Dtan1slave10u110 =     Deltatangentxis[2](1,0);
    const double Dtan1slave10u101 =     Deltatangentxis[1](1,0);
    const double Dtan1slave10u100 =     Deltatangentxis[0](1,0);
    const double Dtan1slave01u111 =     Deltatangentxis[3](0,1);
    const double Dtan1slave01u110 =     Deltatangentxis[2](0,1);
    const double Dtan1slave01u101 =     Deltatangentxis[1](0,1);
    const double Dtan1slave01u100 =     Deltatangentxis[0](0,1);
    const double Dtan1slave00u111 =     Deltatangentxis[3](0,0);
    const double Dtan1slave00u110 =     Deltatangentxis[2](0,0);
    const double Dtan1slave00u101 =     Deltatangentxis[1](0,0);
    const double Dtan1slave00u100 =     Deltatangentxis[0](0,0);
    const double DPhi1u111 =     DeltaPhi[3][1];
    const double DPhi1u110 =     DeltaPhi[2][1];
    const double DPhi1u101 =     DeltaPhi[1][1];
    const double DPhi1u100 =     DeltaPhi[0][1];
    const double DPhi0u111 =     DeltaPhi[3][0];
    const double DPhi0u110 =     DeltaPhi[2][0];
    const double DPhi0u101 =     DeltaPhi[1][0];
    const double DPhi0u100 =     DeltaPhi[0][0];
    const double DN21u211 =     DeltaN2[7][1];
    const double DN21u210 =     DeltaN2[6][1];
    const double DN21u201 =     DeltaN2[5][1];
    const double DN21u200 =     DeltaN2[4][1];
    const double DN21u111 =     DeltaN2[3][1];
    const double DN21u110 =     DeltaN2[2][1];
    const double DN21u101 =     DeltaN2[1][1];
    const double DN21u100 =     DeltaN2[0][1];
    const double DN20u211 =     DeltaN2[7][0];
    const double DN20u210 =     DeltaN2[6][0];
    const double DN20u201 =     DeltaN2[5][0];
    const double DN20u200 =     DeltaN2[4][0];
    const double DN20u111 =     DeltaN2[3][0];
    const double DN20u110 =     DeltaN2[2][0];
    const double DN20u101 =     DeltaN2[1][0];
    const double DN20u100 =     DeltaN2[0][0];
    const double DdetJu111 =     DeltaJs[3];
    const double DdetJu110 =     DeltaJs[2];
    const double DdetJu101 =     DeltaJs[1];
    const double DdetJu100 =     DeltaJs[0];
 
    const double clhs0 =     1.0/Dt;
    const double clhs1 =     Phi[0]; // PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs2 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs3 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs4 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs5 =     clhs0*clhs1*clhs2*(GPnormal[0]*clhs3 + GPnormal[1]*clhs4);
    const double clhs6 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs7 =     N1[0]*clhs4 + N1[1]*clhs6;
    const double clhs8 =     Dt*v2(0,1);
    const double clhs9 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs10 =     DN20u200; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs11 =     Dt*v2(1,1);
    const double clhs12 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs13 =     DN21u200; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs14 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs15 =     N1[0]*clhs3 + N1[1]*clhs14;
    const double clhs16 =     Dt*v2(0,0);
    const double clhs17 =     Dt*v2(1,0);
    const double clhs18 =     clhs15*(clhs10*clhs16 + clhs13*clhs17 + clhs9) + clhs7*(clhs10*clhs8 + clhs11*clhs13);
    const double clhs19 =     -clhs18*clhs5;
    const double clhs20 =     DN20u201; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs21 =     DN21u201; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs22 =     clhs15*(clhs16*clhs20 + clhs17*clhs21) + clhs7*(clhs11*clhs21 + clhs20*clhs8 + clhs9);
    const double clhs23 =     -clhs22*clhs5;
    const double clhs24 =     DN20u210; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs25 =     DN21u210; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs26 =     clhs15*(clhs12 + clhs16*clhs24 + clhs17*clhs25) + clhs7*(clhs11*clhs25 + clhs24*clhs8);
    const double clhs27 =     -clhs26*clhs5;
    const double clhs28 =     DN20u211; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs29 =     DN21u211; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs30 =     clhs15*(clhs16*clhs28 + clhs17*clhs29) + clhs7*(clhs11*clhs29 + clhs12 + clhs28*clhs8);
    const double clhs31 =     -clhs30*clhs5;
    const double clhs32 =     Dtan1slave00u100; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs33 =     N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - clhs12*clhs17 - clhs16*clhs9;
    const double clhs34 =     N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - clhs11*clhs12 - clhs8*clhs9;
    const double clhs35 =     clhs15*clhs33 + clhs34*clhs7;
    const double clhs36 =     clhs1*clhs2*clhs35;
    const double clhs37 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs38 =     clhs1*clhs3*clhs35;
    const double clhs39 =     DPhi0u100; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs40 =     clhs2*clhs3*clhs35;
    const double clhs41 =     DN20u100; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs42 =     DN21u100; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs43 =     -N1[0];
    const double clhs44 =     Dtan1slave10u100; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs45 =     Dtan1slave01u100; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs46 =     Dtan1slave11u100; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs47 =     clhs15*(clhs16*clhs41 + clhs17*clhs42 + clhs43) - clhs33*(N1[0]*clhs32 + N1[1]*clhs44) - clhs34*(N1[0]*clhs45 + N1[1]*clhs46) + clhs7*(clhs11*clhs42 + clhs41*clhs8);
    const double clhs48 =     clhs1*clhs2*clhs47;
    const double clhs49 =     -clhs3*clhs48 + clhs32*clhs36 + clhs37*clhs38 + clhs39*clhs40;
    const double clhs50 =     clhs1*clhs35*clhs4;
    const double clhs51 =     clhs2*clhs35*clhs4;
    const double clhs52 =     clhs36*clhs45 + clhs37*clhs50 + clhs39*clhs51 - clhs4*clhs48;
    const double clhs53 =     clhs0*(GPnormal[0]*clhs49 + GPnormal[1]*clhs52);
    const double clhs54 =     Dtan1slave00u101; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs55 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs56 =     DPhi0u101; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs57 =     DN20u101; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs58 =     DN21u101; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs59 =     Dtan1slave10u101; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs60 =     Dtan1slave01u101; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs61 =     Dtan1slave11u101; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs62 =     clhs15*(clhs16*clhs57 + clhs17*clhs58) - clhs33*(N1[0]*clhs54 + N1[1]*clhs59) - clhs34*(N1[0]*clhs60 + N1[1]*clhs61) + clhs7*(clhs11*clhs58 + clhs43 + clhs57*clhs8);
    const double clhs63 =     clhs1*clhs2*clhs62;
    const double clhs64 =     -clhs3*clhs63 + clhs36*clhs54 + clhs38*clhs55 + clhs40*clhs56;
    const double clhs65 =     clhs36*clhs60 - clhs4*clhs63 + clhs50*clhs55 + clhs51*clhs56;
    const double clhs66 =     clhs0*(GPnormal[0]*clhs64 + GPnormal[1]*clhs65);
    const double clhs67 =     Dtan1slave00u110; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs68 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs69 =     DPhi0u110; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs70 =     DN20u110; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs71 =     DN21u110; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs72 =     -N1[1];
    const double clhs73 =     Dtan1slave10u110; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs74 =     Dtan1slave01u110; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs75 =     Dtan1slave11u110; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs76 =     clhs15*(clhs16*clhs70 + clhs17*clhs71 + clhs72) - clhs33*(N1[0]*clhs67 + N1[1]*clhs73) - clhs34*(N1[0]*clhs74 + N1[1]*clhs75) + clhs7*(clhs11*clhs71 + clhs70*clhs8);
    const double clhs77 =     clhs1*clhs2*clhs76;
    const double clhs78 =     -clhs3*clhs77 + clhs36*clhs67 + clhs38*clhs68 + clhs40*clhs69;
    const double clhs79 =     clhs36*clhs74 - clhs4*clhs77 + clhs50*clhs68 + clhs51*clhs69;
    const double clhs80 =     clhs0*(GPnormal[0]*clhs78 + GPnormal[1]*clhs79);
    const double clhs81 =     Dtan1slave00u111; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs82 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs83 =     DPhi0u111; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs84 =     DN20u111; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs85 =     DN21u111; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs86 =     Dtan1slave10u111; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs87 =     Dtan1slave01u111; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs88 =     Dtan1slave11u111; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs89 =     clhs15*(clhs16*clhs84 + clhs17*clhs85) - clhs33*(N1[0]*clhs81 + N1[1]*clhs86) - clhs34*(N1[0]*clhs87 + N1[1]*clhs88) + clhs7*(clhs11*clhs85 + clhs72 + clhs8*clhs84);
    const double clhs90 =     clhs1*clhs2*clhs89;
    const double clhs91 =     -clhs3*clhs90 + clhs36*clhs81 + clhs38*clhs82 + clhs40*clhs83;
    const double clhs92 =     clhs36*clhs87 - clhs4*clhs90 + clhs50*clhs82 + clhs51*clhs83;
    const double clhs93 =     clhs0*(GPnormal[0]*clhs91 + GPnormal[1]*clhs92);
    const double clhs94 =     clhs0*clhs1*clhs2*(GPtangent1[0]*clhs3 + GPtangent1[1]*clhs4);
    const double clhs95 =     -clhs18*clhs94;
    const double clhs96 =     -clhs22*clhs94;
    const double clhs97 =     -clhs26*clhs94;
    const double clhs98 =     -clhs30*clhs94;
    const double clhs99 =     clhs0*(GPtangent1[0]*clhs49 + GPtangent1[1]*clhs52);
    const double clhs100 =     clhs0*(GPtangent1[0]*clhs64 + GPtangent1[1]*clhs65);
    const double clhs101 =     clhs0*(GPtangent1[0]*clhs78 + GPtangent1[1]*clhs79);
    const double clhs102 =     clhs0*(GPtangent1[0]*clhs91 + GPtangent1[1]*clhs92);
    const double clhs103 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs104 =     clhs0*clhs103*clhs2*(GPnormal[0]*clhs14 + GPnormal[1]*clhs6);
    const double clhs105 =     -clhs104*clhs18;
    const double clhs106 =     -clhs104*clhs22;
    const double clhs107 =     -clhs104*clhs26;
    const double clhs108 =     -clhs104*clhs30;
    const double clhs109 =     clhs103*clhs2*clhs35;
    const double clhs110 =     clhs103*clhs14*clhs35;
    const double clhs111 =     DPhi1u100; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs112 =     clhs14*clhs2*clhs35;
    const double clhs113 =     clhs103*clhs2*clhs47;
    const double clhs114 =     clhs109*clhs44 + clhs110*clhs37 + clhs111*clhs112 - clhs113*clhs14;
    const double clhs115 =     clhs103*clhs35*clhs6;
    const double clhs116 =     clhs2*clhs35*clhs6;
    const double clhs117 =     clhs109*clhs46 + clhs111*clhs116 - clhs113*clhs6 + clhs115*clhs37;
    const double clhs118 =     clhs0*(GPnormal[0]*clhs114 + GPnormal[1]*clhs117);
    const double clhs119 =     DPhi1u101; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs120 =     clhs103*clhs2*clhs62;
    const double clhs121 =     clhs109*clhs59 + clhs110*clhs55 + clhs112*clhs119 - clhs120*clhs14;
    const double clhs122 =     clhs109*clhs61 + clhs115*clhs55 + clhs116*clhs119 - clhs120*clhs6;
    const double clhs123 =     clhs0*(GPnormal[0]*clhs121 + GPnormal[1]*clhs122);
    const double clhs124 =     DPhi1u110; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs125 =     clhs103*clhs2*clhs76;
    const double clhs126 =     clhs109*clhs73 + clhs110*clhs68 + clhs112*clhs124 - clhs125*clhs14;
    const double clhs127 =     clhs109*clhs75 + clhs115*clhs68 + clhs116*clhs124 - clhs125*clhs6;
    const double clhs128 =     clhs0*(GPnormal[0]*clhs126 + GPnormal[1]*clhs127);
    const double clhs129 =     DPhi1u111; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs130 =     clhs103*clhs2*clhs89;
    const double clhs131 =     clhs109*clhs86 + clhs110*clhs82 + clhs112*clhs129 - clhs130*clhs14;
    const double clhs132 =     clhs109*clhs88 + clhs115*clhs82 + clhs116*clhs129 - clhs130*clhs6;
    const double clhs133 =     clhs0*(GPnormal[0]*clhs131 + GPnormal[1]*clhs132);
    const double clhs134 =     clhs0*clhs103*clhs2*(GPtangent1[0]*clhs14 + GPtangent1[1]*clhs6);
    const double clhs135 =     -clhs134*clhs18;
    const double clhs136 =     -clhs134*clhs22;
    const double clhs137 =     -clhs134*clhs26;
    const double clhs138 =     -clhs134*clhs30;
    const double clhs139 =     clhs0*(GPtangent1[0]*clhs114 + GPtangent1[1]*clhs117);
    const double clhs140 =     clhs0*(GPtangent1[0]*clhs121 + GPtangent1[1]*clhs122);
    const double clhs141 =     clhs0*(GPtangent1[0]*clhs126 + GPtangent1[1]*clhs127);
    const double clhs142 =     clhs0*(GPtangent1[0]*clhs131 + GPtangent1[1]*clhs132);

    lhs(0,0)=0;
    lhs(0,1)=0;
    lhs(0,2)=0;
    lhs(0,3)=0;
    lhs(0,4)=0;
    lhs(0,5)=0;
    lhs(0,6)=0;
    lhs(0,7)=0;
    lhs(0,8)=0;
    lhs(0,9)=0;
    lhs(0,10)=0;
    lhs(0,11)=0;
    lhs(0,12)=0;
    lhs(0,13)=0;
    lhs(0,14)=0;
    lhs(0,15)=0;
    lhs(1,0)=0;
    lhs(1,1)=0;
    lhs(1,2)=0;
    lhs(1,3)=0;
    lhs(1,4)=0;
    lhs(1,5)=0;
    lhs(1,6)=0;
    lhs(1,7)=0;
    lhs(1,8)=0;
    lhs(1,9)=0;
    lhs(1,10)=0;
    lhs(1,11)=0;
    lhs(1,12)=0;
    lhs(1,13)=0;
    lhs(1,14)=0;
    lhs(1,15)=0;
    lhs(2,0)=0;
    lhs(2,1)=0;
    lhs(2,2)=0;
    lhs(2,3)=0;
    lhs(2,4)=0;
    lhs(2,5)=0;
    lhs(2,6)=0;
    lhs(2,7)=0;
    lhs(2,8)=0;
    lhs(2,9)=0;
    lhs(2,10)=0;
    lhs(2,11)=0;
    lhs(2,12)=0;
    lhs(2,13)=0;
    lhs(2,14)=0;
    lhs(2,15)=0;
    lhs(3,0)=0;
    lhs(3,1)=0;
    lhs(3,2)=0;
    lhs(3,3)=0;
    lhs(3,4)=0;
    lhs(3,5)=0;
    lhs(3,6)=0;
    lhs(3,7)=0;
    lhs(3,8)=0;
    lhs(3,9)=0;
    lhs(3,10)=0;
    lhs(3,11)=0;
    lhs(3,12)=0;
    lhs(3,13)=0;
    lhs(3,14)=0;
    lhs(3,15)=0;
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=0;
    lhs(4,5)=0;
    lhs(4,6)=0;
    lhs(4,7)=0;
    lhs(4,8)=0;
    lhs(4,9)=0;
    lhs(4,10)=0;
    lhs(4,11)=0;
    lhs(4,12)=0;
    lhs(4,13)=0;
    lhs(4,14)=0;
    lhs(4,15)=0;
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=0;
    lhs(5,5)=0;
    lhs(5,6)=0;
    lhs(5,7)=0;
    lhs(5,8)=0;
    lhs(5,9)=0;
    lhs(5,10)=0;
    lhs(5,11)=0;
    lhs(5,12)=0;
    lhs(5,13)=0;
    lhs(5,14)=0;
    lhs(5,15)=0;
    lhs(6,0)=0;
    lhs(6,1)=0;
    lhs(6,2)=0;
    lhs(6,3)=0;
    lhs(6,4)=0;
    lhs(6,5)=0;
    lhs(6,6)=0;
    lhs(6,7)=0;
    lhs(6,8)=0;
    lhs(6,9)=0;
    lhs(6,10)=0;
    lhs(6,11)=0;
    lhs(6,12)=0;
    lhs(6,13)=0;
    lhs(6,14)=0;
    lhs(6,15)=0;
    lhs(7,0)=0;
    lhs(7,1)=0;
    lhs(7,2)=0;
    lhs(7,3)=0;
    lhs(7,4)=0;
    lhs(7,5)=0;
    lhs(7,6)=0;
    lhs(7,7)=0;
    lhs(7,8)=0;
    lhs(7,9)=0;
    lhs(7,10)=0;
    lhs(7,11)=0;
    lhs(7,12)=0;
    lhs(7,13)=0;
    lhs(7,14)=0;
    lhs(7,15)=0;
    lhs(8,0)=clhs19;
    lhs(8,1)=clhs23;
    lhs(8,2)=clhs27;
    lhs(8,3)=clhs31;
    lhs(8,4)=clhs53;
    lhs(8,5)=clhs66;
    lhs(8,6)=clhs80;
    lhs(8,7)=clhs93;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(8,12)=0;
    lhs(8,13)=0;
    lhs(8,14)=0;
    lhs(8,15)=0;
    lhs(9,0)=clhs95;
    lhs(9,1)=clhs96;
    lhs(9,2)=clhs97;
    lhs(9,3)=clhs98;
    lhs(9,4)=clhs99;
    lhs(9,5)=clhs100;
    lhs(9,6)=clhs101;
    lhs(9,7)=clhs102;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(9,12)=0;
    lhs(9,13)=0;
    lhs(9,14)=0;
    lhs(9,15)=0;
    lhs(10,0)=clhs105;
    lhs(10,1)=clhs106;
    lhs(10,2)=clhs107;
    lhs(10,3)=clhs108;
    lhs(10,4)=clhs118;
    lhs(10,5)=clhs123;
    lhs(10,6)=clhs128;
    lhs(10,7)=clhs133;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(10,12)=0;
    lhs(10,13)=0;
    lhs(10,14)=0;
    lhs(10,15)=0;
    lhs(11,0)=clhs135;
    lhs(11,1)=clhs136;
    lhs(11,2)=clhs137;
    lhs(11,3)=clhs138;
    lhs(11,4)=clhs139;
    lhs(11,5)=clhs140;
    lhs(11,6)=clhs141;
    lhs(11,7)=clhs142;
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;
    lhs(11,12)=0;
    lhs(11,13)=0;
    lhs(11,14)=0;
    lhs(11,15)=0;
    lhs(12,0)=clhs19;
    lhs(12,1)=clhs23;
    lhs(12,2)=clhs27;
    lhs(12,3)=clhs31;
    lhs(12,4)=clhs53;
    lhs(12,5)=clhs66;
    lhs(12,6)=clhs80;
    lhs(12,7)=clhs93;
    lhs(12,8)=0;
    lhs(12,9)=0;
    lhs(12,10)=0;
    lhs(12,11)=0;
    lhs(12,12)=0;
    lhs(12,13)=0;
    lhs(12,14)=0;
    lhs(12,15)=0;
    lhs(13,0)=clhs95;
    lhs(13,1)=clhs96;
    lhs(13,2)=clhs97;
    lhs(13,3)=clhs98;
    lhs(13,4)=clhs99;
    lhs(13,5)=clhs100;
    lhs(13,6)=clhs101;
    lhs(13,7)=clhs102;
    lhs(13,8)=0;
    lhs(13,9)=0;
    lhs(13,10)=0;
    lhs(13,11)=0;
    lhs(13,12)=0;
    lhs(13,13)=0;
    lhs(13,14)=0;
    lhs(13,15)=0;
    lhs(14,0)=clhs105;
    lhs(14,1)=clhs106;
    lhs(14,2)=clhs107;
    lhs(14,3)=clhs108;
    lhs(14,4)=clhs118;
    lhs(14,5)=clhs123;
    lhs(14,6)=clhs128;
    lhs(14,7)=clhs133;
    lhs(14,8)=0;
    lhs(14,9)=0;
    lhs(14,10)=0;
    lhs(14,11)=0;
    lhs(14,12)=0;
    lhs(14,13)=0;
    lhs(14,14)=0;
    lhs(14,15)=0;
    lhs(15,0)=clhs135;
    lhs(15,1)=clhs136;
    lhs(15,2)=clhs137;
    lhs(15,3)=clhs138;
    lhs(15,4)=clhs139;
    lhs(15,5)=clhs140;
    lhs(15,6)=clhs141;
    lhs(15,7)=clhs142;
    lhs(15,8)=0;
    lhs(15,9)=0;
    lhs(15,10)=0;
    lhs(15,11)=0;
    lhs(15,12)=0;
    lhs(15,13)=0;
    lhs(15,14)=0;
    lhs(15,15)=0;

    
    return lhs;
}
    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,16,16> ComputeGaussPointSlipLHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
//         const Matrix DPhi, 
        const double detJ, 
        const double mu, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    bounded_matrix<double,16,16> lhs;
    
    const Matrix normalmaster     = rContactData.Normal_m;
    const Matrix normalslave      = rContactData.Normal_s;
    const Matrix tan1slave        = rContactData.Tangent_xi_s;
    const Matrix lm               = rContactData.LagrangeMultipliers;
    const Matrix dlm              = rContactData.DoubleLagrangeMultipliers;
    const double Dt               = rContactData.Dt;
//     const double epsilon          = rContactData.epsilon;
//     const double epsilon_tangent  = rContactData.epsilon_tangent;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
//     const double sign_tangpress = boost::math::sign(augmented_tangent_lm);
    
    const Matrix v1 = rContactData.v1;
    const Matrix v2 = rContactData.v2;

    const std::vector<double> DeltaJs         = rContactData.DeltaJ_s;
    const std::vector<Matrix> DeltaNormals    = rContactData.Delta_Normal_s;
    const std::vector<Matrix> Deltatangentxis = rContactData.Delta_Tangent_xi_s;
    const std::vector<double> DeltaGap        = rContactData.DeltaGap;
    const std::vector<Vector> DeltaPhi        = rContactData.DeltaPhi;
    const std::vector<Vector> DeltaN2         = rContactData.DeltaN2;
    
    const double Dtan1slave11u111 =     Deltatangentxis[3](1,1);
    const double Dtan1slave11u110 =     Deltatangentxis[2](1,1);
    const double Dtan1slave11u101 =     Deltatangentxis[1](1,1);
    const double Dtan1slave11u100 =     Deltatangentxis[0](1,1);
    const double Dtan1slave10u111 =     Deltatangentxis[3](1,0);
    const double Dtan1slave10u110 =     Deltatangentxis[2](1,0);
    const double Dtan1slave10u101 =     Deltatangentxis[1](1,0);
    const double Dtan1slave10u100 =     Deltatangentxis[0](1,0);
    const double Dtan1slave01u111 =     Deltatangentxis[3](0,1);
    const double Dtan1slave01u110 =     Deltatangentxis[2](0,1);
    const double Dtan1slave01u101 =     Deltatangentxis[1](0,1);
    const double Dtan1slave01u100 =     Deltatangentxis[0](0,1);
    const double Dtan1slave00u111 =     Deltatangentxis[3](0,0);
    const double Dtan1slave00u110 =     Deltatangentxis[2](0,0);
    const double Dtan1slave00u101 =     Deltatangentxis[1](0,0);
    const double Dtan1slave00u100 =     Deltatangentxis[0](0,0);
    const double DPhi1u111 =     DeltaPhi[3][1];
    const double DPhi1u110 =     DeltaPhi[2][1];
    const double DPhi1u101 =     DeltaPhi[1][1];
    const double DPhi1u100 =     DeltaPhi[0][1];
    const double DPhi0u111 =     DeltaPhi[3][0];
    const double DPhi0u110 =     DeltaPhi[2][0];
    const double DPhi0u101 =     DeltaPhi[1][0];
    const double DPhi0u100 =     DeltaPhi[0][0];
    const double DN21u211 =     DeltaN2[7][1];
    const double DN21u210 =     DeltaN2[6][1];
    const double DN21u201 =     DeltaN2[5][1];
    const double DN21u200 =     DeltaN2[4][1];
    const double DN21u111 =     DeltaN2[3][1];
    const double DN21u110 =     DeltaN2[2][1];
    const double DN21u101 =     DeltaN2[1][1];
    const double DN21u100 =     DeltaN2[0][1];
    const double DN20u211 =     DeltaN2[7][0];
    const double DN20u210 =     DeltaN2[6][0];
    const double DN20u201 =     DeltaN2[5][0];
    const double DN20u200 =     DeltaN2[4][0];
    const double DN20u111 =     DeltaN2[3][0];
    const double DN20u110 =     DeltaN2[2][0];
    const double DN20u101 =     DeltaN2[1][0];
    const double DN20u100 =     DeltaN2[0][0];
    const double DdetJu111 =     DeltaJs[3];
    const double DdetJu110 =     DeltaJs[2];
    const double DdetJu101 =     DeltaJs[1];
    const double DdetJu100 =     DeltaJs[0];
 
    const double clhs0 =     1.0/Dt;
    const double clhs1 =     Phi[0]; // PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs2 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs3 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs4 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs5 =     clhs0*clhs1*clhs2*(GPnormal[0]*clhs3 + GPnormal[1]*clhs4);
    const double clhs6 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs7 =     N1[0]*clhs4 + N1[1]*clhs6;
    const double clhs8 =     Dt*v2(0,1);
    const double clhs9 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs10 =     DN20u200; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs11 =     Dt*v2(1,1);
    const double clhs12 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs13 =     DN21u200; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs14 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs15 =     N1[0]*clhs3 + N1[1]*clhs14;
    const double clhs16 =     Dt*v2(0,0);
    const double clhs17 =     Dt*v2(1,0);
    const double clhs18 =     clhs15*(clhs10*clhs16 + clhs13*clhs17 + clhs9) + clhs7*(clhs10*clhs8 + clhs11*clhs13);
    const double clhs19 =     -clhs18*clhs5;
    const double clhs20 =     DN20u201; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs21 =     DN21u201; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs22 =     clhs15*(clhs16*clhs20 + clhs17*clhs21) + clhs7*(clhs11*clhs21 + clhs20*clhs8 + clhs9);
    const double clhs23 =     -clhs22*clhs5;
    const double clhs24 =     DN20u210; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs25 =     DN21u210; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs26 =     clhs15*(clhs12 + clhs16*clhs24 + clhs17*clhs25) + clhs7*(clhs11*clhs25 + clhs24*clhs8);
    const double clhs27 =     -clhs26*clhs5;
    const double clhs28 =     DN20u211; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs29 =     DN21u211; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs30 =     clhs15*(clhs16*clhs28 + clhs17*clhs29) + clhs7*(clhs11*clhs29 + clhs12 + clhs28*clhs8);
    const double clhs31 =     -clhs30*clhs5;
    const double clhs32 =     Dtan1slave00u100; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs33 =     N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - clhs12*clhs17 - clhs16*clhs9;
    const double clhs34 =     N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - clhs11*clhs12 - clhs8*clhs9;
    const double clhs35 =     clhs15*clhs33 + clhs34*clhs7;
    const double clhs36 =     clhs1*clhs2*clhs35;
    const double clhs37 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs38 =     clhs1*clhs3*clhs35;
    const double clhs39 =     DPhi0u100; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs40 =     clhs2*clhs3*clhs35;
    const double clhs41 =     DN20u100; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs42 =     DN21u100; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs43 =     -N1[0];
    const double clhs44 =     Dtan1slave10u100; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs45 =     Dtan1slave01u100; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs46 =     Dtan1slave11u100; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs47 =     clhs15*(clhs16*clhs41 + clhs17*clhs42 + clhs43) - clhs33*(N1[0]*clhs32 + N1[1]*clhs44) - clhs34*(N1[0]*clhs45 + N1[1]*clhs46) + clhs7*(clhs11*clhs42 + clhs41*clhs8);
    const double clhs48 =     clhs1*clhs2*clhs47;
    const double clhs49 =     -clhs3*clhs48 + clhs32*clhs36 + clhs37*clhs38 + clhs39*clhs40;
    const double clhs50 =     clhs1*clhs35*clhs4;
    const double clhs51 =     clhs2*clhs35*clhs4;
    const double clhs52 =     clhs36*clhs45 + clhs37*clhs50 + clhs39*clhs51 - clhs4*clhs48;
    const double clhs53 =     clhs0*(GPnormal[0]*clhs49 + GPnormal[1]*clhs52);
    const double clhs54 =     Dtan1slave00u101; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs55 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs56 =     DPhi0u101; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs57 =     DN20u101; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs58 =     DN21u101; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs59 =     Dtan1slave10u101; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs60 =     Dtan1slave01u101; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs61 =     Dtan1slave11u101; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs62 =     clhs15*(clhs16*clhs57 + clhs17*clhs58) - clhs33*(N1[0]*clhs54 + N1[1]*clhs59) - clhs34*(N1[0]*clhs60 + N1[1]*clhs61) + clhs7*(clhs11*clhs58 + clhs43 + clhs57*clhs8);
    const double clhs63 =     clhs1*clhs2*clhs62;
    const double clhs64 =     -clhs3*clhs63 + clhs36*clhs54 + clhs38*clhs55 + clhs40*clhs56;
    const double clhs65 =     clhs36*clhs60 - clhs4*clhs63 + clhs50*clhs55 + clhs51*clhs56;
    const double clhs66 =     clhs0*(GPnormal[0]*clhs64 + GPnormal[1]*clhs65);
    const double clhs67 =     Dtan1slave00u110; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs68 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs69 =     DPhi0u110; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs70 =     DN20u110; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs71 =     DN21u110; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs72 =     -N1[1];
    const double clhs73 =     Dtan1slave10u110; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs74 =     Dtan1slave01u110; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs75 =     Dtan1slave11u110; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs76 =     clhs15*(clhs16*clhs70 + clhs17*clhs71 + clhs72) - clhs33*(N1[0]*clhs67 + N1[1]*clhs73) - clhs34*(N1[0]*clhs74 + N1[1]*clhs75) + clhs7*(clhs11*clhs71 + clhs70*clhs8);
    const double clhs77 =     clhs1*clhs2*clhs76;
    const double clhs78 =     -clhs3*clhs77 + clhs36*clhs67 + clhs38*clhs68 + clhs40*clhs69;
    const double clhs79 =     clhs36*clhs74 - clhs4*clhs77 + clhs50*clhs68 + clhs51*clhs69;
    const double clhs80 =     clhs0*(GPnormal[0]*clhs78 + GPnormal[1]*clhs79);
    const double clhs81 =     Dtan1slave00u111; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs82 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs83 =     DPhi0u111; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs84 =     DN20u111; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs85 =     DN21u111; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs86 =     Dtan1slave10u111; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs87 =     Dtan1slave01u111; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs88 =     Dtan1slave11u111; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs89 =     clhs15*(clhs16*clhs84 + clhs17*clhs85) - clhs33*(N1[0]*clhs81 + N1[1]*clhs86) - clhs34*(N1[0]*clhs87 + N1[1]*clhs88) + clhs7*(clhs11*clhs85 + clhs72 + clhs8*clhs84);
    const double clhs90 =     clhs1*clhs2*clhs89;
    const double clhs91 =     -clhs3*clhs90 + clhs36*clhs81 + clhs38*clhs82 + clhs40*clhs83;
    const double clhs92 =     clhs36*clhs87 - clhs4*clhs90 + clhs50*clhs82 + clhs51*clhs83;
    const double clhs93 =     clhs0*(GPnormal[0]*clhs91 + GPnormal[1]*clhs92);
    const double clhs94 =     clhs0*clhs1*clhs2*(GPtangent1[0]*clhs3 + GPtangent1[1]*clhs4);
    const double clhs95 =     -clhs18*clhs94;
    const double clhs96 =     -clhs22*clhs94;
    const double clhs97 =     -clhs26*clhs94;
    const double clhs98 =     -clhs30*clhs94;
    const double clhs99 =     clhs0*(GPtangent1[0]*clhs49 + GPtangent1[1]*clhs52);
    const double clhs100 =     clhs0*(GPtangent1[0]*clhs64 + GPtangent1[1]*clhs65);
    const double clhs101 =     clhs0*(GPtangent1[0]*clhs78 + GPtangent1[1]*clhs79);
    const double clhs102 =     clhs0*(GPtangent1[0]*clhs91 + GPtangent1[1]*clhs92);
    const double clhs103 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs104 =     clhs0*clhs103*clhs2*(GPnormal[0]*clhs14 + GPnormal[1]*clhs6);
    const double clhs105 =     -clhs104*clhs18;
    const double clhs106 =     -clhs104*clhs22;
    const double clhs107 =     -clhs104*clhs26;
    const double clhs108 =     -clhs104*clhs30;
    const double clhs109 =     clhs103*clhs2*clhs35;
    const double clhs110 =     clhs103*clhs14*clhs35;
    const double clhs111 =     DPhi1u100; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs112 =     clhs14*clhs2*clhs35;
    const double clhs113 =     clhs103*clhs2*clhs47;
    const double clhs114 =     clhs109*clhs44 + clhs110*clhs37 + clhs111*clhs112 - clhs113*clhs14;
    const double clhs115 =     clhs103*clhs35*clhs6;
    const double clhs116 =     clhs2*clhs35*clhs6;
    const double clhs117 =     clhs109*clhs46 + clhs111*clhs116 - clhs113*clhs6 + clhs115*clhs37;
    const double clhs118 =     clhs0*(GPnormal[0]*clhs114 + GPnormal[1]*clhs117);
    const double clhs119 =     DPhi1u101; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs120 =     clhs103*clhs2*clhs62;
    const double clhs121 =     clhs109*clhs59 + clhs110*clhs55 + clhs112*clhs119 - clhs120*clhs14;
    const double clhs122 =     clhs109*clhs61 + clhs115*clhs55 + clhs116*clhs119 - clhs120*clhs6;
    const double clhs123 =     clhs0*(GPnormal[0]*clhs121 + GPnormal[1]*clhs122);
    const double clhs124 =     DPhi1u110; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs125 =     clhs103*clhs2*clhs76;
    const double clhs126 =     clhs109*clhs73 + clhs110*clhs68 + clhs112*clhs124 - clhs125*clhs14;
    const double clhs127 =     clhs109*clhs75 + clhs115*clhs68 + clhs116*clhs124 - clhs125*clhs6;
    const double clhs128 =     clhs0*(GPnormal[0]*clhs126 + GPnormal[1]*clhs127);
    const double clhs129 =     DPhi1u111; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs130 =     clhs103*clhs2*clhs89;
    const double clhs131 =     clhs109*clhs86 + clhs110*clhs82 + clhs112*clhs129 - clhs130*clhs14;
    const double clhs132 =     clhs109*clhs88 + clhs115*clhs82 + clhs116*clhs129 - clhs130*clhs6;
    const double clhs133 =     clhs0*(GPnormal[0]*clhs131 + GPnormal[1]*clhs132);
    const double clhs134 =     clhs0*clhs103*clhs2*(GPtangent1[0]*clhs14 + GPtangent1[1]*clhs6);
    const double clhs135 =     -clhs134*clhs18;
    const double clhs136 =     -clhs134*clhs22;
    const double clhs137 =     -clhs134*clhs26;
    const double clhs138 =     -clhs134*clhs30;
    const double clhs139 =     clhs0*(GPtangent1[0]*clhs114 + GPtangent1[1]*clhs117);
    const double clhs140 =     clhs0*(GPtangent1[0]*clhs121 + GPtangent1[1]*clhs122);
    const double clhs141 =     clhs0*(GPtangent1[0]*clhs126 + GPtangent1[1]*clhs127);
    const double clhs142 =     clhs0*(GPtangent1[0]*clhs131 + GPtangent1[1]*clhs132);

    lhs(0,0)=0;
    lhs(0,1)=0;
    lhs(0,2)=0;
    lhs(0,3)=0;
    lhs(0,4)=0;
    lhs(0,5)=0;
    lhs(0,6)=0;
    lhs(0,7)=0;
    lhs(0,8)=0;
    lhs(0,9)=0;
    lhs(0,10)=0;
    lhs(0,11)=0;
    lhs(0,12)=0;
    lhs(0,13)=0;
    lhs(0,14)=0;
    lhs(0,15)=0;
    lhs(1,0)=0;
    lhs(1,1)=0;
    lhs(1,2)=0;
    lhs(1,3)=0;
    lhs(1,4)=0;
    lhs(1,5)=0;
    lhs(1,6)=0;
    lhs(1,7)=0;
    lhs(1,8)=0;
    lhs(1,9)=0;
    lhs(1,10)=0;
    lhs(1,11)=0;
    lhs(1,12)=0;
    lhs(1,13)=0;
    lhs(1,14)=0;
    lhs(1,15)=0;
    lhs(2,0)=0;
    lhs(2,1)=0;
    lhs(2,2)=0;
    lhs(2,3)=0;
    lhs(2,4)=0;
    lhs(2,5)=0;
    lhs(2,6)=0;
    lhs(2,7)=0;
    lhs(2,8)=0;
    lhs(2,9)=0;
    lhs(2,10)=0;
    lhs(2,11)=0;
    lhs(2,12)=0;
    lhs(2,13)=0;
    lhs(2,14)=0;
    lhs(2,15)=0;
    lhs(3,0)=0;
    lhs(3,1)=0;
    lhs(3,2)=0;
    lhs(3,3)=0;
    lhs(3,4)=0;
    lhs(3,5)=0;
    lhs(3,6)=0;
    lhs(3,7)=0;
    lhs(3,8)=0;
    lhs(3,9)=0;
    lhs(3,10)=0;
    lhs(3,11)=0;
    lhs(3,12)=0;
    lhs(3,13)=0;
    lhs(3,14)=0;
    lhs(3,15)=0;
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=0;
    lhs(4,5)=0;
    lhs(4,6)=0;
    lhs(4,7)=0;
    lhs(4,8)=0;
    lhs(4,9)=0;
    lhs(4,10)=0;
    lhs(4,11)=0;
    lhs(4,12)=0;
    lhs(4,13)=0;
    lhs(4,14)=0;
    lhs(4,15)=0;
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=0;
    lhs(5,5)=0;
    lhs(5,6)=0;
    lhs(5,7)=0;
    lhs(5,8)=0;
    lhs(5,9)=0;
    lhs(5,10)=0;
    lhs(5,11)=0;
    lhs(5,12)=0;
    lhs(5,13)=0;
    lhs(5,14)=0;
    lhs(5,15)=0;
    lhs(6,0)=0;
    lhs(6,1)=0;
    lhs(6,2)=0;
    lhs(6,3)=0;
    lhs(6,4)=0;
    lhs(6,5)=0;
    lhs(6,6)=0;
    lhs(6,7)=0;
    lhs(6,8)=0;
    lhs(6,9)=0;
    lhs(6,10)=0;
    lhs(6,11)=0;
    lhs(6,12)=0;
    lhs(6,13)=0;
    lhs(6,14)=0;
    lhs(6,15)=0;
    lhs(7,0)=0;
    lhs(7,1)=0;
    lhs(7,2)=0;
    lhs(7,3)=0;
    lhs(7,4)=0;
    lhs(7,5)=0;
    lhs(7,6)=0;
    lhs(7,7)=0;
    lhs(7,8)=0;
    lhs(7,9)=0;
    lhs(7,10)=0;
    lhs(7,11)=0;
    lhs(7,12)=0;
    lhs(7,13)=0;
    lhs(7,14)=0;
    lhs(7,15)=0;
    lhs(8,0)=clhs19;
    lhs(8,1)=clhs23;
    lhs(8,2)=clhs27;
    lhs(8,3)=clhs31;
    lhs(8,4)=clhs53;
    lhs(8,5)=clhs66;
    lhs(8,6)=clhs80;
    lhs(8,7)=clhs93;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(8,12)=0;
    lhs(8,13)=0;
    lhs(8,14)=0;
    lhs(8,15)=0;
    lhs(9,0)=clhs95;
    lhs(9,1)=clhs96;
    lhs(9,2)=clhs97;
    lhs(9,3)=clhs98;
    lhs(9,4)=clhs99;
    lhs(9,5)=clhs100;
    lhs(9,6)=clhs101;
    lhs(9,7)=clhs102;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(9,12)=0;
    lhs(9,13)=0;
    lhs(9,14)=0;
    lhs(9,15)=0;
    lhs(10,0)=clhs105;
    lhs(10,1)=clhs106;
    lhs(10,2)=clhs107;
    lhs(10,3)=clhs108;
    lhs(10,4)=clhs118;
    lhs(10,5)=clhs123;
    lhs(10,6)=clhs128;
    lhs(10,7)=clhs133;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(10,12)=0;
    lhs(10,13)=0;
    lhs(10,14)=0;
    lhs(10,15)=0;
    lhs(11,0)=clhs135;
    lhs(11,1)=clhs136;
    lhs(11,2)=clhs137;
    lhs(11,3)=clhs138;
    lhs(11,4)=clhs139;
    lhs(11,5)=clhs140;
    lhs(11,6)=clhs141;
    lhs(11,7)=clhs142;
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;
    lhs(11,12)=0;
    lhs(11,13)=0;
    lhs(11,14)=0;
    lhs(11,15)=0;
    lhs(12,0)=clhs19;
    lhs(12,1)=clhs23;
    lhs(12,2)=clhs27;
    lhs(12,3)=clhs31;
    lhs(12,4)=clhs53;
    lhs(12,5)=clhs66;
    lhs(12,6)=clhs80;
    lhs(12,7)=clhs93;
    lhs(12,8)=0;
    lhs(12,9)=0;
    lhs(12,10)=0;
    lhs(12,11)=0;
    lhs(12,12)=0;
    lhs(12,13)=0;
    lhs(12,14)=0;
    lhs(12,15)=0;
    lhs(13,0)=clhs95;
    lhs(13,1)=clhs96;
    lhs(13,2)=clhs97;
    lhs(13,3)=clhs98;
    lhs(13,4)=clhs99;
    lhs(13,5)=clhs100;
    lhs(13,6)=clhs101;
    lhs(13,7)=clhs102;
    lhs(13,8)=0;
    lhs(13,9)=0;
    lhs(13,10)=0;
    lhs(13,11)=0;
    lhs(13,12)=0;
    lhs(13,13)=0;
    lhs(13,14)=0;
    lhs(13,15)=0;
    lhs(14,0)=clhs105;
    lhs(14,1)=clhs106;
    lhs(14,2)=clhs107;
    lhs(14,3)=clhs108;
    lhs(14,4)=clhs118;
    lhs(14,5)=clhs123;
    lhs(14,6)=clhs128;
    lhs(14,7)=clhs133;
    lhs(14,8)=0;
    lhs(14,9)=0;
    lhs(14,10)=0;
    lhs(14,11)=0;
    lhs(14,12)=0;
    lhs(14,13)=0;
    lhs(14,14)=0;
    lhs(14,15)=0;
    lhs(15,0)=clhs135;
    lhs(15,1)=clhs136;
    lhs(15,2)=clhs137;
    lhs(15,3)=clhs138;
    lhs(15,4)=clhs139;
    lhs(15,5)=clhs140;
    lhs(15,6)=clhs141;
    lhs(15,7)=clhs142;
    lhs(15,8)=0;
    lhs(15,9)=0;
    lhs(15,10)=0;
    lhs(15,11)=0;
    lhs(15,12)=0;
    lhs(15,13)=0;
    lhs(15,14)=0;
    lhs(15,15)=0;

    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,16,16> ComputeGaussPointInactiveLHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
//         const Matrix DPhi, 
        const double detJ, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    bounded_matrix<double,16,16> lhs;
    
//     const Matrix normalmaster     = rContactData.Normal_m;
//     const Matrix normalslave      = rContactData.Normal_s;
//     const Matrix tan1slave        = rContactData.Tangent_xi_s;
//     const Matrix lm               = rContactData.LagrangeMultipliers;
//     const Matrix dlm              = rContactData.DoubleLagrangeMultipliers;
//     const double Dt               = rContactData.Dt;
//     const double epsilon          = rContactData.epsilon;
//     const double epsilon_tangent  = rContactData.epsilon_tangent;
//     
//     const Vector GPnormal     = prod(trans(normalslave), N1);
//     const Vector GPtangent1   = prod(trans(tan1slave), N1);
//     
//     const Matrix v1 = rContactData.v1;
//     const Matrix v2 = rContactData.v2;
// 
//     const std::vector<double> DeltaJs         = rContactData.DeltaJ_s;
//     const std::vector<Matrix> DeltaNormals    = rContactData.Delta_Normal_s;
//     const std::vector<Matrix> Deltatangentxis = rContactData.Delta_Tangent_xi_s;
//     const std::vector<double> DeltaGap        = rContactData.DeltaGap;
//     const std::vector<Vector> DeltaPhi        = rContactData.DeltaPhi;
//     const std::vector<Vector> DeltaN2         = rContactData.DeltaN2;
//     
//substitute_derivatives_variables_inactive 
//substitute_inactive_lhs
    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,16> ComputeGaussPointActiveRHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    array_1d<double,16> rhs;
    
    const Matrix normalslave     = rContactData.Normal_s;
    const Matrix tan1slave       = rContactData.Tangent_xi_s;
    const Matrix lm              = rContactData.LagrangeMultipliers;
    const Matrix dlm             = rContactData.DoubleLagrangeMultipliers;
//     const double Dt             = rContactData.Dt;
    const double epsilon         = rContactData.epsilon;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const Matrix v1 = rContactData.v1;
    const Matrix v2 = rContactData.v2;
    
    const double crhs0 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs1 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs2 =     Phi[0]; // PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs3 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs4 =     crhs2*dlm(0,0) + crhs3*dlm(1,0);
    const double crhs5 =     crhs2*lm(0,0);
    const double crhs6 =     crhs3*lm(1,0);
    const double crhs7 =     crhs4 + crhs5 + crhs6;
    const double crhs8 =     crhs1*crhs7;
    const double crhs9 =     crhs2*dlm(0,1) + crhs3*dlm(1,1);
    const double crhs10 =     crhs2*lm(0,1);
    const double crhs11 =     crhs3*lm(1,1);
    const double crhs12 =     crhs10 + crhs11 + crhs9;
    const double crhs13 =     crhs1*crhs12;
    const double crhs14 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs15 =     N1[0]*crhs1;
    const double crhs16 =     N1[1]*crhs1;
    const double crhs17 =     crhs1*crhs2;
    const double crhs18 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs19 =     crhs18*normalslave(0,0); // CRHS18*NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs20 =     epsilon*(crhs4 - crhs5 - crhs6);
    const double crhs21 =     -crhs19 + crhs20;
    const double crhs22 =     crhs18*normalslave(0,1); // CRHS18*NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs23 =     epsilon*(-crhs10 - crhs11 + crhs9);
    const double crhs24 =     -crhs22 + crhs23;
    const double crhs25 =     crhs1*crhs3;
    const double crhs26 =     crhs18*normalslave(1,0); // CRHS18*NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs27 =     crhs20 - crhs26;
    const double crhs28 =     crhs18*normalslave(1,1); // CRHS18*NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs29 =     crhs23 - crhs28;
    const double crhs30 =     crhs19 + crhs20;
    const double crhs31 =     crhs22 + crhs23;
    const double crhs32 =     crhs20 + crhs26;
    const double crhs33 =     crhs23 + crhs28;

    rhs[0]=-crhs0*crhs8;
    rhs[1]=-crhs0*crhs13;
    rhs[2]=-crhs14*crhs8;
    rhs[3]=-crhs13*crhs14;
    rhs[4]=crhs15*crhs7;
    rhs[5]=crhs12*crhs15;
    rhs[6]=crhs16*crhs7;
    rhs[7]=crhs12*crhs16;
    rhs[8]=-crhs17*(GPnormal[0]*crhs21 + GPnormal[1]*crhs24);
    rhs[9]=-crhs17*(GPtangent1[0]*crhs21 + GPtangent1[1]*crhs24);
    rhs[10]=-crhs25*(GPnormal[0]*crhs27 + GPnormal[1]*crhs29);
    rhs[11]=-crhs25*(GPtangent1[0]*crhs27 + GPtangent1[1]*crhs29);
    rhs[12]=crhs17*(GPnormal[0]*crhs30 + GPnormal[1]*crhs31);
    rhs[13]=crhs17*(GPtangent1[0]*crhs30 + GPtangent1[1]*crhs31);
    rhs[14]=crhs25*(GPnormal[0]*crhs32 + GPnormal[1]*crhs33);
    rhs[15]=crhs25*(GPtangent1[0]*crhs32 + GPtangent1[1]*crhs33);

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,16> ComputeGaussPointStickRHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const double mu, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    array_1d<double,16> rhs;
    
    const Matrix normalslave      = rContactData.Normal_s;
    const Matrix tan1slave        = rContactData.Tangent_xi_s;
    const Matrix lm               = rContactData.LagrangeMultipliers;
    const Matrix dlm              = rContactData.DoubleLagrangeMultipliers;
    const double Dt               = rContactData.Dt;
//     const double epsilon          = rContactData.epsilon;
//     const double epsilon_tangent  = rContactData.epsilon_tangent;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const Matrix v1 = rContactData.v1;
    const Matrix v2 = rContactData.v2;
    
    const double crhs0 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs1 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs2 =     1.0/Dt;
    const double crhs3 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs4 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs5 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs6 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs7 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs8 =     (N1[0]*crhs0 + N1[1]*crhs4)*(N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - crhs5*(Dt*v2(0,0)) - crhs6*(Dt*v2(1,0))) + (N1[0]*crhs1 + N1[1]*crhs7)*(N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - crhs5*(Dt*v2(0,1)) - crhs6*(Dt*v2(1,1)));
    const double crhs9 =     crhs2*crhs3*crhs8*Phi[0]; // CRHS2*CRHS3*CRHS8*PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs10 =     -crhs9*(GPnormal[0]*crhs0 + GPnormal[1]*crhs1);
    const double crhs11 =     -crhs9*(GPtangent1[0]*crhs0 + GPtangent1[1]*crhs1);
    const double crhs12 =     crhs2*crhs3*crhs8*Phi[1]; // CRHS2*CRHS3*CRHS8*PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs13 =     -crhs12*(GPnormal[0]*crhs4 + GPnormal[1]*crhs7);
    const double crhs14 =     -crhs12*(GPtangent1[0]*crhs4 + GPtangent1[1]*crhs7);

    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=crhs10;
    rhs[9]=crhs11;
    rhs[10]=crhs13;
    rhs[11]=crhs14;
    rhs[12]=crhs10;
    rhs[13]=crhs11;
    rhs[14]=crhs13;
    rhs[15]=crhs14;

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,16> ComputeGaussPointSlipRHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const double mu, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    array_1d<double,16> rhs;
    
    const Matrix normalslave      = rContactData.Normal_s;
    const Matrix tan1slave        = rContactData.Tangent_xi_s;
    const Matrix lm               = rContactData.LagrangeMultipliers;
    const Matrix dlm              = rContactData.DoubleLagrangeMultipliers;
    const double Dt               = rContactData.Dt;
//     const double epsilon          = rContactData.epsilon;
//     const double epsilon_tangent  = rContactData.epsilon_tangent;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
//     const double sign_tangpress = boost::math::sign(augmented_tangent_lm);
    
    const Matrix v1 = rContactData.v1;
    const Matrix v2 = rContactData.v2;
    
    const double crhs0 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs1 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs2 =     1.0/Dt;
    const double crhs3 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs4 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs5 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs6 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs7 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs8 =     (N1[0]*crhs0 + N1[1]*crhs4)*(N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - crhs5*(Dt*v2(0,0)) - crhs6*(Dt*v2(1,0))) + (N1[0]*crhs1 + N1[1]*crhs7)*(N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - crhs5*(Dt*v2(0,1)) - crhs6*(Dt*v2(1,1)));
    const double crhs9 =     crhs2*crhs3*crhs8*Phi[0]; // CRHS2*CRHS3*CRHS8*PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs10 =     -crhs9*(GPnormal[0]*crhs0 + GPnormal[1]*crhs1);
    const double crhs11 =     -crhs9*(GPtangent1[0]*crhs0 + GPtangent1[1]*crhs1);
    const double crhs12 =     crhs2*crhs3*crhs8*Phi[1]; // CRHS2*CRHS3*CRHS8*PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs13 =     -crhs12*(GPnormal[0]*crhs4 + GPnormal[1]*crhs7);
    const double crhs14 =     -crhs12*(GPtangent1[0]*crhs4 + GPtangent1[1]*crhs7);

    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=crhs10;
    rhs[9]=crhs11;
    rhs[10]=crhs13;
    rhs[11]=crhs14;
    rhs[12]=crhs10;
    rhs[13]=crhs11;
    rhs[14]=crhs13;
    rhs[15]=crhs14;

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,16> ComputeGaussPointInactiveRHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    array_1d<double,16> rhs;
    
//     const Matrix normalslave      = rContactData.Normal_s;
//     const Matrix tan1slave        = rContactData.Tangent_xi_s;
//     const Matrix lm               = rContactData.LagrangeMultipliers;
//     const Matrix dlm              = rContactData.DoubleLagrangeMultipliers;
//     const double Dt               = rContactData.Dt;
//     const double epsilon          = rContactData.epsilon;
//     const double epsilon_tangent  = rContactData.epsilon_tangent;
//     
//     const Vector GPnormal     = prod(trans(normalslave), N1);
//     const Vector GPtangent1   = prod(trans(tan1slave), N1);
//     
//     const Matrix v1 = rContactData.v1;
//     const Matrix v2 = rContactData.v2;
//     
//substitute_inactive_rhs
    
    return rhs;
}
    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,2,2> ComputeDeltaDe(
        const Vector N1, 
        const ContactData& rContactData,
        const unsigned int derivative_index
        )
{
    bounded_matrix<double,2,2> DeltaDe;
    
    const double DeltaDetJ = rContactData.DeltaJ_s[derivative_index];
    
    DeltaDe(0,0) = DeltaDetJ * N1[0];
    DeltaDe(0,1) = 0;
    DeltaDe(1,0) = 0;
    DeltaDe(1,1) = DeltaDetJ * N1[1];
    
    return DeltaDe;
}
    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,2,2> ComputeDeltaMe(
        const Vector N1, 
        const ContactData& rContactData,
        const unsigned int derivative_index
        )
{
    bounded_matrix<double,2,2> DeltaMe;

    const double DeltaDetJ = rContactData.DeltaJ_s[derivative_index];
    
    DeltaMe = DeltaDetJ * outer_prod(N1, N1);
    
    return DeltaMe;
}

private:
};// class Contact2D2N2NDLM
}
#endif /* KRATOS_CONTACT2D2N2NDLM defined */
