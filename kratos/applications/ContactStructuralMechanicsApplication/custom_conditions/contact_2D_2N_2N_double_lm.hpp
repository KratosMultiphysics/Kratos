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
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const ContactData<2,2>& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    bounded_matrix<double,16,16> lhs;
    
    const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
    const bounded_matrix<double, 2, 2> dlm             = rContactData.DoubleLagrangeMultipliers;
//     const double Dt              = rContactData.Dt;
    const double epsilon = rContactData.epsilon;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
    
    const array_1d<double,2> GPnormal     = prod(trans(normalslave), N1);
    const array_1d<double,2> GPtangent1   = prod(trans(tan1slave), N1);
    
    const bounded_matrix<double, 2, 2> v1 = rContactData.v1;
    const bounded_matrix<double, 2, 2> v2 = rContactData.v2;

    const std::vector<double> DeltaJs  = rContactData.DeltaJ_s;
    const std::vector<double> DeltaGap = rContactData.DeltaGap;
    const std::vector<array_1d<double,2>> DeltaPhi = rContactData.DeltaPhi;
    const std::vector<array_1d<double,2>> DeltaN1  = rContactData.DeltaN1;
    const std::vector<array_1d<double,2>> DeltaN2  = rContactData.DeltaN2;
    const std::vector<bounded_matrix<double, 2, 2>> DeltaNormals = rContactData.Delta_Normal_s;
    
    const double clhs0 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs1 =     DeltaN2[4][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
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
    const double clhs12 =     DeltaN2[5][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs13 =     DeltaN2[6][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs14 =     DeltaN2[7][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs15 =     DeltaJs[0]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs16 =     clhs0*clhs7;
    const double clhs17 =     DeltaN2[0][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs18 =     clhs2*clhs7;
    const double clhs19 =     clhs0*clhs10;
    const double clhs20 =     clhs10*clhs2;
    const double clhs21 =     clhs0*clhs2;
    const double clhs22 =     DeltaPhi[0][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs23 =     clhs22*dlm(0,0);
    const double clhs24 =     DeltaPhi[0][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs25 =     clhs24*dlm(1,0);
    const double clhs26 =     clhs23 + clhs25;
    const double clhs27 =     clhs22*lm(0,0);
    const double clhs28 =     clhs24*lm(1,0);
    const double clhs29 =     clhs27 + clhs28;
    const double clhs30 =     DeltaJs[1]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs31 =     DeltaN2[1][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs32 =     DeltaPhi[1][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs33 =     clhs32*dlm(0,0);
    const double clhs34 =     DeltaPhi[1][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs35 =     clhs34*dlm(1,0);
    const double clhs36 =     clhs33 + clhs35;
    const double clhs37 =     clhs32*lm(0,0);
    const double clhs38 =     clhs34*lm(1,0);
    const double clhs39 =     clhs37 + clhs38;
    const double clhs40 =     DeltaJs[2]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs41 =     DeltaN2[2][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs42 =     DeltaPhi[2][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs43 =     clhs42*dlm(0,0);
    const double clhs44 =     DeltaPhi[2][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs45 =     clhs44*dlm(1,0);
    const double clhs46 =     clhs43 + clhs45;
    const double clhs47 =     clhs42*lm(0,0);
    const double clhs48 =     clhs44*lm(1,0);
    const double clhs49 =     clhs47 + clhs48;
    const double clhs50 =     DeltaJs[3]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs51 =     DeltaN2[3][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs52 =     DeltaPhi[3][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs53 =     clhs52*dlm(0,0);
    const double clhs54 =     DeltaPhi[3][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs55 =     clhs54*dlm(1,0);
    const double clhs56 =     clhs53 + clhs55;
    const double clhs57 =     clhs52*lm(0,0);
    const double clhs58 =     clhs54*lm(1,0);
    const double clhs59 =     clhs57 + clhs58;
    const double clhs60 =     clhs21*clhs3;
    const double clhs61 =     clhs21*clhs5;
    const double clhs62 =     clhs3*dlm(0,1);
    const double clhs63 =     clhs5*dlm(1,1);
    const double clhs64 =     clhs62 + clhs63;
    const double clhs65 =     clhs3*lm(0,1);
    const double clhs66 =     clhs5*lm(1,1);
    const double clhs67 =     clhs65 + clhs66;
    const double clhs68 =     clhs2*(clhs64 + clhs67);
    const double clhs69 =     clhs0*clhs64;
    const double clhs70 =     clhs2*clhs64;
    const double clhs71 =     clhs0*clhs67;
    const double clhs72 =     clhs2*clhs67;
    const double clhs73 =     clhs22*dlm(0,1);
    const double clhs74 =     clhs24*dlm(1,1);
    const double clhs75 =     clhs73 + clhs74;
    const double clhs76 =     clhs22*lm(0,1);
    const double clhs77 =     clhs24*lm(1,1);
    const double clhs78 =     clhs76 + clhs77;
    const double clhs79 =     clhs32*dlm(0,1);
    const double clhs80 =     clhs34*dlm(1,1);
    const double clhs81 =     clhs79 + clhs80;
    const double clhs82 =     clhs32*lm(0,1);
    const double clhs83 =     clhs34*lm(1,1);
    const double clhs84 =     clhs82 + clhs83;
    const double clhs85 =     clhs42*dlm(0,1);
    const double clhs86 =     clhs44*dlm(1,1);
    const double clhs87 =     clhs85 + clhs86;
    const double clhs88 =     clhs42*lm(0,1);
    const double clhs89 =     clhs44*lm(1,1);
    const double clhs90 =     clhs88 + clhs89;
    const double clhs91 =     clhs52*dlm(0,1);
    const double clhs92 =     clhs54*dlm(1,1);
    const double clhs93 =     clhs91 + clhs92;
    const double clhs94 =     clhs52*lm(0,1);
    const double clhs95 =     clhs54*lm(1,1);
    const double clhs96 =     clhs94 + clhs95;
    const double clhs97 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs98 =     DeltaN2[4][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs99 =     DeltaN2[5][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs100 =     DeltaN2[6][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs101 =     DeltaN2[7][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs102 =     clhs7*clhs97;
    const double clhs103 =     DeltaN2[0][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs104 =     clhs10*clhs97;
    const double clhs105 =     clhs2*clhs97;
    const double clhs106 =     DeltaN2[1][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs107 =     DeltaN2[2][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs108 =     DeltaN2[3][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs109 =     clhs105*clhs3;
    const double clhs110 =     clhs105*clhs5;
    const double clhs111 =     clhs64*clhs97;
    const double clhs112 =     clhs67*clhs97;
    const double clhs113 =     N1[0]; // N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs114 =     DeltaN1[4][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs115 =     DeltaN1[5][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs116 =     DeltaN1[6][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs117 =     DeltaN1[7][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs118 =     clhs113*clhs7;
    const double clhs119 =     DeltaN1[0][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs120 =     clhs10*clhs113;
    const double clhs121 =     clhs113*clhs2;
    const double clhs122 =     DeltaN1[1][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs123 =     DeltaN1[2][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs124 =     DeltaN1[3][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs125 =     -clhs121*clhs3;
    const double clhs126 =     -clhs121*clhs5;
    const double clhs127 =     clhs113*clhs64;
    const double clhs128 =     clhs113*clhs67;
    const double clhs129 =     N1[1]; // N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs130 =     DeltaN1[4][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs131 =     DeltaN1[5][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs132 =     DeltaN1[6][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs133 =     DeltaN1[7][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs134 =     clhs129*clhs7;
    const double clhs135 =     DeltaN1[0][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs136 =     clhs10*clhs129;
    const double clhs137 =     clhs129*clhs2;
    const double clhs138 =     DeltaN1[1][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs139 =     DeltaN1[2][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs140 =     DeltaN1[3][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs141 =     -clhs137*clhs3;
    const double clhs142 =     -clhs137*clhs5;
    const double clhs143 =     clhs129*clhs64;
    const double clhs144 =     clhs129*clhs67;
    const double clhs145 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs146 =     DeltaGap[4]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs147 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs148 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs149 =     clhs2*clhs3*(GPnormal[0]*clhs147 + GPnormal[1]*clhs148);
    const double clhs150 =     -clhs146*clhs149;
    const double clhs151 =     DeltaGap[5]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs152 =     -clhs149*clhs151;
    const double clhs153 =     DeltaGap[6]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs154 =     -clhs149*clhs153;
    const double clhs155 =     DeltaGap[7]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs156 =     -clhs149*clhs155;
    const double clhs157 =     clhs145*clhs2*clhs3;
    const double clhs158 =     DeltaGap[0]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs159 =     clhs147*clhs2*clhs3;
    const double clhs160 =     clhs145*clhs147*clhs3;
    const double clhs161 =     clhs145*clhs147*clhs2;
    const double clhs162 =     clhs15*clhs160 + clhs157*DeltaNormals[0](0,0) + clhs158*clhs159 + clhs161*clhs22; // CLHS15*CLHS160 + CLHS157*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS158*CLHS159 + CLHS161*CLHS22
    const double clhs163 =     clhs10 - clhs4 - clhs6;
    const double clhs164 =     clhs163*clhs3*epsilon;
    const double clhs165 =     clhs163*clhs2*epsilon;
    const double clhs166 =     clhs2*epsilon*(-clhs23 - clhs25 + clhs29);
    const double clhs167 =     clhs15*clhs164 + clhs162 + clhs165*clhs22 + clhs166*clhs3;
    const double clhs168 =     clhs148*clhs2*clhs3;
    const double clhs169 =     clhs145*clhs148*clhs3;
    const double clhs170 =     clhs145*clhs148*clhs2;
    const double clhs171 =     clhs15*clhs169 + clhs157*DeltaNormals[0](0,1) + clhs158*clhs168 + clhs170*clhs22; // CLHS15*CLHS169 + CLHS157*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS158*CLHS168 + CLHS170*CLHS22
    const double clhs172 =     -clhs62 - clhs63 + clhs67;
    const double clhs173 =     clhs172*clhs3*epsilon;
    const double clhs174 =     clhs172*clhs2*epsilon;
    const double clhs175 =     clhs2*epsilon*(-clhs73 - clhs74 + clhs78);
    const double clhs176 =     clhs15*clhs173 + clhs171 + clhs174*clhs22 + clhs175*clhs3;
    const double clhs177 =     DeltaGap[1]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs178 =     clhs157*DeltaNormals[1](0,0) + clhs159*clhs177 + clhs160*clhs30 + clhs161*clhs32; // CLHS157*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS159*CLHS177 + CLHS160*CLHS30 + CLHS161*CLHS32
    const double clhs179 =     clhs2*epsilon*(-clhs33 - clhs35 + clhs39);
    const double clhs180 =     clhs164*clhs30 + clhs165*clhs32 + clhs178 + clhs179*clhs3;
    const double clhs181 =     clhs157*DeltaNormals[1](0,1) + clhs168*clhs177 + clhs169*clhs30 + clhs170*clhs32; // CLHS157*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS168*CLHS177 + CLHS169*CLHS30 + CLHS170*CLHS32
    const double clhs182 =     clhs2*epsilon*(-clhs79 - clhs80 + clhs84);
    const double clhs183 =     clhs173*clhs30 + clhs174*clhs32 + clhs181 + clhs182*clhs3;
    const double clhs184 =     DeltaGap[2]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs185 =     clhs157*DeltaNormals[2](0,0) + clhs159*clhs184 + clhs160*clhs40 + clhs161*clhs42; // CLHS157*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS159*CLHS184 + CLHS160*CLHS40 + CLHS161*CLHS42
    const double clhs186 =     clhs2*epsilon*(-clhs43 - clhs45 + clhs49);
    const double clhs187 =     clhs164*clhs40 + clhs165*clhs42 + clhs185 + clhs186*clhs3;
    const double clhs188 =     clhs157*DeltaNormals[2](0,1) + clhs168*clhs184 + clhs169*clhs40 + clhs170*clhs42; // CLHS157*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS168*CLHS184 + CLHS169*CLHS40 + CLHS170*CLHS42
    const double clhs189 =     clhs2*epsilon*(-clhs85 - clhs86 + clhs90);
    const double clhs190 =     clhs173*clhs40 + clhs174*clhs42 + clhs188 + clhs189*clhs3;
    const double clhs191 =     DeltaGap[3]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs192 =     clhs157*DeltaNormals[3](0,0) + clhs159*clhs191 + clhs160*clhs50 + clhs161*clhs52; // CLHS157*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS159*CLHS191 + CLHS160*CLHS50 + CLHS161*CLHS52
    const double clhs193 =     clhs2*epsilon*(-clhs53 - clhs55 + clhs59);
    const double clhs194 =     clhs164*clhs50 + clhs165*clhs52 + clhs192 + clhs193*clhs3;
    const double clhs195 =     clhs157*DeltaNormals[3](0,1) + clhs168*clhs191 + clhs169*clhs50 + clhs170*clhs52; // CLHS157*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS168*CLHS191 + CLHS169*CLHS50 + CLHS170*CLHS52
    const double clhs196 =     clhs2*epsilon*(-clhs91 - clhs92 + clhs96);
    const double clhs197 =     clhs173*clhs50 + clhs174*clhs52 + clhs195 + clhs196*clhs3;
    const double clhs198 =     std::pow(clhs3, 2);
    const double clhs199 =     GPnormal[0]*clhs2*epsilon;
    const double clhs200 =     clhs198*clhs199;
    const double clhs201 =     -clhs200;
    const double clhs202 =     GPnormal[1]*clhs2*epsilon;
    const double clhs203 =     clhs198*clhs202;
    const double clhs204 =     -clhs203;
    const double clhs205 =     clhs3*clhs5;
    const double clhs206 =     clhs199*clhs205;
    const double clhs207 =     -clhs206;
    const double clhs208 =     clhs202*clhs205;
    const double clhs209 =     -clhs208;
    const double clhs210 =     clhs2*clhs3*(GPtangent1[0]*clhs147 + GPtangent1[1]*clhs148);
    const double clhs211 =     -clhs146*clhs210;
    const double clhs212 =     -clhs151*clhs210;
    const double clhs213 =     -clhs153*clhs210;
    const double clhs214 =     -clhs155*clhs210;
    const double clhs215 =     GPtangent1[0]*clhs2*epsilon;
    const double clhs216 =     clhs198*clhs215;
    const double clhs217 =     -clhs216;
    const double clhs218 =     GPtangent1[1]*clhs2*epsilon;
    const double clhs219 =     clhs198*clhs218;
    const double clhs220 =     -clhs219;
    const double clhs221 =     clhs205*clhs215;
    const double clhs222 =     -clhs221;
    const double clhs223 =     clhs205*clhs218;
    const double clhs224 =     -clhs223;
    const double clhs225 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs226 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs227 =     clhs2*clhs5*(GPnormal[0]*clhs225 + GPnormal[1]*clhs226);
    const double clhs228 =     -clhs146*clhs227;
    const double clhs229 =     -clhs151*clhs227;
    const double clhs230 =     -clhs153*clhs227;
    const double clhs231 =     -clhs155*clhs227;
    const double clhs232 =     clhs145*clhs2*clhs5;
    const double clhs233 =     clhs2*clhs225*clhs5;
    const double clhs234 =     clhs145*clhs225*clhs5;
    const double clhs235 =     clhs145*clhs2*clhs225;
    const double clhs236 =     clhs15*clhs234 + clhs158*clhs233 + clhs232*DeltaNormals[0](1,0) + clhs235*clhs24; // CLHS15*CLHS234 + CLHS158*CLHS233 + CLHS232*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS235*CLHS24
    const double clhs237 =     clhs163*clhs5*epsilon;
    const double clhs238 =     clhs15*clhs237 + clhs165*clhs24 + clhs166*clhs5 + clhs236;
    const double clhs239 =     clhs2*clhs226*clhs5;
    const double clhs240 =     clhs145*clhs226*clhs5;
    const double clhs241 =     clhs145*clhs2*clhs226;
    const double clhs242 =     clhs15*clhs240 + clhs158*clhs239 + clhs232*DeltaNormals[0](1,1) + clhs24*clhs241; // CLHS15*CLHS240 + CLHS158*CLHS239 + CLHS232*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS24*CLHS241
    const double clhs243 =     clhs172*clhs5*epsilon;
    const double clhs244 =     clhs15*clhs243 + clhs174*clhs24 + clhs175*clhs5 + clhs242;
    const double clhs245 =     clhs177*clhs233 + clhs232*DeltaNormals[1](1,0) + clhs234*clhs30 + clhs235*clhs34; // CLHS177*CLHS233 + CLHS232*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS234*CLHS30 + CLHS235*CLHS34
    const double clhs246 =     clhs165*clhs34 + clhs179*clhs5 + clhs237*clhs30 + clhs245;
    const double clhs247 =     clhs177*clhs239 + clhs232*DeltaNormals[1](1,1) + clhs240*clhs30 + clhs241*clhs34; // CLHS177*CLHS239 + CLHS232*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS240*CLHS30 + CLHS241*CLHS34
    const double clhs248 =     clhs174*clhs34 + clhs182*clhs5 + clhs243*clhs30 + clhs247;
    const double clhs249 =     clhs184*clhs233 + clhs232*DeltaNormals[2](1,0) + clhs234*clhs40 + clhs235*clhs44; // CLHS184*CLHS233 + CLHS232*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS234*CLHS40 + CLHS235*CLHS44
    const double clhs250 =     clhs165*clhs44 + clhs186*clhs5 + clhs237*clhs40 + clhs249;
    const double clhs251 =     clhs184*clhs239 + clhs232*DeltaNormals[2](1,1) + clhs240*clhs40 + clhs241*clhs44; // CLHS184*CLHS239 + CLHS232*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS240*CLHS40 + CLHS241*CLHS44
    const double clhs252 =     clhs174*clhs44 + clhs189*clhs5 + clhs243*clhs40 + clhs251;
    const double clhs253 =     clhs191*clhs233 + clhs232*DeltaNormals[3](1,0) + clhs234*clhs50 + clhs235*clhs54; // CLHS191*CLHS233 + CLHS232*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS234*CLHS50 + CLHS235*CLHS54
    const double clhs254 =     clhs165*clhs54 + clhs193*clhs5 + clhs237*clhs50 + clhs253;
    const double clhs255 =     clhs191*clhs239 + clhs232*DeltaNormals[3](1,1) + clhs240*clhs50 + clhs241*clhs54; // CLHS191*CLHS239 + CLHS232*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS240*CLHS50 + CLHS241*CLHS54
    const double clhs256 =     clhs174*clhs54 + clhs196*clhs5 + clhs243*clhs50 + clhs255;
    const double clhs257 =     std::pow(clhs5, 2);
    const double clhs258 =     clhs199*clhs257;
    const double clhs259 =     -clhs258;
    const double clhs260 =     clhs202*clhs257;
    const double clhs261 =     -clhs260;
    const double clhs262 =     clhs2*clhs5*(GPtangent1[0]*clhs225 + GPtangent1[1]*clhs226);
    const double clhs263 =     -clhs146*clhs262;
    const double clhs264 =     -clhs151*clhs262;
    const double clhs265 =     -clhs153*clhs262;
    const double clhs266 =     -clhs155*clhs262;
    const double clhs267 =     clhs215*clhs257;
    const double clhs268 =     -clhs267;
    const double clhs269 =     clhs218*clhs257;
    const double clhs270 =     -clhs269;
    const double clhs271 =     clhs7 - clhs8 - clhs9;
    const double clhs272 =     clhs271*clhs3*epsilon;
    const double clhs273 =     clhs2*clhs271*epsilon;
    const double clhs274 =     clhs2*epsilon*(clhs26 - clhs27 - clhs28);
    const double clhs275 =     clhs15*clhs272 + clhs162 + clhs22*clhs273 + clhs274*clhs3;
    const double clhs276 =     clhs64 - clhs65 - clhs66;
    const double clhs277 =     clhs276*clhs3*epsilon;
    const double clhs278 =     clhs2*clhs276*epsilon;
    const double clhs279 =     clhs2*epsilon*(clhs75 - clhs76 - clhs77);
    const double clhs280 =     clhs15*clhs277 + clhs171 + clhs22*clhs278 + clhs279*clhs3;
    const double clhs281 =     clhs2*epsilon*(clhs36 - clhs37 - clhs38);
    const double clhs282 =     clhs178 + clhs272*clhs30 + clhs273*clhs32 + clhs281*clhs3;
    const double clhs283 =     clhs2*epsilon*(clhs81 - clhs82 - clhs83);
    const double clhs284 =     clhs181 + clhs277*clhs30 + clhs278*clhs32 + clhs283*clhs3;
    const double clhs285 =     clhs2*epsilon*(clhs46 - clhs47 - clhs48);
    const double clhs286 =     clhs185 + clhs272*clhs40 + clhs273*clhs42 + clhs285*clhs3;
    const double clhs287 =     clhs2*epsilon*(clhs87 - clhs88 - clhs89);
    const double clhs288 =     clhs188 + clhs277*clhs40 + clhs278*clhs42 + clhs287*clhs3;
    const double clhs289 =     clhs2*epsilon*(clhs56 - clhs57 - clhs58);
    const double clhs290 =     clhs192 + clhs272*clhs50 + clhs273*clhs52 + clhs289*clhs3;
    const double clhs291 =     clhs2*epsilon*(clhs93 - clhs94 - clhs95);
    const double clhs292 =     clhs195 + clhs277*clhs50 + clhs278*clhs52 + clhs291*clhs3;
    const double clhs293 =     clhs271*clhs5*epsilon;
    const double clhs294 =     clhs15*clhs293 + clhs236 + clhs24*clhs273 + clhs274*clhs5;
    const double clhs295 =     clhs276*clhs5*epsilon;
    const double clhs296 =     clhs15*clhs295 + clhs24*clhs278 + clhs242 + clhs279*clhs5;
    const double clhs297 =     clhs245 + clhs273*clhs34 + clhs281*clhs5 + clhs293*clhs30;
    const double clhs298 =     clhs247 + clhs278*clhs34 + clhs283*clhs5 + clhs295*clhs30;
    const double clhs299 =     clhs249 + clhs273*clhs44 + clhs285*clhs5 + clhs293*clhs40;
    const double clhs300 =     clhs251 + clhs278*clhs44 + clhs287*clhs5 + clhs295*clhs40;
    const double clhs301 =     clhs253 + clhs273*clhs54 + clhs289*clhs5 + clhs293*clhs50;
    const double clhs302 =     clhs255 + clhs278*clhs54 + clhs291*clhs5 + clhs295*clhs50;

    lhs(0,0)=clhs1*clhs11;
    lhs(0,1)=clhs11*clhs12;
    lhs(0,2)=clhs11*clhs13;
    lhs(0,3)=clhs11*clhs14;
    lhs(0,4)=clhs15*clhs16 + clhs15*clhs19 + clhs17*clhs18 + clhs17*clhs20 + clhs21*clhs26 + clhs21*clhs29;
    lhs(0,5)=clhs16*clhs30 + clhs18*clhs31 + clhs19*clhs30 + clhs20*clhs31 + clhs21*clhs36 + clhs21*clhs39;
    lhs(0,6)=clhs16*clhs40 + clhs18*clhs41 + clhs19*clhs40 + clhs20*clhs41 + clhs21*clhs46 + clhs21*clhs49;
    lhs(0,7)=clhs16*clhs50 + clhs18*clhs51 + clhs19*clhs50 + clhs20*clhs51 + clhs21*clhs56 + clhs21*clhs59;
    lhs(0,8)=clhs60;
    lhs(0,9)=0;
    lhs(0,10)=clhs61;
    lhs(0,11)=0;
    lhs(0,12)=clhs60;
    lhs(0,13)=0;
    lhs(0,14)=clhs61;
    lhs(0,15)=0;
    lhs(1,0)=clhs1*clhs68;
    lhs(1,1)=clhs12*clhs68;
    lhs(1,2)=clhs13*clhs68;
    lhs(1,3)=clhs14*clhs68;
    lhs(1,4)=clhs15*clhs69 + clhs15*clhs71 + clhs17*clhs70 + clhs17*clhs72 + clhs21*clhs75 + clhs21*clhs78;
    lhs(1,5)=clhs21*clhs81 + clhs21*clhs84 + clhs30*clhs69 + clhs30*clhs71 + clhs31*clhs70 + clhs31*clhs72;
    lhs(1,6)=clhs21*clhs87 + clhs21*clhs90 + clhs40*clhs69 + clhs40*clhs71 + clhs41*clhs70 + clhs41*clhs72;
    lhs(1,7)=clhs21*clhs93 + clhs21*clhs96 + clhs50*clhs69 + clhs50*clhs71 + clhs51*clhs70 + clhs51*clhs72;
    lhs(1,8)=0;
    lhs(1,9)=clhs60;
    lhs(1,10)=0;
    lhs(1,11)=clhs61;
    lhs(1,12)=0;
    lhs(1,13)=clhs60;
    lhs(1,14)=0;
    lhs(1,15)=clhs61;
    lhs(2,0)=clhs11*clhs98;
    lhs(2,1)=clhs11*clhs99;
    lhs(2,2)=clhs100*clhs11;
    lhs(2,3)=clhs101*clhs11;
    lhs(2,4)=clhs102*clhs15 + clhs103*clhs18 + clhs103*clhs20 + clhs104*clhs15 + clhs105*clhs26 + clhs105*clhs29;
    lhs(2,5)=clhs102*clhs30 + clhs104*clhs30 + clhs105*clhs36 + clhs105*clhs39 + clhs106*clhs18 + clhs106*clhs20;
    lhs(2,6)=clhs102*clhs40 + clhs104*clhs40 + clhs105*clhs46 + clhs105*clhs49 + clhs107*clhs18 + clhs107*clhs20;
    lhs(2,7)=clhs102*clhs50 + clhs104*clhs50 + clhs105*clhs56 + clhs105*clhs59 + clhs108*clhs18 + clhs108*clhs20;
    lhs(2,8)=clhs109;
    lhs(2,9)=0;
    lhs(2,10)=clhs110;
    lhs(2,11)=0;
    lhs(2,12)=clhs109;
    lhs(2,13)=0;
    lhs(2,14)=clhs110;
    lhs(2,15)=0;
    lhs(3,0)=clhs68*clhs98;
    lhs(3,1)=clhs68*clhs99;
    lhs(3,2)=clhs100*clhs68;
    lhs(3,3)=clhs101*clhs68;
    lhs(3,4)=clhs103*clhs70 + clhs103*clhs72 + clhs105*clhs75 + clhs105*clhs78 + clhs111*clhs15 + clhs112*clhs15;
    lhs(3,5)=clhs105*clhs81 + clhs105*clhs84 + clhs106*clhs70 + clhs106*clhs72 + clhs111*clhs30 + clhs112*clhs30;
    lhs(3,6)=clhs105*clhs87 + clhs105*clhs90 + clhs107*clhs70 + clhs107*clhs72 + clhs111*clhs40 + clhs112*clhs40;
    lhs(3,7)=clhs105*clhs93 + clhs105*clhs96 + clhs108*clhs70 + clhs108*clhs72 + clhs111*clhs50 + clhs112*clhs50;
    lhs(3,8)=0;
    lhs(3,9)=clhs109;
    lhs(3,10)=0;
    lhs(3,11)=clhs110;
    lhs(3,12)=0;
    lhs(3,13)=clhs109;
    lhs(3,14)=0;
    lhs(3,15)=clhs110;
    lhs(4,0)=-clhs11*clhs114;
    lhs(4,1)=-clhs11*clhs115;
    lhs(4,2)=-clhs11*clhs116;
    lhs(4,3)=-clhs11*clhs117;
    lhs(4,4)=-clhs118*clhs15 - clhs119*clhs18 - clhs119*clhs20 - clhs120*clhs15 - clhs121*clhs26 - clhs121*clhs29;
    lhs(4,5)=-clhs118*clhs30 - clhs120*clhs30 - clhs121*clhs36 - clhs121*clhs39 - clhs122*clhs18 - clhs122*clhs20;
    lhs(4,6)=-clhs118*clhs40 - clhs120*clhs40 - clhs121*clhs46 - clhs121*clhs49 - clhs123*clhs18 - clhs123*clhs20;
    lhs(4,7)=-clhs118*clhs50 - clhs120*clhs50 - clhs121*clhs56 - clhs121*clhs59 - clhs124*clhs18 - clhs124*clhs20;
    lhs(4,8)=clhs125;
    lhs(4,9)=0;
    lhs(4,10)=clhs126;
    lhs(4,11)=0;
    lhs(4,12)=clhs125;
    lhs(4,13)=0;
    lhs(4,14)=clhs126;
    lhs(4,15)=0;
    lhs(5,0)=-clhs114*clhs68;
    lhs(5,1)=-clhs115*clhs68;
    lhs(5,2)=-clhs116*clhs68;
    lhs(5,3)=-clhs117*clhs68;
    lhs(5,4)=-clhs119*clhs70 - clhs119*clhs72 - clhs121*clhs75 - clhs121*clhs78 - clhs127*clhs15 - clhs128*clhs15;
    lhs(5,5)=-clhs121*clhs81 - clhs121*clhs84 - clhs122*clhs70 - clhs122*clhs72 - clhs127*clhs30 - clhs128*clhs30;
    lhs(5,6)=-clhs121*clhs87 - clhs121*clhs90 - clhs123*clhs70 - clhs123*clhs72 - clhs127*clhs40 - clhs128*clhs40;
    lhs(5,7)=-clhs121*clhs93 - clhs121*clhs96 - clhs124*clhs70 - clhs124*clhs72 - clhs127*clhs50 - clhs128*clhs50;
    lhs(5,8)=0;
    lhs(5,9)=clhs125;
    lhs(5,10)=0;
    lhs(5,11)=clhs126;
    lhs(5,12)=0;
    lhs(5,13)=clhs125;
    lhs(5,14)=0;
    lhs(5,15)=clhs126;
    lhs(6,0)=-clhs11*clhs130;
    lhs(6,1)=-clhs11*clhs131;
    lhs(6,2)=-clhs11*clhs132;
    lhs(6,3)=-clhs11*clhs133;
    lhs(6,4)=-clhs134*clhs15 - clhs135*clhs18 - clhs135*clhs20 - clhs136*clhs15 - clhs137*clhs26 - clhs137*clhs29;
    lhs(6,5)=-clhs134*clhs30 - clhs136*clhs30 - clhs137*clhs36 - clhs137*clhs39 - clhs138*clhs18 - clhs138*clhs20;
    lhs(6,6)=-clhs134*clhs40 - clhs136*clhs40 - clhs137*clhs46 - clhs137*clhs49 - clhs139*clhs18 - clhs139*clhs20;
    lhs(6,7)=-clhs134*clhs50 - clhs136*clhs50 - clhs137*clhs56 - clhs137*clhs59 - clhs140*clhs18 - clhs140*clhs20;
    lhs(6,8)=clhs141;
    lhs(6,9)=0;
    lhs(6,10)=clhs142;
    lhs(6,11)=0;
    lhs(6,12)=clhs141;
    lhs(6,13)=0;
    lhs(6,14)=clhs142;
    lhs(6,15)=0;
    lhs(7,0)=-clhs130*clhs68;
    lhs(7,1)=-clhs131*clhs68;
    lhs(7,2)=-clhs132*clhs68;
    lhs(7,3)=-clhs133*clhs68;
    lhs(7,4)=-clhs135*clhs70 - clhs135*clhs72 - clhs137*clhs75 - clhs137*clhs78 - clhs143*clhs15 - clhs144*clhs15;
    lhs(7,5)=-clhs137*clhs81 - clhs137*clhs84 - clhs138*clhs70 - clhs138*clhs72 - clhs143*clhs30 - clhs144*clhs30;
    lhs(7,6)=-clhs137*clhs87 - clhs137*clhs90 - clhs139*clhs70 - clhs139*clhs72 - clhs143*clhs40 - clhs144*clhs40;
    lhs(7,7)=-clhs137*clhs93 - clhs137*clhs96 - clhs140*clhs70 - clhs140*clhs72 - clhs143*clhs50 - clhs144*clhs50;
    lhs(7,8)=0;
    lhs(7,9)=clhs141;
    lhs(7,10)=0;
    lhs(7,11)=clhs142;
    lhs(7,12)=0;
    lhs(7,13)=clhs141;
    lhs(7,14)=0;
    lhs(7,15)=clhs142;
    lhs(8,0)=clhs150;
    lhs(8,1)=clhs152;
    lhs(8,2)=clhs154;
    lhs(8,3)=clhs156;
    lhs(8,4)=-GPnormal[0]*clhs167 - GPnormal[1]*clhs176;
    lhs(8,5)=-GPnormal[0]*clhs180 - GPnormal[1]*clhs183;
    lhs(8,6)=-GPnormal[0]*clhs187 - GPnormal[1]*clhs190;
    lhs(8,7)=-GPnormal[0]*clhs194 - GPnormal[1]*clhs197;
    lhs(8,8)=clhs201;
    lhs(8,9)=clhs204;
    lhs(8,10)=clhs207;
    lhs(8,11)=clhs209;
    lhs(8,12)=clhs200;
    lhs(8,13)=clhs203;
    lhs(8,14)=clhs206;
    lhs(8,15)=clhs208;
    lhs(9,0)=clhs211;
    lhs(9,1)=clhs212;
    lhs(9,2)=clhs213;
    lhs(9,3)=clhs214;
    lhs(9,4)=-GPtangent1[0]*clhs167 - GPtangent1[1]*clhs176;
    lhs(9,5)=-GPtangent1[0]*clhs180 - GPtangent1[1]*clhs183;
    lhs(9,6)=-GPtangent1[0]*clhs187 - GPtangent1[1]*clhs190;
    lhs(9,7)=-GPtangent1[0]*clhs194 - GPtangent1[1]*clhs197;
    lhs(9,8)=clhs217;
    lhs(9,9)=clhs220;
    lhs(9,10)=clhs222;
    lhs(9,11)=clhs224;
    lhs(9,12)=clhs216;
    lhs(9,13)=clhs219;
    lhs(9,14)=clhs221;
    lhs(9,15)=clhs223;
    lhs(10,0)=clhs228;
    lhs(10,1)=clhs229;
    lhs(10,2)=clhs230;
    lhs(10,3)=clhs231;
    lhs(10,4)=-GPnormal[0]*clhs238 - GPnormal[1]*clhs244;
    lhs(10,5)=-GPnormal[0]*clhs246 - GPnormal[1]*clhs248;
    lhs(10,6)=-GPnormal[0]*clhs250 - GPnormal[1]*clhs252;
    lhs(10,7)=-GPnormal[0]*clhs254 - GPnormal[1]*clhs256;
    lhs(10,8)=clhs207;
    lhs(10,9)=clhs209;
    lhs(10,10)=clhs259;
    lhs(10,11)=clhs261;
    lhs(10,12)=clhs206;
    lhs(10,13)=clhs208;
    lhs(10,14)=clhs258;
    lhs(10,15)=clhs260;
    lhs(11,0)=clhs263;
    lhs(11,1)=clhs264;
    lhs(11,2)=clhs265;
    lhs(11,3)=clhs266;
    lhs(11,4)=-GPtangent1[0]*clhs238 - GPtangent1[1]*clhs244;
    lhs(11,5)=-GPtangent1[0]*clhs246 - GPtangent1[1]*clhs248;
    lhs(11,6)=-GPtangent1[0]*clhs250 - GPtangent1[1]*clhs252;
    lhs(11,7)=-GPtangent1[0]*clhs254 - GPtangent1[1]*clhs256;
    lhs(11,8)=clhs222;
    lhs(11,9)=clhs224;
    lhs(11,10)=clhs268;
    lhs(11,11)=clhs270;
    lhs(11,12)=clhs221;
    lhs(11,13)=clhs223;
    lhs(11,14)=clhs267;
    lhs(11,15)=clhs269;
    lhs(12,0)=clhs150;
    lhs(12,1)=clhs152;
    lhs(12,2)=clhs154;
    lhs(12,3)=clhs156;
    lhs(12,4)=-GPnormal[0]*clhs275 - GPnormal[1]*clhs280;
    lhs(12,5)=-GPnormal[0]*clhs282 - GPnormal[1]*clhs284;
    lhs(12,6)=-GPnormal[0]*clhs286 - GPnormal[1]*clhs288;
    lhs(12,7)=-GPnormal[0]*clhs290 - GPnormal[1]*clhs292;
    lhs(12,8)=clhs200;
    lhs(12,9)=clhs203;
    lhs(12,10)=clhs206;
    lhs(12,11)=clhs208;
    lhs(12,12)=clhs201;
    lhs(12,13)=clhs204;
    lhs(12,14)=clhs207;
    lhs(12,15)=clhs209;
    lhs(13,0)=clhs211;
    lhs(13,1)=clhs212;
    lhs(13,2)=clhs213;
    lhs(13,3)=clhs214;
    lhs(13,4)=-GPtangent1[0]*clhs275 - GPtangent1[1]*clhs280;
    lhs(13,5)=-GPtangent1[0]*clhs282 - GPtangent1[1]*clhs284;
    lhs(13,6)=-GPtangent1[0]*clhs286 - GPtangent1[1]*clhs288;
    lhs(13,7)=-GPtangent1[0]*clhs290 - GPtangent1[1]*clhs292;
    lhs(13,8)=clhs216;
    lhs(13,9)=clhs219;
    lhs(13,10)=clhs221;
    lhs(13,11)=clhs223;
    lhs(13,12)=clhs217;
    lhs(13,13)=clhs220;
    lhs(13,14)=clhs222;
    lhs(13,15)=clhs224;
    lhs(14,0)=clhs228;
    lhs(14,1)=clhs229;
    lhs(14,2)=clhs230;
    lhs(14,3)=clhs231;
    lhs(14,4)=-GPnormal[0]*clhs294 - GPnormal[1]*clhs296;
    lhs(14,5)=-GPnormal[0]*clhs297 - GPnormal[1]*clhs298;
    lhs(14,6)=-GPnormal[0]*clhs299 - GPnormal[1]*clhs300;
    lhs(14,7)=-GPnormal[0]*clhs301 - GPnormal[1]*clhs302;
    lhs(14,8)=clhs206;
    lhs(14,9)=clhs208;
    lhs(14,10)=clhs258;
    lhs(14,11)=clhs260;
    lhs(14,12)=clhs207;
    lhs(14,13)=clhs209;
    lhs(14,14)=clhs259;
    lhs(14,15)=clhs261;
    lhs(15,0)=clhs263;
    lhs(15,1)=clhs264;
    lhs(15,2)=clhs265;
    lhs(15,3)=clhs266;
    lhs(15,4)=-GPtangent1[0]*clhs294 - GPtangent1[1]*clhs296;
    lhs(15,5)=-GPtangent1[0]*clhs297 - GPtangent1[1]*clhs298;
    lhs(15,6)=-GPtangent1[0]*clhs299 - GPtangent1[1]*clhs300;
    lhs(15,7)=-GPtangent1[0]*clhs301 - GPtangent1[1]*clhs302;
    lhs(15,8)=clhs221;
    lhs(15,9)=clhs223;
    lhs(15,10)=clhs267;
    lhs(15,11)=clhs269;
    lhs(15,12)=clhs222;
    lhs(15,13)=clhs224;
    lhs(15,14)=clhs268;
    lhs(15,15)=clhs270;

    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,16,16> ComputeGaussPointStickLHS(
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const double mu, 
        const ContactData<2,2>& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    bounded_matrix<double,16,16> lhs;
    
    const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
    const bounded_matrix<double, 2, 2> dlm             = rContactData.DoubleLagrangeMultipliers;
    const double Dt              = rContactData.Dt;
//     const double epsilon  = rContactData.epsilon;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
    
    const array_1d<double,2> GPnormal     = prod(trans(normalslave), N1);
    const array_1d<double,2> GPtangent1   = prod(trans(tan1slave), N1);
    
    const bounded_matrix<double, 2, 2> v1 = rContactData.v1;
    const bounded_matrix<double, 2, 2> v2 = rContactData.v2;
    
    const std::vector<double> DeltaJs = rContactData.DeltaJ_s;
    const std::vector<array_1d<double,2>> DeltaPhi = rContactData.DeltaPhi;
    const std::vector<array_1d<double,2>> DeltaN1  = rContactData.DeltaN1;
    const std::vector<array_1d<double,2>> DeltaN2  = rContactData.DeltaN2;
    const std::vector<bounded_matrix<double, 2, 2>> DeltaNormals    = rContactData.Delta_Normal_s;
    const std::vector<bounded_matrix<double, 2, 2>> Deltatangentxis = rContactData.Delta_Tangent_xi_s;
    
    const double clhs0 =     1.0/Dt;
    const double clhs1 =     Phi[0]; // PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs2 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs3 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs4 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs5 =     clhs0*clhs1*clhs2*(GPnormal[0]*clhs3 + GPnormal[1]*clhs4);
    const double clhs6 =     N1[0]; // N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs7 =     DeltaN1[4][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs8 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs9 =     N1[1]; // N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs10 =     DeltaN1[4][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs11 =     Dt*v1(0,0);
    const double clhs12 =     Dt*v1(1,0);
    const double clhs13 =     Dt*v2(0,0);
    const double clhs14 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs15 =     Dt*v2(1,0);
    const double clhs16 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs17 =     clhs11*clhs6 + clhs12*clhs9 - clhs13*clhs14 - clhs15*clhs16;
    const double clhs18 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs19 =     Dt*v1(0,1);
    const double clhs20 =     Dt*v1(1,1);
    const double clhs21 =     Dt*v2(0,1);
    const double clhs22 =     Dt*v2(1,1);
    const double clhs23 =     -clhs14*clhs21 - clhs16*clhs22 + clhs19*clhs6 + clhs20*clhs9;
    const double clhs24 =     clhs18*clhs9 + clhs4*clhs6;
    const double clhs25 =     DeltaN2[4][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs26 =     DeltaN2[4][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs27 =     clhs3*clhs6 + clhs8*clhs9;
    const double clhs28 =     clhs17*(clhs10*clhs8 + clhs3*clhs7) + clhs23*(clhs10*clhs18 + clhs4*clhs7) + clhs24*(clhs10*clhs20 + clhs19*clhs7 - clhs21*clhs25 - clhs22*clhs26) - clhs27*(-clhs10*clhs12 - clhs11*clhs7 + clhs13*clhs25 + clhs14 + clhs15*clhs26);
    const double clhs29 =     clhs28*clhs5;
    const double clhs30 =     DeltaN1[5][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs31 =     DeltaN1[5][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs32 =     DeltaN2[5][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs33 =     DeltaN2[5][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs34 =     clhs17*(clhs3*clhs30 + clhs31*clhs8) + clhs23*(clhs18*clhs31 + clhs30*clhs4) - clhs24*(clhs14 - clhs19*clhs30 - clhs20*clhs31 + clhs21*clhs32 + clhs22*clhs33) + clhs27*(clhs11*clhs30 + clhs12*clhs31 - clhs13*clhs32 - clhs15*clhs33);
    const double clhs35 =     clhs34*clhs5;
    const double clhs36 =     DeltaN1[6][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs37 =     DeltaN1[6][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs38 =     DeltaN2[6][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs39 =     DeltaN2[6][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs40 =     clhs17*(clhs3*clhs36 + clhs37*clhs8) + clhs23*(clhs18*clhs37 + clhs36*clhs4) + clhs24*(clhs19*clhs36 + clhs20*clhs37 - clhs21*clhs38 - clhs22*clhs39) - clhs27*(-clhs11*clhs36 - clhs12*clhs37 + clhs13*clhs38 + clhs15*clhs39 + clhs16);
    const double clhs41 =     clhs40*clhs5;
    const double clhs42 =     DeltaN1[7][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs43 =     DeltaN1[7][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs44 =     DeltaN2[7][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs45 =     DeltaN2[7][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs46 =     clhs17*(clhs3*clhs42 + clhs43*clhs8) + clhs23*(clhs18*clhs43 + clhs4*clhs42) - clhs24*(clhs16 - clhs19*clhs42 - clhs20*clhs43 + clhs21*clhs44 + clhs22*clhs45) + clhs27*(clhs11*clhs42 + clhs12*clhs43 - clhs13*clhs44 - clhs15*clhs45);
    const double clhs47 =     clhs46*clhs5;
    const double clhs48 =     Deltatangentxis[0](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs49 =     clhs17*clhs27 + clhs23*clhs24;
    const double clhs50 =     clhs1*clhs2*clhs49;
    const double clhs51 =     DeltaJs[0]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs52 =     clhs1*clhs3*clhs49;
    const double clhs53 =     DeltaPhi[0][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs54 =     clhs2*clhs3*clhs49;
    const double clhs55 =     DeltaN1[0][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs56 =     DeltaN1[0][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs57 =     DeltaN2[0][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs58 =     DeltaN2[0][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs59 =     Deltatangentxis[0](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs60 =     Deltatangentxis[0](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs61 =     Deltatangentxis[0](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs62 =     clhs17*(clhs3*clhs55 + clhs48*clhs6 + clhs56*clhs8 + clhs59*clhs9) + clhs23*(clhs18*clhs56 + clhs4*clhs55 + clhs6*clhs60 + clhs61*clhs9) + clhs24*(clhs19*clhs55 + clhs20*clhs56 - clhs21*clhs57 - clhs22*clhs58) + clhs27*(clhs11*clhs55 + clhs12*clhs56 - clhs13*clhs57 - clhs15*clhs58 + clhs6);
    const double clhs63 =     clhs1*clhs2*clhs62;
    const double clhs64 =     clhs3*clhs63 + clhs48*clhs50 + clhs51*clhs52 + clhs53*clhs54;
    const double clhs65 =     clhs1*clhs4*clhs49;
    const double clhs66 =     clhs2*clhs4*clhs49;
    const double clhs67 =     clhs4*clhs63 + clhs50*clhs60 + clhs51*clhs65 + clhs53*clhs66;
    const double clhs68 =     clhs0*(GPnormal[0]*clhs64 + GPnormal[1]*clhs67);
    const double clhs69 =     Deltatangentxis[1](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs70 =     DeltaJs[1]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs71 =     DeltaPhi[1][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs72 =     DeltaN1[1][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs73 =     DeltaN1[1][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs74 =     DeltaN2[1][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs75 =     DeltaN2[1][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs76 =     Deltatangentxis[1](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs77 =     Deltatangentxis[1](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs78 =     Deltatangentxis[1](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs79 =     clhs17*(clhs3*clhs72 + clhs6*clhs69 + clhs73*clhs8 + clhs76*clhs9) + clhs23*(clhs18*clhs73 + clhs4*clhs72 + clhs6*clhs77 + clhs78*clhs9) + clhs24*(clhs19*clhs72 + clhs20*clhs73 - clhs21*clhs74 - clhs22*clhs75 + clhs6) + clhs27*(clhs11*clhs72 + clhs12*clhs73 - clhs13*clhs74 - clhs15*clhs75);
    const double clhs80 =     clhs1*clhs2*clhs79;
    const double clhs81 =     clhs3*clhs80 + clhs50*clhs69 + clhs52*clhs70 + clhs54*clhs71;
    const double clhs82 =     clhs4*clhs80 + clhs50*clhs77 + clhs65*clhs70 + clhs66*clhs71;
    const double clhs83 =     clhs0*(GPnormal[0]*clhs81 + GPnormal[1]*clhs82);
    const double clhs84 =     Deltatangentxis[2](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs85 =     DeltaJs[2]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs86 =     DeltaPhi[2][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs87 =     DeltaN1[2][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs88 =     DeltaN1[2][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs89 =     DeltaN2[2][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs90 =     DeltaN2[2][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs91 =     Deltatangentxis[2](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs92 =     Deltatangentxis[2](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs93 =     Deltatangentxis[2](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs94 =     clhs17*(clhs3*clhs87 + clhs6*clhs84 + clhs8*clhs88 + clhs9*clhs91) + clhs23*(clhs18*clhs88 + clhs4*clhs87 + clhs6*clhs92 + clhs9*clhs93) + clhs24*(clhs19*clhs87 + clhs20*clhs88 - clhs21*clhs89 - clhs22*clhs90) + clhs27*(clhs11*clhs87 + clhs12*clhs88 - clhs13*clhs89 - clhs15*clhs90 + clhs9);
    const double clhs95 =     clhs1*clhs2*clhs94;
    const double clhs96 =     clhs3*clhs95 + clhs50*clhs84 + clhs52*clhs85 + clhs54*clhs86;
    const double clhs97 =     clhs4*clhs95 + clhs50*clhs92 + clhs65*clhs85 + clhs66*clhs86;
    const double clhs98 =     clhs0*(GPnormal[0]*clhs96 + GPnormal[1]*clhs97);
    const double clhs99 =     Deltatangentxis[3](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs100 =     DeltaJs[3]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs101 =     DeltaPhi[3][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs102 =     DeltaN1[3][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs103 =     DeltaN1[3][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs104 =     DeltaN2[3][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs105 =     DeltaN2[3][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs106 =     Deltatangentxis[3](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs107 =     Deltatangentxis[3](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs108 =     Deltatangentxis[3](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs109 =     clhs17*(clhs102*clhs3 + clhs103*clhs8 + clhs106*clhs9 + clhs6*clhs99) + clhs23*(clhs102*clhs4 + clhs103*clhs18 + clhs107*clhs6 + clhs108*clhs9) + clhs24*(clhs102*clhs19 + clhs103*clhs20 - clhs104*clhs21 - clhs105*clhs22 + clhs9) + clhs27*(clhs102*clhs11 + clhs103*clhs12 - clhs104*clhs13 - clhs105*clhs15);
    const double clhs110 =     clhs1*clhs109*clhs2;
    const double clhs111 =     clhs100*clhs52 + clhs101*clhs54 + clhs110*clhs3 + clhs50*clhs99;
    const double clhs112 =     clhs100*clhs65 + clhs101*clhs66 + clhs107*clhs50 + clhs110*clhs4;
    const double clhs113 =     clhs0*(GPnormal[0]*clhs111 + GPnormal[1]*clhs112);
    const double clhs114 =     clhs0*clhs1*clhs2*(GPtangent1[0]*clhs3 + GPtangent1[1]*clhs4);
    const double clhs115 =     clhs114*clhs28;
    const double clhs116 =     clhs114*clhs34;
    const double clhs117 =     clhs114*clhs40;
    const double clhs118 =     clhs114*clhs46;
    const double clhs119 =     clhs0*(GPtangent1[0]*clhs64 + GPtangent1[1]*clhs67);
    const double clhs120 =     clhs0*(GPtangent1[0]*clhs81 + GPtangent1[1]*clhs82);
    const double clhs121 =     clhs0*(GPtangent1[0]*clhs96 + GPtangent1[1]*clhs97);
    const double clhs122 =     clhs0*(GPtangent1[0]*clhs111 + GPtangent1[1]*clhs112);
    const double clhs123 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs124 =     clhs0*clhs123*clhs2*(GPnormal[0]*clhs8 + GPnormal[1]*clhs18);
    const double clhs125 =     clhs124*clhs28;
    const double clhs126 =     clhs124*clhs34;
    const double clhs127 =     clhs124*clhs40;
    const double clhs128 =     clhs124*clhs46;
    const double clhs129 =     clhs123*clhs2*clhs49;
    const double clhs130 =     clhs123*clhs49*clhs8;
    const double clhs131 =     DeltaPhi[0][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs132 =     clhs2*clhs49*clhs8;
    const double clhs133 =     clhs123*clhs2*clhs62;
    const double clhs134 =     clhs129*clhs59 + clhs130*clhs51 + clhs131*clhs132 + clhs133*clhs8;
    const double clhs135 =     clhs123*clhs18*clhs49;
    const double clhs136 =     clhs18*clhs2*clhs49;
    const double clhs137 =     clhs129*clhs61 + clhs131*clhs136 + clhs133*clhs18 + clhs135*clhs51;
    const double clhs138 =     clhs0*(GPnormal[0]*clhs134 + GPnormal[1]*clhs137);
    const double clhs139 =     DeltaPhi[1][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs140 =     clhs123*clhs2*clhs79;
    const double clhs141 =     clhs129*clhs76 + clhs130*clhs70 + clhs132*clhs139 + clhs140*clhs8;
    const double clhs142 =     clhs129*clhs78 + clhs135*clhs70 + clhs136*clhs139 + clhs140*clhs18;
    const double clhs143 =     clhs0*(GPnormal[0]*clhs141 + GPnormal[1]*clhs142);
    const double clhs144 =     DeltaPhi[2][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs145 =     clhs123*clhs2*clhs94;
    const double clhs146 =     clhs129*clhs91 + clhs130*clhs85 + clhs132*clhs144 + clhs145*clhs8;
    const double clhs147 =     clhs129*clhs93 + clhs135*clhs85 + clhs136*clhs144 + clhs145*clhs18;
    const double clhs148 =     clhs0*(GPnormal[0]*clhs146 + GPnormal[1]*clhs147);
    const double clhs149 =     DeltaPhi[3][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs150 =     clhs109*clhs123*clhs2;
    const double clhs151 =     clhs100*clhs130 + clhs106*clhs129 + clhs132*clhs149 + clhs150*clhs8;
    const double clhs152 =     clhs100*clhs135 + clhs108*clhs129 + clhs136*clhs149 + clhs150*clhs18;
    const double clhs153 =     clhs0*(GPnormal[0]*clhs151 + GPnormal[1]*clhs152);
    const double clhs154 =     clhs0*clhs123*clhs2*(GPtangent1[0]*clhs8 + GPtangent1[1]*clhs18);
    const double clhs155 =     clhs154*clhs28;
    const double clhs156 =     clhs154*clhs34;
    const double clhs157 =     clhs154*clhs40;
    const double clhs158 =     clhs154*clhs46;
    const double clhs159 =     clhs0*(GPtangent1[0]*clhs134 + GPtangent1[1]*clhs137);
    const double clhs160 =     clhs0*(GPtangent1[0]*clhs141 + GPtangent1[1]*clhs142);
    const double clhs161 =     clhs0*(GPtangent1[0]*clhs146 + GPtangent1[1]*clhs147);
    const double clhs162 =     clhs0*(GPtangent1[0]*clhs151 + GPtangent1[1]*clhs152);

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
    lhs(8,0)=clhs29;
    lhs(8,1)=clhs35;
    lhs(8,2)=clhs41;
    lhs(8,3)=clhs47;
    lhs(8,4)=clhs68;
    lhs(8,5)=clhs83;
    lhs(8,6)=clhs98;
    lhs(8,7)=clhs113;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(8,12)=0;
    lhs(8,13)=0;
    lhs(8,14)=0;
    lhs(8,15)=0;
    lhs(9,0)=clhs115;
    lhs(9,1)=clhs116;
    lhs(9,2)=clhs117;
    lhs(9,3)=clhs118;
    lhs(9,4)=clhs119;
    lhs(9,5)=clhs120;
    lhs(9,6)=clhs121;
    lhs(9,7)=clhs122;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(9,12)=0;
    lhs(9,13)=0;
    lhs(9,14)=0;
    lhs(9,15)=0;
    lhs(10,0)=clhs125;
    lhs(10,1)=clhs126;
    lhs(10,2)=clhs127;
    lhs(10,3)=clhs128;
    lhs(10,4)=clhs138;
    lhs(10,5)=clhs143;
    lhs(10,6)=clhs148;
    lhs(10,7)=clhs153;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(10,12)=0;
    lhs(10,13)=0;
    lhs(10,14)=0;
    lhs(10,15)=0;
    lhs(11,0)=clhs155;
    lhs(11,1)=clhs156;
    lhs(11,2)=clhs157;
    lhs(11,3)=clhs158;
    lhs(11,4)=clhs159;
    lhs(11,5)=clhs160;
    lhs(11,6)=clhs161;
    lhs(11,7)=clhs162;
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;
    lhs(11,12)=0;
    lhs(11,13)=0;
    lhs(11,14)=0;
    lhs(11,15)=0;
    lhs(12,0)=clhs29;
    lhs(12,1)=clhs35;
    lhs(12,2)=clhs41;
    lhs(12,3)=clhs47;
    lhs(12,4)=clhs68;
    lhs(12,5)=clhs83;
    lhs(12,6)=clhs98;
    lhs(12,7)=clhs113;
    lhs(12,8)=0;
    lhs(12,9)=0;
    lhs(12,10)=0;
    lhs(12,11)=0;
    lhs(12,12)=0;
    lhs(12,13)=0;
    lhs(12,14)=0;
    lhs(12,15)=0;
    lhs(13,0)=clhs115;
    lhs(13,1)=clhs116;
    lhs(13,2)=clhs117;
    lhs(13,3)=clhs118;
    lhs(13,4)=clhs119;
    lhs(13,5)=clhs120;
    lhs(13,6)=clhs121;
    lhs(13,7)=clhs122;
    lhs(13,8)=0;
    lhs(13,9)=0;
    lhs(13,10)=0;
    lhs(13,11)=0;
    lhs(13,12)=0;
    lhs(13,13)=0;
    lhs(13,14)=0;
    lhs(13,15)=0;
    lhs(14,0)=clhs125;
    lhs(14,1)=clhs126;
    lhs(14,2)=clhs127;
    lhs(14,3)=clhs128;
    lhs(14,4)=clhs138;
    lhs(14,5)=clhs143;
    lhs(14,6)=clhs148;
    lhs(14,7)=clhs153;
    lhs(14,8)=0;
    lhs(14,9)=0;
    lhs(14,10)=0;
    lhs(14,11)=0;
    lhs(14,12)=0;
    lhs(14,13)=0;
    lhs(14,14)=0;
    lhs(14,15)=0;
    lhs(15,0)=clhs155;
    lhs(15,1)=clhs156;
    lhs(15,2)=clhs157;
    lhs(15,3)=clhs158;
    lhs(15,4)=clhs159;
    lhs(15,5)=clhs160;
    lhs(15,6)=clhs161;
    lhs(15,7)=clhs162;
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
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const double mu, 
        const ContactData<2,2>& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    bounded_matrix<double,16,16> lhs;
    
    const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
    const double Dt              = rContactData.Dt;
//     const double epsilon  = rContactData.epsilon;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
    
    const array_1d<double,2> GPnormal     = prod(trans(normalslave), N1);
    const array_1d<double,2> GPtangent1   = prod(trans(tan1slave), N1);
    
//     const double sign_tangpress = boost::math::sign(augmented_tangent_lm);
    
    const bounded_matrix<double, 2, 2> v1 = rContactData.v1;
    const bounded_matrix<double, 2, 2> v2 = rContactData.v2;
    
    const std::vector<double> DeltaJs = rContactData.DeltaJ_s;
    const std::vector<array_1d<double,2>> DeltaPhi = rContactData.DeltaPhi;
    const std::vector<array_1d<double,2>> DeltaN1  = rContactData.DeltaN1;
    const std::vector<array_1d<double,2>> DeltaN2  = rContactData.DeltaN2;
    const std::vector<bounded_matrix<double, 2, 2>> DeltaNormals    = rContactData.Delta_Normal_s;
    const std::vector<bounded_matrix<double, 2, 2>> Deltatangentxis = rContactData.Delta_Tangent_xi_s;
    
    const double clhs0 =     1.0/Dt;
    const double clhs1 =     Phi[0]; // PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs2 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs3 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs4 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs5 =     clhs0*clhs1*clhs2*(GPnormal[0]*clhs3 + GPnormal[1]*clhs4);
    const double clhs6 =     N1[0]; // N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs7 =     DeltaN1[4][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs8 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs9 =     N1[1]; // N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs10 =     DeltaN1[4][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs11 =     Dt*v1(0,0);
    const double clhs12 =     Dt*v1(1,0);
    const double clhs13 =     Dt*v2(0,0);
    const double clhs14 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs15 =     Dt*v2(1,0);
    const double clhs16 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs17 =     clhs11*clhs6 + clhs12*clhs9 - clhs13*clhs14 - clhs15*clhs16;
    const double clhs18 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs19 =     Dt*v1(0,1);
    const double clhs20 =     Dt*v1(1,1);
    const double clhs21 =     Dt*v2(0,1);
    const double clhs22 =     Dt*v2(1,1);
    const double clhs23 =     -clhs14*clhs21 - clhs16*clhs22 + clhs19*clhs6 + clhs20*clhs9;
    const double clhs24 =     clhs18*clhs9 + clhs4*clhs6;
    const double clhs25 =     DeltaN2[4][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs26 =     DeltaN2[4][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs27 =     clhs3*clhs6 + clhs8*clhs9;
    const double clhs28 =     clhs17*(clhs10*clhs8 + clhs3*clhs7) + clhs23*(clhs10*clhs18 + clhs4*clhs7) + clhs24*(clhs10*clhs20 + clhs19*clhs7 - clhs21*clhs25 - clhs22*clhs26) - clhs27*(-clhs10*clhs12 - clhs11*clhs7 + clhs13*clhs25 + clhs14 + clhs15*clhs26);
    const double clhs29 =     clhs28*clhs5;
    const double clhs30 =     DeltaN1[5][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs31 =     DeltaN1[5][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs32 =     DeltaN2[5][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs33 =     DeltaN2[5][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs34 =     clhs17*(clhs3*clhs30 + clhs31*clhs8) + clhs23*(clhs18*clhs31 + clhs30*clhs4) - clhs24*(clhs14 - clhs19*clhs30 - clhs20*clhs31 + clhs21*clhs32 + clhs22*clhs33) + clhs27*(clhs11*clhs30 + clhs12*clhs31 - clhs13*clhs32 - clhs15*clhs33);
    const double clhs35 =     clhs34*clhs5;
    const double clhs36 =     DeltaN1[6][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs37 =     DeltaN1[6][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs38 =     DeltaN2[6][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs39 =     DeltaN2[6][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs40 =     clhs17*(clhs3*clhs36 + clhs37*clhs8) + clhs23*(clhs18*clhs37 + clhs36*clhs4) + clhs24*(clhs19*clhs36 + clhs20*clhs37 - clhs21*clhs38 - clhs22*clhs39) - clhs27*(-clhs11*clhs36 - clhs12*clhs37 + clhs13*clhs38 + clhs15*clhs39 + clhs16);
    const double clhs41 =     clhs40*clhs5;
    const double clhs42 =     DeltaN1[7][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs43 =     DeltaN1[7][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs44 =     DeltaN2[7][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs45 =     DeltaN2[7][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs46 =     clhs17*(clhs3*clhs42 + clhs43*clhs8) + clhs23*(clhs18*clhs43 + clhs4*clhs42) - clhs24*(clhs16 - clhs19*clhs42 - clhs20*clhs43 + clhs21*clhs44 + clhs22*clhs45) + clhs27*(clhs11*clhs42 + clhs12*clhs43 - clhs13*clhs44 - clhs15*clhs45);
    const double clhs47 =     clhs46*clhs5;
    const double clhs48 =     Deltatangentxis[0](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs49 =     clhs17*clhs27 + clhs23*clhs24;
    const double clhs50 =     clhs1*clhs2*clhs49;
    const double clhs51 =     DeltaJs[0]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs52 =     clhs1*clhs3*clhs49;
    const double clhs53 =     DeltaPhi[0][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs54 =     clhs2*clhs3*clhs49;
    const double clhs55 =     DeltaN1[0][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs56 =     DeltaN1[0][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs57 =     DeltaN2[0][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs58 =     DeltaN2[0][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs59 =     Deltatangentxis[0](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs60 =     Deltatangentxis[0](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs61 =     Deltatangentxis[0](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs62 =     clhs17*(clhs3*clhs55 + clhs48*clhs6 + clhs56*clhs8 + clhs59*clhs9) + clhs23*(clhs18*clhs56 + clhs4*clhs55 + clhs6*clhs60 + clhs61*clhs9) + clhs24*(clhs19*clhs55 + clhs20*clhs56 - clhs21*clhs57 - clhs22*clhs58) + clhs27*(clhs11*clhs55 + clhs12*clhs56 - clhs13*clhs57 - clhs15*clhs58 + clhs6);
    const double clhs63 =     clhs1*clhs2*clhs62;
    const double clhs64 =     clhs3*clhs63 + clhs48*clhs50 + clhs51*clhs52 + clhs53*clhs54;
    const double clhs65 =     clhs1*clhs4*clhs49;
    const double clhs66 =     clhs2*clhs4*clhs49;
    const double clhs67 =     clhs4*clhs63 + clhs50*clhs60 + clhs51*clhs65 + clhs53*clhs66;
    const double clhs68 =     clhs0*(GPnormal[0]*clhs64 + GPnormal[1]*clhs67);
    const double clhs69 =     Deltatangentxis[1](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs70 =     DeltaJs[1]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs71 =     DeltaPhi[1][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs72 =     DeltaN1[1][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs73 =     DeltaN1[1][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs74 =     DeltaN2[1][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs75 =     DeltaN2[1][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs76 =     Deltatangentxis[1](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs77 =     Deltatangentxis[1](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs78 =     Deltatangentxis[1](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs79 =     clhs17*(clhs3*clhs72 + clhs6*clhs69 + clhs73*clhs8 + clhs76*clhs9) + clhs23*(clhs18*clhs73 + clhs4*clhs72 + clhs6*clhs77 + clhs78*clhs9) + clhs24*(clhs19*clhs72 + clhs20*clhs73 - clhs21*clhs74 - clhs22*clhs75 + clhs6) + clhs27*(clhs11*clhs72 + clhs12*clhs73 - clhs13*clhs74 - clhs15*clhs75);
    const double clhs80 =     clhs1*clhs2*clhs79;
    const double clhs81 =     clhs3*clhs80 + clhs50*clhs69 + clhs52*clhs70 + clhs54*clhs71;
    const double clhs82 =     clhs4*clhs80 + clhs50*clhs77 + clhs65*clhs70 + clhs66*clhs71;
    const double clhs83 =     clhs0*(GPnormal[0]*clhs81 + GPnormal[1]*clhs82);
    const double clhs84 =     Deltatangentxis[2](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs85 =     DeltaJs[2]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs86 =     DeltaPhi[2][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs87 =     DeltaN1[2][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs88 =     DeltaN1[2][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs89 =     DeltaN2[2][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs90 =     DeltaN2[2][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs91 =     Deltatangentxis[2](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs92 =     Deltatangentxis[2](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs93 =     Deltatangentxis[2](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs94 =     clhs17*(clhs3*clhs87 + clhs6*clhs84 + clhs8*clhs88 + clhs9*clhs91) + clhs23*(clhs18*clhs88 + clhs4*clhs87 + clhs6*clhs92 + clhs9*clhs93) + clhs24*(clhs19*clhs87 + clhs20*clhs88 - clhs21*clhs89 - clhs22*clhs90) + clhs27*(clhs11*clhs87 + clhs12*clhs88 - clhs13*clhs89 - clhs15*clhs90 + clhs9);
    const double clhs95 =     clhs1*clhs2*clhs94;
    const double clhs96 =     clhs3*clhs95 + clhs50*clhs84 + clhs52*clhs85 + clhs54*clhs86;
    const double clhs97 =     clhs4*clhs95 + clhs50*clhs92 + clhs65*clhs85 + clhs66*clhs86;
    const double clhs98 =     clhs0*(GPnormal[0]*clhs96 + GPnormal[1]*clhs97);
    const double clhs99 =     Deltatangentxis[3](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs100 =     DeltaJs[3]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs101 =     DeltaPhi[3][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs102 =     DeltaN1[3][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs103 =     DeltaN1[3][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs104 =     DeltaN2[3][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs105 =     DeltaN2[3][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs106 =     Deltatangentxis[3](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs107 =     Deltatangentxis[3](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs108 =     Deltatangentxis[3](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs109 =     clhs17*(clhs102*clhs3 + clhs103*clhs8 + clhs106*clhs9 + clhs6*clhs99) + clhs23*(clhs102*clhs4 + clhs103*clhs18 + clhs107*clhs6 + clhs108*clhs9) + clhs24*(clhs102*clhs19 + clhs103*clhs20 - clhs104*clhs21 - clhs105*clhs22 + clhs9) + clhs27*(clhs102*clhs11 + clhs103*clhs12 - clhs104*clhs13 - clhs105*clhs15);
    const double clhs110 =     clhs1*clhs109*clhs2;
    const double clhs111 =     clhs100*clhs52 + clhs101*clhs54 + clhs110*clhs3 + clhs50*clhs99;
    const double clhs112 =     clhs100*clhs65 + clhs101*clhs66 + clhs107*clhs50 + clhs110*clhs4;
    const double clhs113 =     clhs0*(GPnormal[0]*clhs111 + GPnormal[1]*clhs112);
    const double clhs114 =     clhs0*clhs1*clhs2*(GPtangent1[0]*clhs3 + GPtangent1[1]*clhs4);
    const double clhs115 =     clhs114*clhs28;
    const double clhs116 =     clhs114*clhs34;
    const double clhs117 =     clhs114*clhs40;
    const double clhs118 =     clhs114*clhs46;
    const double clhs119 =     clhs0*(GPtangent1[0]*clhs64 + GPtangent1[1]*clhs67);
    const double clhs120 =     clhs0*(GPtangent1[0]*clhs81 + GPtangent1[1]*clhs82);
    const double clhs121 =     clhs0*(GPtangent1[0]*clhs96 + GPtangent1[1]*clhs97);
    const double clhs122 =     clhs0*(GPtangent1[0]*clhs111 + GPtangent1[1]*clhs112);
    const double clhs123 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs124 =     clhs0*clhs123*clhs2*(GPnormal[0]*clhs8 + GPnormal[1]*clhs18);
    const double clhs125 =     clhs124*clhs28;
    const double clhs126 =     clhs124*clhs34;
    const double clhs127 =     clhs124*clhs40;
    const double clhs128 =     clhs124*clhs46;
    const double clhs129 =     clhs123*clhs2*clhs49;
    const double clhs130 =     clhs123*clhs49*clhs8;
    const double clhs131 =     DeltaPhi[0][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs132 =     clhs2*clhs49*clhs8;
    const double clhs133 =     clhs123*clhs2*clhs62;
    const double clhs134 =     clhs129*clhs59 + clhs130*clhs51 + clhs131*clhs132 + clhs133*clhs8;
    const double clhs135 =     clhs123*clhs18*clhs49;
    const double clhs136 =     clhs18*clhs2*clhs49;
    const double clhs137 =     clhs129*clhs61 + clhs131*clhs136 + clhs133*clhs18 + clhs135*clhs51;
    const double clhs138 =     clhs0*(GPnormal[0]*clhs134 + GPnormal[1]*clhs137);
    const double clhs139 =     DeltaPhi[1][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs140 =     clhs123*clhs2*clhs79;
    const double clhs141 =     clhs129*clhs76 + clhs130*clhs70 + clhs132*clhs139 + clhs140*clhs8;
    const double clhs142 =     clhs129*clhs78 + clhs135*clhs70 + clhs136*clhs139 + clhs140*clhs18;
    const double clhs143 =     clhs0*(GPnormal[0]*clhs141 + GPnormal[1]*clhs142);
    const double clhs144 =     DeltaPhi[2][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs145 =     clhs123*clhs2*clhs94;
    const double clhs146 =     clhs129*clhs91 + clhs130*clhs85 + clhs132*clhs144 + clhs145*clhs8;
    const double clhs147 =     clhs129*clhs93 + clhs135*clhs85 + clhs136*clhs144 + clhs145*clhs18;
    const double clhs148 =     clhs0*(GPnormal[0]*clhs146 + GPnormal[1]*clhs147);
    const double clhs149 =     DeltaPhi[3][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs150 =     clhs109*clhs123*clhs2;
    const double clhs151 =     clhs100*clhs130 + clhs106*clhs129 + clhs132*clhs149 + clhs150*clhs8;
    const double clhs152 =     clhs100*clhs135 + clhs108*clhs129 + clhs136*clhs149 + clhs150*clhs18;
    const double clhs153 =     clhs0*(GPnormal[0]*clhs151 + GPnormal[1]*clhs152);
    const double clhs154 =     clhs0*clhs123*clhs2*(GPtangent1[0]*clhs8 + GPtangent1[1]*clhs18);
    const double clhs155 =     clhs154*clhs28;
    const double clhs156 =     clhs154*clhs34;
    const double clhs157 =     clhs154*clhs40;
    const double clhs158 =     clhs154*clhs46;
    const double clhs159 =     clhs0*(GPtangent1[0]*clhs134 + GPtangent1[1]*clhs137);
    const double clhs160 =     clhs0*(GPtangent1[0]*clhs141 + GPtangent1[1]*clhs142);
    const double clhs161 =     clhs0*(GPtangent1[0]*clhs146 + GPtangent1[1]*clhs147);
    const double clhs162 =     clhs0*(GPtangent1[0]*clhs151 + GPtangent1[1]*clhs152);

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
    lhs(8,0)=clhs29;
    lhs(8,1)=clhs35;
    lhs(8,2)=clhs41;
    lhs(8,3)=clhs47;
    lhs(8,4)=clhs68;
    lhs(8,5)=clhs83;
    lhs(8,6)=clhs98;
    lhs(8,7)=clhs113;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(8,12)=0;
    lhs(8,13)=0;
    lhs(8,14)=0;
    lhs(8,15)=0;
    lhs(9,0)=clhs115;
    lhs(9,1)=clhs116;
    lhs(9,2)=clhs117;
    lhs(9,3)=clhs118;
    lhs(9,4)=clhs119;
    lhs(9,5)=clhs120;
    lhs(9,6)=clhs121;
    lhs(9,7)=clhs122;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(9,12)=0;
    lhs(9,13)=0;
    lhs(9,14)=0;
    lhs(9,15)=0;
    lhs(10,0)=clhs125;
    lhs(10,1)=clhs126;
    lhs(10,2)=clhs127;
    lhs(10,3)=clhs128;
    lhs(10,4)=clhs138;
    lhs(10,5)=clhs143;
    lhs(10,6)=clhs148;
    lhs(10,7)=clhs153;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(10,12)=0;
    lhs(10,13)=0;
    lhs(10,14)=0;
    lhs(10,15)=0;
    lhs(11,0)=clhs155;
    lhs(11,1)=clhs156;
    lhs(11,2)=clhs157;
    lhs(11,3)=clhs158;
    lhs(11,4)=clhs159;
    lhs(11,5)=clhs160;
    lhs(11,6)=clhs161;
    lhs(11,7)=clhs162;
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;
    lhs(11,12)=0;
    lhs(11,13)=0;
    lhs(11,14)=0;
    lhs(11,15)=0;
    lhs(12,0)=clhs29;
    lhs(12,1)=clhs35;
    lhs(12,2)=clhs41;
    lhs(12,3)=clhs47;
    lhs(12,4)=clhs68;
    lhs(12,5)=clhs83;
    lhs(12,6)=clhs98;
    lhs(12,7)=clhs113;
    lhs(12,8)=0;
    lhs(12,9)=0;
    lhs(12,10)=0;
    lhs(12,11)=0;
    lhs(12,12)=0;
    lhs(12,13)=0;
    lhs(12,14)=0;
    lhs(12,15)=0;
    lhs(13,0)=clhs115;
    lhs(13,1)=clhs116;
    lhs(13,2)=clhs117;
    lhs(13,3)=clhs118;
    lhs(13,4)=clhs119;
    lhs(13,5)=clhs120;
    lhs(13,6)=clhs121;
    lhs(13,7)=clhs122;
    lhs(13,8)=0;
    lhs(13,9)=0;
    lhs(13,10)=0;
    lhs(13,11)=0;
    lhs(13,12)=0;
    lhs(13,13)=0;
    lhs(13,14)=0;
    lhs(13,15)=0;
    lhs(14,0)=clhs125;
    lhs(14,1)=clhs126;
    lhs(14,2)=clhs127;
    lhs(14,3)=clhs128;
    lhs(14,4)=clhs138;
    lhs(14,5)=clhs143;
    lhs(14,6)=clhs148;
    lhs(14,7)=clhs153;
    lhs(14,8)=0;
    lhs(14,9)=0;
    lhs(14,10)=0;
    lhs(14,11)=0;
    lhs(14,12)=0;
    lhs(14,13)=0;
    lhs(14,14)=0;
    lhs(14,15)=0;
    lhs(15,0)=clhs155;
    lhs(15,1)=clhs156;
    lhs(15,2)=clhs157;
    lhs(15,3)=clhs158;
    lhs(15,4)=clhs159;
    lhs(15,5)=clhs160;
    lhs(15,6)=clhs161;
    lhs(15,7)=clhs162;
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
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
//         const bounded_matrix<double, 2, 2> DPhi, 
        const double detJ, 
        const ContactData<2,2>& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    bounded_matrix<double,16,16> lhs;
    
//     const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
//     const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
//     const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
//     const double Dt              = rContactData.Dt;
//     const double epsilon  = rContactData.epsilon;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
//     
//     const array_1d<double,2> GPnormal     = prod(trans(normalslave), N1);
//     const array_1d<double,2> GPtangent1   = prod(trans(tan1slave), N1);
//     
//     const bounded_matrix<double, 2, 2> v1 = rContactData.v1;
//     const bounded_matrix<double, 2, 2> v2 = rContactData.v2;
// 
//     const std::vector<double> DeltaJs  = rContactData.DeltaJ_s;
//     const std::vector<double> DeltaGap = rContactData.DeltaGap;
//     const std::vector<array_1d<double,2>> DeltaPhi = rContactData.DeltaPhi;
//     const std::vector<array_1d<double,2>> DeltaN1  = rContactData.DeltaN1;
//     const std::vector<array_1d<double,2>> DeltaN2  = rContactData.DeltaN2;
//     const std::vector<bounded_matrix<double, 2, 2>> DeltaNormals    = rContactData.Delta_Normal_s;
//     const std::vector<bounded_matrix<double, 2, 2>> Deltatangentxis = rContactData.Delta_Tangent_xi_s;
//
//substitute_inactive_lhs
    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,16> ComputeGaussPointActiveRHS(
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const ContactData<2,2>& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    array_1d<double,16> rhs;
    
    const bounded_matrix<double, 2, 2> normalslave    = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave      = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm             = rContactData.LagrangeMultipliers;
    const bounded_matrix<double, 2, 2> dlm             = rContactData.DoubleLagrangeMultipliers;
//     const double Dt             = rContactData.Dt;
    const double epsilon = rContactData.epsilon;
    
    const array_1d<double,2> GPnormal     = prod(trans(normalslave), N1);
    const array_1d<double,2> GPtangent1   = prod(trans(tan1slave), N1);
    
    const bounded_matrix<double, 2, 2> v1 = rContactData.v1;
    const bounded_matrix<double, 2, 2> v2 = rContactData.v2;
    
    const double crhs0 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs1 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs2 =     Phi[0]; // PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs3 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs4 =     crhs2*dlm(0,0) + crhs3*dlm(1,0);
    const double crhs5 =     crhs2*lm(0,0);
    const double crhs6 =     crhs3*lm(1,0);
    const double crhs7 =     crhs1*(crhs4 + crhs5 + crhs6);
    const double crhs8 =     crhs2*dlm(0,1) + crhs3*dlm(1,1);
    const double crhs9 =     crhs2*lm(0,1);
    const double crhs10 =     crhs3*lm(1,1);
    const double crhs11 =     crhs1*(crhs10 + crhs8 + crhs9);
    const double crhs12 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs13 =     N1[0]; // N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs14 =     N1[1]; // N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs15 =     crhs1*crhs2;
    const double crhs16 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs17 =     crhs16*normalslave(0,0); // CRHS16*NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs18 =     epsilon*(crhs4 - crhs5 - crhs6);
    const double crhs19 =     -crhs17 + crhs18;
    const double crhs20 =     crhs16*normalslave(0,1); // CRHS16*NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs21 =     epsilon*(-crhs10 + crhs8 - crhs9);
    const double crhs22 =     -crhs20 + crhs21;
    const double crhs23 =     crhs1*crhs3;
    const double crhs24 =     crhs16*normalslave(1,0); // CRHS16*NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs25 =     crhs18 - crhs24;
    const double crhs26 =     crhs16*normalslave(1,1); // CRHS16*NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs27 =     crhs21 - crhs26;
    const double crhs28 =     crhs17 + crhs18;
    const double crhs29 =     crhs20 + crhs21;
    const double crhs30 =     crhs18 + crhs24;
    const double crhs31 =     crhs21 + crhs26;

    rhs[0]=-crhs0*crhs7;
    rhs[1]=-crhs0*crhs11;
    rhs[2]=-crhs12*crhs7;
    rhs[3]=-crhs11*crhs12;
    rhs[4]=crhs13*crhs7;
    rhs[5]=crhs11*crhs13;
    rhs[6]=crhs14*crhs7;
    rhs[7]=crhs11*crhs14;
    rhs[8]=-crhs15*(GPnormal[0]*crhs19 + GPnormal[1]*crhs22);
    rhs[9]=-crhs15*(GPtangent1[0]*crhs19 + GPtangent1[1]*crhs22);
    rhs[10]=-crhs23*(GPnormal[0]*crhs25 + GPnormal[1]*crhs27);
    rhs[11]=-crhs23*(GPtangent1[0]*crhs25 + GPtangent1[1]*crhs27);
    rhs[12]=crhs15*(GPnormal[0]*crhs28 + GPnormal[1]*crhs29);
    rhs[13]=crhs15*(GPtangent1[0]*crhs28 + GPtangent1[1]*crhs29);
    rhs[14]=crhs23*(GPnormal[0]*crhs30 + GPnormal[1]*crhs31);
    rhs[15]=crhs23*(GPtangent1[0]*crhs30 + GPtangent1[1]*crhs31);

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,16> ComputeGaussPointStickRHS(
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const double mu, 
        const ContactData<2,2>& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    array_1d<double,16> rhs;
    
    const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
    const double Dt              = rContactData.Dt;
//     const double epsilon  = rContactData.epsilon;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
    
    const array_1d<double,2> GPnormal     = prod(trans(normalslave), N1);
    const array_1d<double,2> GPtangent1   = prod(trans(tan1slave), N1);
    
    const bounded_matrix<double, 2, 2> v1 = rContactData.v1;
    const bounded_matrix<double, 2, 2> v2 = rContactData.v2;
    
    const double crhs0 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs1 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs2 =     1.0/Dt;
    const double crhs3 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs4 =     N1[0]; // N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs5 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs6 =     N1[1]; // N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs7 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs8 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs9 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs10 =     (crhs0*crhs4 + crhs5*crhs6)*(crhs4*(Dt*v1(0,0)) + crhs6*(Dt*v1(1,0)) - crhs7*(Dt*v2(0,0)) - crhs8*(Dt*v2(1,0))) + (crhs1*crhs4 + crhs6*crhs9)*(crhs4*(Dt*v1(0,1)) + crhs6*(Dt*v1(1,1)) - crhs7*(Dt*v2(0,1)) - crhs8*(Dt*v2(1,1)));
    const double crhs11 =     crhs10*crhs2*crhs3*Phi[0]; // CRHS10*CRHS2*CRHS3*PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs12 =     -crhs11*(GPnormal[0]*crhs0 + GPnormal[1]*crhs1);
    const double crhs13 =     -crhs11*(GPtangent1[0]*crhs0 + GPtangent1[1]*crhs1);
    const double crhs14 =     crhs10*crhs2*crhs3*Phi[1]; // CRHS10*CRHS2*CRHS3*PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs15 =     -crhs14*(GPnormal[0]*crhs5 + GPnormal[1]*crhs9);
    const double crhs16 =     -crhs14*(GPtangent1[0]*crhs5 + GPtangent1[1]*crhs9);

    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=crhs12;
    rhs[9]=crhs13;
    rhs[10]=crhs15;
    rhs[11]=crhs16;
    rhs[12]=crhs12;
    rhs[13]=crhs13;
    rhs[14]=crhs15;
    rhs[15]=crhs16;

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,16> ComputeGaussPointSlipRHS(
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const double mu, 
        const ContactData<2,2>& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    array_1d<double,16> rhs;
    
    const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
    const double Dt              = rContactData.Dt;
//     const double epsilon  = rContactData.epsilon;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
    
    const array_1d<double,2> GPnormal     = prod(trans(normalslave), N1);
    const array_1d<double,2> GPtangent1   = prod(trans(tan1slave), N1);
    
//     const double sign_tangpress = boost::math::sign(augmented_tangent_lm);
    
    const bounded_matrix<double, 2, 2> v1 = rContactData.v1;
    const bounded_matrix<double, 2, 2> v2 = rContactData.v2;
    
    const double crhs0 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs1 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs2 =     1.0/Dt;
    const double crhs3 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs4 =     N1[0]; // N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs5 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs6 =     N1[1]; // N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs7 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs8 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs9 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs10 =     (crhs0*crhs4 + crhs5*crhs6)*(crhs4*(Dt*v1(0,0)) + crhs6*(Dt*v1(1,0)) - crhs7*(Dt*v2(0,0)) - crhs8*(Dt*v2(1,0))) + (crhs1*crhs4 + crhs6*crhs9)*(crhs4*(Dt*v1(0,1)) + crhs6*(Dt*v1(1,1)) - crhs7*(Dt*v2(0,1)) - crhs8*(Dt*v2(1,1)));
    const double crhs11 =     crhs10*crhs2*crhs3*Phi[0]; // CRHS10*CRHS2*CRHS3*PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs12 =     -crhs11*(GPnormal[0]*crhs0 + GPnormal[1]*crhs1);
    const double crhs13 =     -crhs11*(GPtangent1[0]*crhs0 + GPtangent1[1]*crhs1);
    const double crhs14 =     crhs10*crhs2*crhs3*Phi[1]; // CRHS10*CRHS2*CRHS3*PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs15 =     -crhs14*(GPnormal[0]*crhs5 + GPnormal[1]*crhs9);
    const double crhs16 =     -crhs14*(GPtangent1[0]*crhs5 + GPtangent1[1]*crhs9);

    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=crhs12;
    rhs[9]=crhs13;
    rhs[10]=crhs15;
    rhs[11]=crhs16;
    rhs[12]=crhs12;
    rhs[13]=crhs13;
    rhs[14]=crhs15;
    rhs[15]=crhs16;

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,16> ComputeGaussPointInactiveRHS(
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const ContactData<2,2>& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    array_1d<double,16> rhs;
    
//     const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
//     const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
//     const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
//     const double Dt              = rContactData.Dt;
//     const double epsilon  = rContactData.epsilon;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
//     
//     const array_1d<double,2> GPnormal     = prod(trans(normalslave), N1);
//     const array_1d<double,2> GPtangent1   = prod(trans(tan1slave), N1);
//     
//     const bounded_matrix<double, 2, 2> v1 = rContactData.v1;
//     const bounded_matrix<double, 2, 2> v2 = rContactData.v2;
//     
//substitute_inactive_rhs
    
    return rhs;
}
private:
};// class Contact2D2N2NDLM
}
#endif /* KRATOS_CONTACT2D2N2NDLM defined */
