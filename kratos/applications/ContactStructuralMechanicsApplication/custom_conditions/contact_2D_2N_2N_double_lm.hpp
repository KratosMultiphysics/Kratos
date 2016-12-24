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
    const double epsilon  = rContactData.epsilon;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
    
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
    const double clhs4 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs5 =     clhs3*dlm(0,0) + clhs4*dlm(1,0);
    const double clhs6 =     clhs3*lm(0,0);
    const double clhs7 =     clhs4*lm(1,0);
    const double clhs8 =     clhs6 + clhs7;
    const double clhs9 =     clhs2*(clhs5 + clhs8);
    const double clhs10 =     DeltaN2[5][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs11 =     DeltaN2[6][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs12 =     DeltaN2[7][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs13 =     DeltaJs[0]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs14 =     clhs0*clhs5;
    const double clhs15 =     DeltaN2[0][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs16 =     clhs2*clhs5;
    const double clhs17 =     clhs0*clhs8;
    const double clhs18 =     clhs2*clhs8;
    const double clhs19 =     clhs0*clhs2;
    const double clhs20 =     DeltaPhi[0][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs21 =     DeltaPhi[0][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs22 =     clhs20*dlm(0,0) + clhs21*dlm(1,0);
    const double clhs23 =     clhs20*lm(0,0);
    const double clhs24 =     clhs21*lm(1,0);
    const double clhs25 =     clhs23 + clhs24;
    const double clhs26 =     DeltaJs[1]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs27 =     DeltaN2[1][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs28 =     DeltaPhi[1][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs29 =     DeltaPhi[1][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs30 =     clhs28*dlm(0,0) + clhs29*dlm(1,0);
    const double clhs31 =     clhs28*lm(0,0);
    const double clhs32 =     clhs29*lm(1,0);
    const double clhs33 =     clhs31 + clhs32;
    const double clhs34 =     DeltaJs[2]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs35 =     DeltaN2[2][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs36 =     DeltaPhi[2][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs37 =     DeltaPhi[2][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs38 =     clhs36*dlm(0,0) + clhs37*dlm(1,0);
    const double clhs39 =     clhs36*lm(0,0);
    const double clhs40 =     clhs37*lm(1,0);
    const double clhs41 =     clhs39 + clhs40;
    const double clhs42 =     DeltaJs[3]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs43 =     DeltaN2[3][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs44 =     DeltaPhi[3][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs45 =     DeltaPhi[3][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs46 =     clhs44*dlm(0,0) + clhs45*dlm(1,0);
    const double clhs47 =     clhs44*lm(0,0);
    const double clhs48 =     clhs45*lm(1,0);
    const double clhs49 =     clhs47 + clhs48;
    const double clhs50 =     clhs19*clhs3;
    const double clhs51 =     clhs19*clhs4;
    const double clhs52 =     clhs3*dlm(0,1) + clhs4*dlm(1,1);
    const double clhs53 =     clhs3*lm(0,1);
    const double clhs54 =     clhs4*lm(1,1);
    const double clhs55 =     clhs53 + clhs54;
    const double clhs56 =     clhs2*(clhs52 + clhs55);
    const double clhs57 =     clhs0*clhs52;
    const double clhs58 =     clhs2*clhs52;
    const double clhs59 =     clhs0*clhs55;
    const double clhs60 =     clhs2*clhs55;
    const double clhs61 =     clhs20*dlm(0,1) + clhs21*dlm(1,1);
    const double clhs62 =     clhs20*lm(0,1);
    const double clhs63 =     clhs21*lm(1,1);
    const double clhs64 =     clhs62 + clhs63;
    const double clhs65 =     clhs28*dlm(0,1) + clhs29*dlm(1,1);
    const double clhs66 =     clhs28*lm(0,1);
    const double clhs67 =     clhs29*lm(1,1);
    const double clhs68 =     clhs66 + clhs67;
    const double clhs69 =     clhs36*dlm(0,1) + clhs37*dlm(1,1);
    const double clhs70 =     clhs36*lm(0,1);
    const double clhs71 =     clhs37*lm(1,1);
    const double clhs72 =     clhs70 + clhs71;
    const double clhs73 =     clhs44*dlm(0,1) + clhs45*dlm(1,1);
    const double clhs74 =     clhs44*lm(0,1);
    const double clhs75 =     clhs45*lm(1,1);
    const double clhs76 =     clhs74 + clhs75;
    const double clhs77 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs78 =     DeltaN2[4][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs79 =     DeltaN2[5][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs80 =     DeltaN2[6][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs81 =     DeltaN2[7][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs82 =     clhs5*clhs77;
    const double clhs83 =     DeltaN2[0][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs84 =     clhs77*clhs8;
    const double clhs85 =     clhs2*clhs77;
    const double clhs86 =     DeltaN2[1][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs87 =     DeltaN2[2][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs88 =     DeltaN2[3][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs89 =     clhs3*clhs85;
    const double clhs90 =     clhs4*clhs85;
    const double clhs91 =     clhs52*clhs77;
    const double clhs92 =     clhs55*clhs77;
    const double clhs93 =     N1[0]; // N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs94 =     DeltaN1[4][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs95 =     DeltaN1[5][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs96 =     DeltaN1[6][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs97 =     DeltaN1[7][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs98 =     clhs5*clhs93;
    const double clhs99 =     DeltaN1[0][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs100 =     clhs8*clhs93;
    const double clhs101 =     clhs2*clhs93;
    const double clhs102 =     DeltaN1[1][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs103 =     DeltaN1[2][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs104 =     DeltaN1[3][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs105 =     -clhs101*clhs3;
    const double clhs106 =     -clhs101*clhs4;
    const double clhs107 =     clhs52*clhs93;
    const double clhs108 =     clhs55*clhs93;
    const double clhs109 =     N1[1]; // N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs110 =     DeltaN1[4][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs111 =     DeltaN1[5][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs112 =     DeltaN1[6][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs113 =     DeltaN1[7][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs114 =     clhs109*clhs5;
    const double clhs115 =     DeltaN1[0][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs116 =     clhs109*clhs8;
    const double clhs117 =     clhs109*clhs2;
    const double clhs118 =     DeltaN1[1][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs119 =     DeltaN1[2][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs120 =     DeltaN1[3][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs121 =     -clhs117*clhs3;
    const double clhs122 =     -clhs117*clhs4;
    const double clhs123 =     clhs109*clhs52;
    const double clhs124 =     clhs109*clhs55;
    const double clhs125 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs126 =     DeltaGap[4]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs127 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs128 =     clhs127*clhs2*clhs3;
    const double clhs129 =     -clhs126*clhs128;
    const double clhs130 =     DeltaGap[5]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs131 =     -clhs128*clhs130;
    const double clhs132 =     DeltaGap[6]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs133 =     -clhs128*clhs132;
    const double clhs134 =     DeltaGap[7]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs135 =     -clhs128*clhs134;
    const double clhs136 =     clhs125*clhs2*clhs3;
    const double clhs137 =     clhs136*DeltaNormals[0](0,0); // CLHS136*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs138 =     DeltaGap[0]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs139 =     clhs128*clhs138;
    const double clhs140 =     clhs125*clhs127*clhs3;
    const double clhs141 =     clhs13*clhs140;
    const double clhs142 =     clhs125*clhs127*clhs2;
    const double clhs143 =     clhs142*clhs20;
    const double clhs144 =     clhs5 - clhs6 - clhs7;
    const double clhs145 =     clhs144*clhs3*epsilon;
    const double clhs146 =     clhs144*clhs2*epsilon;
    const double clhs147 =     clhs2*clhs3*epsilon;
    const double clhs148 =     clhs22 - clhs23 - clhs24;
    const double clhs149 =     clhs13*clhs145 + clhs146*clhs20 + clhs147*clhs148;
    const double clhs150 =     clhs136*DeltaNormals[1](0,0); // CLHS136*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs151 =     DeltaGap[1]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs152 =     clhs128*clhs151;
    const double clhs153 =     clhs140*clhs26;
    const double clhs154 =     clhs142*clhs28;
    const double clhs155 =     clhs30 - clhs31 - clhs32;
    const double clhs156 =     clhs145*clhs26 + clhs146*clhs28 + clhs147*clhs155;
    const double clhs157 =     clhs136*DeltaNormals[2](0,0); // CLHS136*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs158 =     DeltaGap[2]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs159 =     clhs128*clhs158;
    const double clhs160 =     clhs140*clhs34;
    const double clhs161 =     clhs142*clhs36;
    const double clhs162 =     clhs38 - clhs39 - clhs40;
    const double clhs163 =     clhs145*clhs34 + clhs146*clhs36 + clhs147*clhs162;
    const double clhs164 =     clhs136*DeltaNormals[3](0,0); // CLHS136*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs165 =     DeltaGap[3]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs166 =     clhs128*clhs165;
    const double clhs167 =     clhs140*clhs42;
    const double clhs168 =     clhs142*clhs44;
    const double clhs169 =     clhs46 - clhs47 - clhs48;
    const double clhs170 =     clhs145*clhs42 + clhs146*clhs44 + clhs147*clhs169;
    const double clhs171 =     clhs2*epsilon;
    const double clhs172 =     clhs171*std::pow(clhs3, 2);
    const double clhs173 =     -clhs172;
    const double clhs174 =     clhs147*clhs4;
    const double clhs175 =     -clhs174;
    const double clhs176 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs177 =     clhs176*clhs2*clhs3;
    const double clhs178 =     -clhs126*clhs177;
    const double clhs179 =     -clhs130*clhs177;
    const double clhs180 =     -clhs132*clhs177;
    const double clhs181 =     -clhs134*clhs177;
    const double clhs182 =     clhs136*DeltaNormals[0](0,1); // CLHS136*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs183 =     clhs138*clhs177;
    const double clhs184 =     clhs125*clhs176*clhs3;
    const double clhs185 =     clhs13*clhs184;
    const double clhs186 =     clhs125*clhs176*clhs2;
    const double clhs187 =     clhs186*clhs20;
    const double clhs188 =     clhs52 - clhs53 - clhs54;
    const double clhs189 =     clhs188*clhs3*epsilon;
    const double clhs190 =     clhs188*clhs2*epsilon;
    const double clhs191 =     clhs61 - clhs62 - clhs63;
    const double clhs192 =     clhs13*clhs189 + clhs147*clhs191 + clhs190*clhs20;
    const double clhs193 =     clhs136*DeltaNormals[1](0,1); // CLHS136*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs194 =     clhs151*clhs177;
    const double clhs195 =     clhs184*clhs26;
    const double clhs196 =     clhs186*clhs28;
    const double clhs197 =     clhs65 - clhs66 - clhs67;
    const double clhs198 =     clhs147*clhs197 + clhs189*clhs26 + clhs190*clhs28;
    const double clhs199 =     clhs136*DeltaNormals[2](0,1); // CLHS136*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs200 =     clhs158*clhs177;
    const double clhs201 =     clhs184*clhs34;
    const double clhs202 =     clhs186*clhs36;
    const double clhs203 =     clhs69 - clhs70 - clhs71;
    const double clhs204 =     clhs147*clhs203 + clhs189*clhs34 + clhs190*clhs36;
    const double clhs205 =     clhs136*DeltaNormals[3](0,1); // CLHS136*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs206 =     clhs165*clhs177;
    const double clhs207 =     clhs184*clhs42;
    const double clhs208 =     clhs186*clhs44;
    const double clhs209 =     clhs73 - clhs74 - clhs75;
    const double clhs210 =     clhs147*clhs209 + clhs189*clhs42 + clhs190*clhs44;
    const double clhs211 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs212 =     clhs2*clhs211*clhs4;
    const double clhs213 =     -clhs126*clhs212;
    const double clhs214 =     -clhs130*clhs212;
    const double clhs215 =     -clhs132*clhs212;
    const double clhs216 =     -clhs134*clhs212;
    const double clhs217 =     clhs125*clhs2*clhs4;
    const double clhs218 =     clhs217*DeltaNormals[0](1,0); // CLHS217*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs219 =     clhs138*clhs212;
    const double clhs220 =     clhs125*clhs211*clhs4;
    const double clhs221 =     clhs13*clhs220;
    const double clhs222 =     clhs125*clhs2*clhs211;
    const double clhs223 =     clhs21*clhs222;
    const double clhs224 =     clhs144*clhs4*epsilon;
    const double clhs225 =     clhs2*clhs4*epsilon;
    const double clhs226 =     clhs13*clhs224 + clhs146*clhs21 + clhs148*clhs225;
    const double clhs227 =     clhs217*DeltaNormals[1](1,0); // CLHS217*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs228 =     clhs151*clhs212;
    const double clhs229 =     clhs220*clhs26;
    const double clhs230 =     clhs222*clhs29;
    const double clhs231 =     clhs146*clhs29 + clhs155*clhs225 + clhs224*clhs26;
    const double clhs232 =     clhs217*DeltaNormals[2](1,0); // CLHS217*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs233 =     clhs158*clhs212;
    const double clhs234 =     clhs220*clhs34;
    const double clhs235 =     clhs222*clhs37;
    const double clhs236 =     clhs146*clhs37 + clhs162*clhs225 + clhs224*clhs34;
    const double clhs237 =     clhs217*DeltaNormals[3](1,0); // CLHS217*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs238 =     clhs165*clhs212;
    const double clhs239 =     clhs220*clhs42;
    const double clhs240 =     clhs222*clhs45;
    const double clhs241 =     clhs146*clhs45 + clhs169*clhs225 + clhs224*clhs42;
    const double clhs242 =     clhs171*std::pow(clhs4, 2);
    const double clhs243 =     -clhs242;
    const double clhs244 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs245 =     clhs2*clhs244*clhs4;
    const double clhs246 =     -clhs126*clhs245;
    const double clhs247 =     -clhs130*clhs245;
    const double clhs248 =     -clhs132*clhs245;
    const double clhs249 =     -clhs134*clhs245;
    const double clhs250 =     clhs217*DeltaNormals[0](1,1); // CLHS217*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs251 =     clhs138*clhs245;
    const double clhs252 =     clhs125*clhs244*clhs4;
    const double clhs253 =     clhs13*clhs252;
    const double clhs254 =     clhs125*clhs2*clhs244;
    const double clhs255 =     clhs21*clhs254;
    const double clhs256 =     clhs188*clhs4*epsilon;
    const double clhs257 =     clhs13*clhs256 + clhs190*clhs21 + clhs191*clhs225;
    const double clhs258 =     clhs217*DeltaNormals[1](1,1); // CLHS217*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs259 =     clhs151*clhs245;
    const double clhs260 =     clhs252*clhs26;
    const double clhs261 =     clhs254*clhs29;
    const double clhs262 =     clhs190*clhs29 + clhs197*clhs225 + clhs256*clhs26;
    const double clhs263 =     clhs217*DeltaNormals[2](1,1); // CLHS217*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs264 =     clhs158*clhs245;
    const double clhs265 =     clhs252*clhs34;
    const double clhs266 =     clhs254*clhs37;
    const double clhs267 =     clhs190*clhs37 + clhs203*clhs225 + clhs256*clhs34;
    const double clhs268 =     clhs217*DeltaNormals[3](1,1); // CLHS217*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs269 =     clhs165*clhs245;
    const double clhs270 =     clhs252*clhs42;
    const double clhs271 =     clhs254*clhs45;
    const double clhs272 =     clhs190*clhs45 + clhs209*clhs225 + clhs256*clhs42;

    lhs(0,0)=clhs1*clhs9;
    lhs(0,1)=clhs10*clhs9;
    lhs(0,2)=clhs11*clhs9;
    lhs(0,3)=clhs12*clhs9;
    lhs(0,4)=clhs13*clhs14 + clhs13*clhs17 + clhs15*clhs16 + clhs15*clhs18 + clhs19*clhs22 + clhs19*clhs25;
    lhs(0,5)=clhs14*clhs26 + clhs16*clhs27 + clhs17*clhs26 + clhs18*clhs27 + clhs19*clhs30 + clhs19*clhs33;
    lhs(0,6)=clhs14*clhs34 + clhs16*clhs35 + clhs17*clhs34 + clhs18*clhs35 + clhs19*clhs38 + clhs19*clhs41;
    lhs(0,7)=clhs14*clhs42 + clhs16*clhs43 + clhs17*clhs42 + clhs18*clhs43 + clhs19*clhs46 + clhs19*clhs49;
    lhs(0,8)=clhs50;
    lhs(0,9)=0;
    lhs(0,10)=clhs51;
    lhs(0,11)=0;
    lhs(0,12)=clhs50;
    lhs(0,13)=0;
    lhs(0,14)=clhs51;
    lhs(0,15)=0;
    lhs(1,0)=clhs1*clhs56;
    lhs(1,1)=clhs10*clhs56;
    lhs(1,2)=clhs11*clhs56;
    lhs(1,3)=clhs12*clhs56;
    lhs(1,4)=clhs13*clhs57 + clhs13*clhs59 + clhs15*clhs58 + clhs15*clhs60 + clhs19*clhs61 + clhs19*clhs64;
    lhs(1,5)=clhs19*clhs65 + clhs19*clhs68 + clhs26*clhs57 + clhs26*clhs59 + clhs27*clhs58 + clhs27*clhs60;
    lhs(1,6)=clhs19*clhs69 + clhs19*clhs72 + clhs34*clhs57 + clhs34*clhs59 + clhs35*clhs58 + clhs35*clhs60;
    lhs(1,7)=clhs19*clhs73 + clhs19*clhs76 + clhs42*clhs57 + clhs42*clhs59 + clhs43*clhs58 + clhs43*clhs60;
    lhs(1,8)=0;
    lhs(1,9)=clhs50;
    lhs(1,10)=0;
    lhs(1,11)=clhs51;
    lhs(1,12)=0;
    lhs(1,13)=clhs50;
    lhs(1,14)=0;
    lhs(1,15)=clhs51;
    lhs(2,0)=clhs78*clhs9;
    lhs(2,1)=clhs79*clhs9;
    lhs(2,2)=clhs80*clhs9;
    lhs(2,3)=clhs81*clhs9;
    lhs(2,4)=clhs13*clhs82 + clhs13*clhs84 + clhs16*clhs83 + clhs18*clhs83 + clhs22*clhs85 + clhs25*clhs85;
    lhs(2,5)=clhs16*clhs86 + clhs18*clhs86 + clhs26*clhs82 + clhs26*clhs84 + clhs30*clhs85 + clhs33*clhs85;
    lhs(2,6)=clhs16*clhs87 + clhs18*clhs87 + clhs34*clhs82 + clhs34*clhs84 + clhs38*clhs85 + clhs41*clhs85;
    lhs(2,7)=clhs16*clhs88 + clhs18*clhs88 + clhs42*clhs82 + clhs42*clhs84 + clhs46*clhs85 + clhs49*clhs85;
    lhs(2,8)=clhs89;
    lhs(2,9)=0;
    lhs(2,10)=clhs90;
    lhs(2,11)=0;
    lhs(2,12)=clhs89;
    lhs(2,13)=0;
    lhs(2,14)=clhs90;
    lhs(2,15)=0;
    lhs(3,0)=clhs56*clhs78;
    lhs(3,1)=clhs56*clhs79;
    lhs(3,2)=clhs56*clhs80;
    lhs(3,3)=clhs56*clhs81;
    lhs(3,4)=clhs13*clhs91 + clhs13*clhs92 + clhs58*clhs83 + clhs60*clhs83 + clhs61*clhs85 + clhs64*clhs85;
    lhs(3,5)=clhs26*clhs91 + clhs26*clhs92 + clhs58*clhs86 + clhs60*clhs86 + clhs65*clhs85 + clhs68*clhs85;
    lhs(3,6)=clhs34*clhs91 + clhs34*clhs92 + clhs58*clhs87 + clhs60*clhs87 + clhs69*clhs85 + clhs72*clhs85;
    lhs(3,7)=clhs42*clhs91 + clhs42*clhs92 + clhs58*clhs88 + clhs60*clhs88 + clhs73*clhs85 + clhs76*clhs85;
    lhs(3,8)=0;
    lhs(3,9)=clhs89;
    lhs(3,10)=0;
    lhs(3,11)=clhs90;
    lhs(3,12)=0;
    lhs(3,13)=clhs89;
    lhs(3,14)=0;
    lhs(3,15)=clhs90;
    lhs(4,0)=-clhs9*clhs94;
    lhs(4,1)=-clhs9*clhs95;
    lhs(4,2)=-clhs9*clhs96;
    lhs(4,3)=-clhs9*clhs97;
    lhs(4,4)=-clhs100*clhs13 - clhs101*clhs22 - clhs101*clhs25 - clhs13*clhs98 - clhs16*clhs99 - clhs18*clhs99;
    lhs(4,5)=-clhs100*clhs26 - clhs101*clhs30 - clhs101*clhs33 - clhs102*clhs16 - clhs102*clhs18 - clhs26*clhs98;
    lhs(4,6)=-clhs100*clhs34 - clhs101*clhs38 - clhs101*clhs41 - clhs103*clhs16 - clhs103*clhs18 - clhs34*clhs98;
    lhs(4,7)=-clhs100*clhs42 - clhs101*clhs46 - clhs101*clhs49 - clhs104*clhs16 - clhs104*clhs18 - clhs42*clhs98;
    lhs(4,8)=clhs105;
    lhs(4,9)=0;
    lhs(4,10)=clhs106;
    lhs(4,11)=0;
    lhs(4,12)=clhs105;
    lhs(4,13)=0;
    lhs(4,14)=clhs106;
    lhs(4,15)=0;
    lhs(5,0)=-clhs56*clhs94;
    lhs(5,1)=-clhs56*clhs95;
    lhs(5,2)=-clhs56*clhs96;
    lhs(5,3)=-clhs56*clhs97;
    lhs(5,4)=-clhs101*clhs61 - clhs101*clhs64 - clhs107*clhs13 - clhs108*clhs13 - clhs58*clhs99 - clhs60*clhs99;
    lhs(5,5)=-clhs101*clhs65 - clhs101*clhs68 - clhs102*clhs58 - clhs102*clhs60 - clhs107*clhs26 - clhs108*clhs26;
    lhs(5,6)=-clhs101*clhs69 - clhs101*clhs72 - clhs103*clhs58 - clhs103*clhs60 - clhs107*clhs34 - clhs108*clhs34;
    lhs(5,7)=-clhs101*clhs73 - clhs101*clhs76 - clhs104*clhs58 - clhs104*clhs60 - clhs107*clhs42 - clhs108*clhs42;
    lhs(5,8)=0;
    lhs(5,9)=clhs105;
    lhs(5,10)=0;
    lhs(5,11)=clhs106;
    lhs(5,12)=0;
    lhs(5,13)=clhs105;
    lhs(5,14)=0;
    lhs(5,15)=clhs106;
    lhs(6,0)=-clhs110*clhs9;
    lhs(6,1)=-clhs111*clhs9;
    lhs(6,2)=-clhs112*clhs9;
    lhs(6,3)=-clhs113*clhs9;
    lhs(6,4)=-clhs114*clhs13 - clhs115*clhs16 - clhs115*clhs18 - clhs116*clhs13 - clhs117*clhs22 - clhs117*clhs25;
    lhs(6,5)=-clhs114*clhs26 - clhs116*clhs26 - clhs117*clhs30 - clhs117*clhs33 - clhs118*clhs16 - clhs118*clhs18;
    lhs(6,6)=-clhs114*clhs34 - clhs116*clhs34 - clhs117*clhs38 - clhs117*clhs41 - clhs119*clhs16 - clhs119*clhs18;
    lhs(6,7)=-clhs114*clhs42 - clhs116*clhs42 - clhs117*clhs46 - clhs117*clhs49 - clhs120*clhs16 - clhs120*clhs18;
    lhs(6,8)=clhs121;
    lhs(6,9)=0;
    lhs(6,10)=clhs122;
    lhs(6,11)=0;
    lhs(6,12)=clhs121;
    lhs(6,13)=0;
    lhs(6,14)=clhs122;
    lhs(6,15)=0;
    lhs(7,0)=-clhs110*clhs56;
    lhs(7,1)=-clhs111*clhs56;
    lhs(7,2)=-clhs112*clhs56;
    lhs(7,3)=-clhs113*clhs56;
    lhs(7,4)=-clhs115*clhs58 - clhs115*clhs60 - clhs117*clhs61 - clhs117*clhs64 - clhs123*clhs13 - clhs124*clhs13;
    lhs(7,5)=-clhs117*clhs65 - clhs117*clhs68 - clhs118*clhs58 - clhs118*clhs60 - clhs123*clhs26 - clhs124*clhs26;
    lhs(7,6)=-clhs117*clhs69 - clhs117*clhs72 - clhs119*clhs58 - clhs119*clhs60 - clhs123*clhs34 - clhs124*clhs34;
    lhs(7,7)=-clhs117*clhs73 - clhs117*clhs76 - clhs120*clhs58 - clhs120*clhs60 - clhs123*clhs42 - clhs124*clhs42;
    lhs(7,8)=0;
    lhs(7,9)=clhs121;
    lhs(7,10)=0;
    lhs(7,11)=clhs122;
    lhs(7,12)=0;
    lhs(7,13)=clhs121;
    lhs(7,14)=0;
    lhs(7,15)=clhs122;
    lhs(8,0)=clhs129;
    lhs(8,1)=clhs131;
    lhs(8,2)=clhs133;
    lhs(8,3)=clhs135;
    lhs(8,4)=-clhs137 - clhs139 - clhs141 - clhs143 + clhs149;
    lhs(8,5)=-clhs150 - clhs152 - clhs153 - clhs154 + clhs156;
    lhs(8,6)=-clhs157 - clhs159 - clhs160 - clhs161 + clhs163;
    lhs(8,7)=-clhs164 - clhs166 - clhs167 - clhs168 + clhs170;
    lhs(8,8)=clhs173;
    lhs(8,9)=0;
    lhs(8,10)=clhs175;
    lhs(8,11)=0;
    lhs(8,12)=clhs172;
    lhs(8,13)=0;
    lhs(8,14)=clhs174;
    lhs(8,15)=0;
    lhs(9,0)=clhs178;
    lhs(9,1)=clhs179;
    lhs(9,2)=clhs180;
    lhs(9,3)=clhs181;
    lhs(9,4)=-clhs182 - clhs183 - clhs185 - clhs187 + clhs192;
    lhs(9,5)=-clhs193 - clhs194 - clhs195 - clhs196 + clhs198;
    lhs(9,6)=-clhs199 - clhs200 - clhs201 - clhs202 + clhs204;
    lhs(9,7)=-clhs205 - clhs206 - clhs207 - clhs208 + clhs210;
    lhs(9,8)=0;
    lhs(9,9)=clhs173;
    lhs(9,10)=0;
    lhs(9,11)=clhs175;
    lhs(9,12)=0;
    lhs(9,13)=clhs172;
    lhs(9,14)=0;
    lhs(9,15)=clhs174;
    lhs(10,0)=clhs213;
    lhs(10,1)=clhs214;
    lhs(10,2)=clhs215;
    lhs(10,3)=clhs216;
    lhs(10,4)=-clhs218 - clhs219 - clhs221 - clhs223 + clhs226;
    lhs(10,5)=-clhs227 - clhs228 - clhs229 - clhs230 + clhs231;
    lhs(10,6)=-clhs232 - clhs233 - clhs234 - clhs235 + clhs236;
    lhs(10,7)=-clhs237 - clhs238 - clhs239 - clhs240 + clhs241;
    lhs(10,8)=clhs175;
    lhs(10,9)=0;
    lhs(10,10)=clhs243;
    lhs(10,11)=0;
    lhs(10,12)=clhs174;
    lhs(10,13)=0;
    lhs(10,14)=clhs242;
    lhs(10,15)=0;
    lhs(11,0)=clhs246;
    lhs(11,1)=clhs247;
    lhs(11,2)=clhs248;
    lhs(11,3)=clhs249;
    lhs(11,4)=-clhs250 - clhs251 - clhs253 - clhs255 + clhs257;
    lhs(11,5)=-clhs258 - clhs259 - clhs260 - clhs261 + clhs262;
    lhs(11,6)=-clhs263 - clhs264 - clhs265 - clhs266 + clhs267;
    lhs(11,7)=-clhs268 - clhs269 - clhs270 - clhs271 + clhs272;
    lhs(11,8)=0;
    lhs(11,9)=clhs175;
    lhs(11,10)=0;
    lhs(11,11)=clhs243;
    lhs(11,12)=0;
    lhs(11,13)=clhs174;
    lhs(11,14)=0;
    lhs(11,15)=clhs242;
    lhs(12,0)=clhs129;
    lhs(12,1)=clhs131;
    lhs(12,2)=clhs133;
    lhs(12,3)=clhs135;
    lhs(12,4)=-clhs137 - clhs139 - clhs141 - clhs143 - clhs149;
    lhs(12,5)=-clhs150 - clhs152 - clhs153 - clhs154 - clhs156;
    lhs(12,6)=-clhs157 - clhs159 - clhs160 - clhs161 - clhs163;
    lhs(12,7)=-clhs164 - clhs166 - clhs167 - clhs168 - clhs170;
    lhs(12,8)=clhs172;
    lhs(12,9)=0;
    lhs(12,10)=clhs174;
    lhs(12,11)=0;
    lhs(12,12)=clhs173;
    lhs(12,13)=0;
    lhs(12,14)=clhs175;
    lhs(12,15)=0;
    lhs(13,0)=clhs178;
    lhs(13,1)=clhs179;
    lhs(13,2)=clhs180;
    lhs(13,3)=clhs181;
    lhs(13,4)=-clhs182 - clhs183 - clhs185 - clhs187 - clhs192;
    lhs(13,5)=-clhs193 - clhs194 - clhs195 - clhs196 - clhs198;
    lhs(13,6)=-clhs199 - clhs200 - clhs201 - clhs202 - clhs204;
    lhs(13,7)=-clhs205 - clhs206 - clhs207 - clhs208 - clhs210;
    lhs(13,8)=0;
    lhs(13,9)=clhs172;
    lhs(13,10)=0;
    lhs(13,11)=clhs174;
    lhs(13,12)=0;
    lhs(13,13)=clhs173;
    lhs(13,14)=0;
    lhs(13,15)=clhs175;
    lhs(14,0)=clhs213;
    lhs(14,1)=clhs214;
    lhs(14,2)=clhs215;
    lhs(14,3)=clhs216;
    lhs(14,4)=-clhs218 - clhs219 - clhs221 - clhs223 - clhs226;
    lhs(14,5)=-clhs227 - clhs228 - clhs229 - clhs230 - clhs231;
    lhs(14,6)=-clhs232 - clhs233 - clhs234 - clhs235 - clhs236;
    lhs(14,7)=-clhs237 - clhs238 - clhs239 - clhs240 - clhs241;
    lhs(14,8)=clhs174;
    lhs(14,9)=0;
    lhs(14,10)=clhs242;
    lhs(14,11)=0;
    lhs(14,12)=clhs175;
    lhs(14,13)=0;
    lhs(14,14)=clhs243;
    lhs(14,15)=0;
    lhs(15,0)=clhs246;
    lhs(15,1)=clhs247;
    lhs(15,2)=clhs248;
    lhs(15,3)=clhs249;
    lhs(15,4)=-clhs250 - clhs251 - clhs253 - clhs255 - clhs257;
    lhs(15,5)=-clhs258 - clhs259 - clhs260 - clhs261 - clhs262;
    lhs(15,6)=-clhs263 - clhs264 - clhs265 - clhs266 - clhs267;
    lhs(15,7)=-clhs268 - clhs269 - clhs270 - clhs271 - clhs272;
    lhs(15,8)=0;
    lhs(15,9)=clhs174;
    lhs(15,10)=0;
    lhs(15,11)=clhs242;
    lhs(15,12)=0;
    lhs(15,13)=clhs175;
    lhs(15,14)=0;
    lhs(15,15)=clhs243;

    
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
    
    const bounded_matrix<double, 2, 2> v1 = rContactData.v1;
    const bounded_matrix<double, 2, 2> v2 = rContactData.v2;
    
    const std::vector<double> DeltaJs = rContactData.DeltaJ_s;
    const std::vector<array_1d<double,2>> DeltaPhi = rContactData.DeltaPhi;
    const std::vector<array_1d<double,2>> DeltaN1  = rContactData.DeltaN1;
    const std::vector<array_1d<double,2>> DeltaN2  = rContactData.DeltaN2;
    const std::vector<bounded_matrix<double, 2, 2>> DeltaNormals    = rContactData.Delta_Normal_s;
    const std::vector<bounded_matrix<double, 2, 2>> Deltatangentxis = rContactData.Delta_Tangent_xi_s;
    
    const double clhs0 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs1 =     1.0/Dt;
    const double clhs2 =     Phi[0]; // PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs3 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs4 =     N1[0]; // N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs5 =     DeltaN1[4][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs6 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs7 =     N1[1]; // N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs8 =     DeltaN1[4][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs9 =     Dt*v1(0,0);
    const double clhs10 =     Dt*v1(1,0);
    const double clhs11 =     Dt*v2(0,0);
    const double clhs12 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs13 =     Dt*v2(1,0);
    const double clhs14 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs15 =     clhs10*clhs7 - clhs11*clhs12 - clhs13*clhs14 + clhs4*clhs9;
    const double clhs16 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs17 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs18 =     Dt*v1(0,1);
    const double clhs19 =     Dt*v1(1,1);
    const double clhs20 =     Dt*v2(0,1);
    const double clhs21 =     Dt*v2(1,1);
    const double clhs22 =     -clhs12*clhs20 - clhs14*clhs21 + clhs18*clhs4 + clhs19*clhs7;
    const double clhs23 =     clhs16*clhs4 + clhs17*clhs7;
    const double clhs24 =     DeltaN2[4][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs25 =     DeltaN2[4][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs26 =     clhs0*clhs4 + clhs6*clhs7;
    const double clhs27 =     clhs15*(clhs0*clhs5 + clhs6*clhs8) + clhs22*(clhs16*clhs5 + clhs17*clhs8) + clhs23*(clhs18*clhs5 + clhs19*clhs8 - clhs20*clhs24 - clhs21*clhs25) - clhs26*(-clhs10*clhs8 + clhs11*clhs24 + clhs12 + clhs13*clhs25 - clhs5*clhs9);
    const double clhs28 =     clhs1*clhs2*clhs27*clhs3;
    const double clhs29 =     clhs0*clhs28;
    const double clhs30 =     DeltaN1[5][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs31 =     DeltaN1[5][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs32 =     DeltaN2[5][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs33 =     DeltaN2[5][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs34 =     clhs15*(clhs0*clhs30 + clhs31*clhs6) + clhs22*(clhs16*clhs30 + clhs17*clhs31) - clhs23*(clhs12 - clhs18*clhs30 - clhs19*clhs31 + clhs20*clhs32 + clhs21*clhs33) + clhs26*(clhs10*clhs31 - clhs11*clhs32 - clhs13*clhs33 + clhs30*clhs9);
    const double clhs35 =     clhs1*clhs2*clhs3*clhs34;
    const double clhs36 =     clhs0*clhs35;
    const double clhs37 =     DeltaN1[6][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs38 =     DeltaN1[6][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs39 =     DeltaN2[6][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs40 =     DeltaN2[6][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs41 =     clhs15*(clhs0*clhs37 + clhs38*clhs6) + clhs22*(clhs16*clhs37 + clhs17*clhs38) + clhs23*(clhs18*clhs37 + clhs19*clhs38 - clhs20*clhs39 - clhs21*clhs40) - clhs26*(-clhs10*clhs38 + clhs11*clhs39 + clhs13*clhs40 + clhs14 - clhs37*clhs9);
    const double clhs42 =     clhs1*clhs2*clhs3*clhs41;
    const double clhs43 =     clhs0*clhs42;
    const double clhs44 =     DeltaN1[7][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs45 =     DeltaN1[7][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs46 =     DeltaN2[7][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs47 =     DeltaN2[7][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs48 =     clhs15*(clhs0*clhs44 + clhs45*clhs6) + clhs22*(clhs16*clhs44 + clhs17*clhs45) - clhs23*(clhs14 - clhs18*clhs44 - clhs19*clhs45 + clhs20*clhs46 + clhs21*clhs47) + clhs26*(clhs10*clhs45 - clhs11*clhs46 - clhs13*clhs47 + clhs44*clhs9);
    const double clhs49 =     clhs1*clhs2*clhs3*clhs48;
    const double clhs50 =     clhs0*clhs49;
    const double clhs51 =     Deltatangentxis[0](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs52 =     clhs15*clhs26 + clhs22*clhs23;
    const double clhs53 =     clhs2*clhs3*clhs52;
    const double clhs54 =     DeltaJs[0]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs55 =     clhs0*clhs2*clhs52;
    const double clhs56 =     DeltaPhi[0][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs57 =     clhs0*clhs3*clhs52;
    const double clhs58 =     DeltaN1[0][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs59 =     DeltaN1[0][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs60 =     DeltaN2[0][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs61 =     DeltaN2[0][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs62 =     Deltatangentxis[0](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs63 =     Deltatangentxis[0](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs64 =     Deltatangentxis[0](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs65 =     clhs15*(clhs0*clhs58 + clhs4*clhs51 + clhs59*clhs6 + clhs62*clhs7) + clhs22*(clhs16*clhs58 + clhs17*clhs59 + clhs4*clhs63 + clhs64*clhs7) + clhs23*(clhs18*clhs58 + clhs19*clhs59 - clhs20*clhs60 - clhs21*clhs61) + clhs26*(clhs10*clhs59 - clhs11*clhs60 - clhs13*clhs61 + clhs4 + clhs58*clhs9);
    const double clhs66 =     clhs2*clhs3*clhs65;
    const double clhs67 =     clhs1*(clhs0*clhs66 + clhs51*clhs53 + clhs54*clhs55 + clhs56*clhs57);
    const double clhs68 =     Deltatangentxis[1](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs69 =     DeltaJs[1]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs70 =     DeltaPhi[1][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs71 =     DeltaN1[1][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs72 =     DeltaN1[1][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs73 =     DeltaN2[1][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs74 =     DeltaN2[1][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs75 =     Deltatangentxis[1](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs76 =     Deltatangentxis[1](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs77 =     Deltatangentxis[1](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs78 =     clhs15*(clhs0*clhs71 + clhs4*clhs68 + clhs6*clhs72 + clhs7*clhs75) + clhs22*(clhs16*clhs71 + clhs17*clhs72 + clhs4*clhs76 + clhs7*clhs77) + clhs23*(clhs18*clhs71 + clhs19*clhs72 - clhs20*clhs73 - clhs21*clhs74 + clhs4) + clhs26*(clhs10*clhs72 - clhs11*clhs73 - clhs13*clhs74 + clhs71*clhs9);
    const double clhs79 =     clhs2*clhs3*clhs78;
    const double clhs80 =     clhs1*(clhs0*clhs79 + clhs53*clhs68 + clhs55*clhs69 + clhs57*clhs70);
    const double clhs81 =     Deltatangentxis[2](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs82 =     DeltaJs[2]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs83 =     DeltaPhi[2][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs84 =     DeltaN1[2][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs85 =     DeltaN1[2][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs86 =     DeltaN2[2][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs87 =     DeltaN2[2][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs88 =     Deltatangentxis[2](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs89 =     Deltatangentxis[2](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs90 =     Deltatangentxis[2](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs91 =     clhs15*(clhs0*clhs84 + clhs4*clhs81 + clhs6*clhs85 + clhs7*clhs88) + clhs22*(clhs16*clhs84 + clhs17*clhs85 + clhs4*clhs89 + clhs7*clhs90) + clhs23*(clhs18*clhs84 + clhs19*clhs85 - clhs20*clhs86 - clhs21*clhs87) + clhs26*(clhs10*clhs85 - clhs11*clhs86 - clhs13*clhs87 + clhs7 + clhs84*clhs9);
    const double clhs92 =     clhs2*clhs3*clhs91;
    const double clhs93 =     clhs1*(clhs0*clhs92 + clhs53*clhs81 + clhs55*clhs82 + clhs57*clhs83);
    const double clhs94 =     Deltatangentxis[3](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs95 =     DeltaJs[3]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs96 =     DeltaPhi[3][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs97 =     DeltaN1[3][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs98 =     DeltaN1[3][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs99 =     DeltaN2[3][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs100 =     DeltaN2[3][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs101 =     Deltatangentxis[3](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs102 =     Deltatangentxis[3](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs103 =     Deltatangentxis[3](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs104 =     clhs15*(clhs0*clhs97 + clhs101*clhs7 + clhs4*clhs94 + clhs6*clhs98) + clhs22*(clhs102*clhs4 + clhs103*clhs7 + clhs16*clhs97 + clhs17*clhs98) + clhs23*(-clhs100*clhs21 + clhs18*clhs97 + clhs19*clhs98 - clhs20*clhs99 + clhs7) + clhs26*(clhs10*clhs98 - clhs100*clhs13 - clhs11*clhs99 + clhs9*clhs97);
    const double clhs105 =     clhs104*clhs2*clhs3;
    const double clhs106 =     clhs1*(clhs0*clhs105 + clhs53*clhs94 + clhs55*clhs95 + clhs57*clhs96);
    const double clhs107 =     clhs16*clhs28;
    const double clhs108 =     clhs16*clhs35;
    const double clhs109 =     clhs16*clhs42;
    const double clhs110 =     clhs16*clhs49;
    const double clhs111 =     clhs16*clhs2*clhs52;
    const double clhs112 =     clhs16*clhs3*clhs52;
    const double clhs113 =     clhs1*(clhs111*clhs54 + clhs112*clhs56 + clhs16*clhs66 + clhs53*clhs63);
    const double clhs114 =     clhs1*(clhs111*clhs69 + clhs112*clhs70 + clhs16*clhs79 + clhs53*clhs76);
    const double clhs115 =     clhs1*(clhs111*clhs82 + clhs112*clhs83 + clhs16*clhs92 + clhs53*clhs89);
    const double clhs116 =     clhs1*(clhs102*clhs53 + clhs105*clhs16 + clhs111*clhs95 + clhs112*clhs96);
    const double clhs117 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs118 =     clhs1*clhs117*clhs27*clhs3;
    const double clhs119 =     clhs118*clhs6;
    const double clhs120 =     clhs1*clhs117*clhs3*clhs34;
    const double clhs121 =     clhs120*clhs6;
    const double clhs122 =     clhs1*clhs117*clhs3*clhs41;
    const double clhs123 =     clhs122*clhs6;
    const double clhs124 =     clhs1*clhs117*clhs3*clhs48;
    const double clhs125 =     clhs124*clhs6;
    const double clhs126 =     clhs117*clhs3*clhs52;
    const double clhs127 =     clhs117*clhs52*clhs6;
    const double clhs128 =     DeltaPhi[0][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs129 =     clhs3*clhs52*clhs6;
    const double clhs130 =     clhs117*clhs3*clhs65;
    const double clhs131 =     clhs1*(clhs126*clhs62 + clhs127*clhs54 + clhs128*clhs129 + clhs130*clhs6);
    const double clhs132 =     DeltaPhi[1][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs133 =     clhs117*clhs3*clhs78;
    const double clhs134 =     clhs1*(clhs126*clhs75 + clhs127*clhs69 + clhs129*clhs132 + clhs133*clhs6);
    const double clhs135 =     DeltaPhi[2][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs136 =     clhs117*clhs3*clhs91;
    const double clhs137 =     clhs1*(clhs126*clhs88 + clhs127*clhs82 + clhs129*clhs135 + clhs136*clhs6);
    const double clhs138 =     DeltaPhi[3][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs139 =     clhs104*clhs117*clhs3;
    const double clhs140 =     clhs1*(clhs101*clhs126 + clhs127*clhs95 + clhs129*clhs138 + clhs139*clhs6);
    const double clhs141 =     clhs118*clhs17;
    const double clhs142 =     clhs120*clhs17;
    const double clhs143 =     clhs122*clhs17;
    const double clhs144 =     clhs124*clhs17;
    const double clhs145 =     clhs117*clhs17*clhs52;
    const double clhs146 =     clhs17*clhs3*clhs52;
    const double clhs147 =     clhs1*(clhs126*clhs64 + clhs128*clhs146 + clhs130*clhs17 + clhs145*clhs54);
    const double clhs148 =     clhs1*(clhs126*clhs77 + clhs132*clhs146 + clhs133*clhs17 + clhs145*clhs69);
    const double clhs149 =     clhs1*(clhs126*clhs90 + clhs135*clhs146 + clhs136*clhs17 + clhs145*clhs82);
    const double clhs150 =     clhs1*(clhs103*clhs126 + clhs138*clhs146 + clhs139*clhs17 + clhs145*clhs95);

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
    lhs(8,1)=clhs36;
    lhs(8,2)=clhs43;
    lhs(8,3)=clhs50;
    lhs(8,4)=clhs67;
    lhs(8,5)=clhs80;
    lhs(8,6)=clhs93;
    lhs(8,7)=clhs106;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(8,12)=0;
    lhs(8,13)=0;
    lhs(8,14)=0;
    lhs(8,15)=0;
    lhs(9,0)=clhs107;
    lhs(9,1)=clhs108;
    lhs(9,2)=clhs109;
    lhs(9,3)=clhs110;
    lhs(9,4)=clhs113;
    lhs(9,5)=clhs114;
    lhs(9,6)=clhs115;
    lhs(9,7)=clhs116;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(9,12)=0;
    lhs(9,13)=0;
    lhs(9,14)=0;
    lhs(9,15)=0;
    lhs(10,0)=clhs119;
    lhs(10,1)=clhs121;
    lhs(10,2)=clhs123;
    lhs(10,3)=clhs125;
    lhs(10,4)=clhs131;
    lhs(10,5)=clhs134;
    lhs(10,6)=clhs137;
    lhs(10,7)=clhs140;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(10,12)=0;
    lhs(10,13)=0;
    lhs(10,14)=0;
    lhs(10,15)=0;
    lhs(11,0)=clhs141;
    lhs(11,1)=clhs142;
    lhs(11,2)=clhs143;
    lhs(11,3)=clhs144;
    lhs(11,4)=clhs147;
    lhs(11,5)=clhs148;
    lhs(11,6)=clhs149;
    lhs(11,7)=clhs150;
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;
    lhs(11,12)=0;
    lhs(11,13)=0;
    lhs(11,14)=0;
    lhs(11,15)=0;
    lhs(12,0)=clhs29;
    lhs(12,1)=clhs36;
    lhs(12,2)=clhs43;
    lhs(12,3)=clhs50;
    lhs(12,4)=clhs67;
    lhs(12,5)=clhs80;
    lhs(12,6)=clhs93;
    lhs(12,7)=clhs106;
    lhs(12,8)=0;
    lhs(12,9)=0;
    lhs(12,10)=0;
    lhs(12,11)=0;
    lhs(12,12)=0;
    lhs(12,13)=0;
    lhs(12,14)=0;
    lhs(12,15)=0;
    lhs(13,0)=clhs107;
    lhs(13,1)=clhs108;
    lhs(13,2)=clhs109;
    lhs(13,3)=clhs110;
    lhs(13,4)=clhs113;
    lhs(13,5)=clhs114;
    lhs(13,6)=clhs115;
    lhs(13,7)=clhs116;
    lhs(13,8)=0;
    lhs(13,9)=0;
    lhs(13,10)=0;
    lhs(13,11)=0;
    lhs(13,12)=0;
    lhs(13,13)=0;
    lhs(13,14)=0;
    lhs(13,15)=0;
    lhs(14,0)=clhs119;
    lhs(14,1)=clhs121;
    lhs(14,2)=clhs123;
    lhs(14,3)=clhs125;
    lhs(14,4)=clhs131;
    lhs(14,5)=clhs134;
    lhs(14,6)=clhs137;
    lhs(14,7)=clhs140;
    lhs(14,8)=0;
    lhs(14,9)=0;
    lhs(14,10)=0;
    lhs(14,11)=0;
    lhs(14,12)=0;
    lhs(14,13)=0;
    lhs(14,14)=0;
    lhs(14,15)=0;
    lhs(15,0)=clhs141;
    lhs(15,1)=clhs142;
    lhs(15,2)=clhs143;
    lhs(15,3)=clhs144;
    lhs(15,4)=clhs147;
    lhs(15,5)=clhs148;
    lhs(15,6)=clhs149;
    lhs(15,7)=clhs150;
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
    const bounded_matrix<double, 2, 2> dlm             = rContactData.DoubleLagrangeMultipliers;
    const double Dt              = rContactData.Dt;
//     const double epsilon  = rContactData.epsilon;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
    
//     const double sign_tangpress = boost::math::sign(augmented_tangent_lm);
    
    const bounded_matrix<double, 2, 2> v1 = rContactData.v1;
    const bounded_matrix<double, 2, 2> v2 = rContactData.v2;
    
    const std::vector<double> DeltaJs = rContactData.DeltaJ_s;
    const std::vector<array_1d<double,2>> DeltaPhi = rContactData.DeltaPhi;
    const std::vector<array_1d<double,2>> DeltaN1  = rContactData.DeltaN1;
    const std::vector<array_1d<double,2>> DeltaN2  = rContactData.DeltaN2;
    const std::vector<bounded_matrix<double, 2, 2>> DeltaNormals    = rContactData.Delta_Normal_s;
    const std::vector<bounded_matrix<double, 2, 2>> Deltatangentxis = rContactData.Delta_Tangent_xi_s;
    
    const double clhs0 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs1 =     1.0/Dt;
    const double clhs2 =     Phi[0]; // PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs3 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs4 =     N1[0]; // N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs5 =     DeltaN1[4][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs6 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs7 =     N1[1]; // N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs8 =     DeltaN1[4][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs9 =     Dt*v1(0,0);
    const double clhs10 =     Dt*v1(1,0);
    const double clhs11 =     Dt*v2(0,0);
    const double clhs12 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs13 =     Dt*v2(1,0);
    const double clhs14 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs15 =     clhs10*clhs7 - clhs11*clhs12 - clhs13*clhs14 + clhs4*clhs9;
    const double clhs16 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs17 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs18 =     Dt*v1(0,1);
    const double clhs19 =     Dt*v1(1,1);
    const double clhs20 =     Dt*v2(0,1);
    const double clhs21 =     Dt*v2(1,1);
    const double clhs22 =     -clhs12*clhs20 - clhs14*clhs21 + clhs18*clhs4 + clhs19*clhs7;
    const double clhs23 =     clhs16*clhs4 + clhs17*clhs7;
    const double clhs24 =     DeltaN2[4][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs25 =     DeltaN2[4][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs26 =     clhs0*clhs4 + clhs6*clhs7;
    const double clhs27 =     clhs15*(clhs0*clhs5 + clhs6*clhs8) + clhs22*(clhs16*clhs5 + clhs17*clhs8) + clhs23*(clhs18*clhs5 + clhs19*clhs8 - clhs20*clhs24 - clhs21*clhs25) - clhs26*(-clhs10*clhs8 + clhs11*clhs24 + clhs12 + clhs13*clhs25 - clhs5*clhs9);
    const double clhs28 =     clhs1*clhs2*clhs27*clhs3;
    const double clhs29 =     clhs0*clhs28;
    const double clhs30 =     DeltaN1[5][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs31 =     DeltaN1[5][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs32 =     DeltaN2[5][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs33 =     DeltaN2[5][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs34 =     clhs15*(clhs0*clhs30 + clhs31*clhs6) + clhs22*(clhs16*clhs30 + clhs17*clhs31) - clhs23*(clhs12 - clhs18*clhs30 - clhs19*clhs31 + clhs20*clhs32 + clhs21*clhs33) + clhs26*(clhs10*clhs31 - clhs11*clhs32 - clhs13*clhs33 + clhs30*clhs9);
    const double clhs35 =     clhs1*clhs2*clhs3*clhs34;
    const double clhs36 =     clhs0*clhs35;
    const double clhs37 =     DeltaN1[6][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs38 =     DeltaN1[6][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs39 =     DeltaN2[6][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs40 =     DeltaN2[6][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs41 =     clhs15*(clhs0*clhs37 + clhs38*clhs6) + clhs22*(clhs16*clhs37 + clhs17*clhs38) + clhs23*(clhs18*clhs37 + clhs19*clhs38 - clhs20*clhs39 - clhs21*clhs40) - clhs26*(-clhs10*clhs38 + clhs11*clhs39 + clhs13*clhs40 + clhs14 - clhs37*clhs9);
    const double clhs42 =     clhs1*clhs2*clhs3*clhs41;
    const double clhs43 =     clhs0*clhs42;
    const double clhs44 =     DeltaN1[7][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs45 =     DeltaN1[7][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs46 =     DeltaN2[7][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs47 =     DeltaN2[7][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs48 =     clhs15*(clhs0*clhs44 + clhs45*clhs6) + clhs22*(clhs16*clhs44 + clhs17*clhs45) - clhs23*(clhs14 - clhs18*clhs44 - clhs19*clhs45 + clhs20*clhs46 + clhs21*clhs47) + clhs26*(clhs10*clhs45 - clhs11*clhs46 - clhs13*clhs47 + clhs44*clhs9);
    const double clhs49 =     clhs1*clhs2*clhs3*clhs48;
    const double clhs50 =     clhs0*clhs49;
    const double clhs51 =     Deltatangentxis[0](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs52 =     clhs15*clhs26 + clhs22*clhs23;
    const double clhs53 =     clhs2*clhs3*clhs52;
    const double clhs54 =     DeltaJs[0]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs55 =     clhs0*clhs2*clhs52;
    const double clhs56 =     DeltaPhi[0][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs57 =     clhs0*clhs3*clhs52;
    const double clhs58 =     DeltaN1[0][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs59 =     DeltaN1[0][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs60 =     DeltaN2[0][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs61 =     DeltaN2[0][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs62 =     Deltatangentxis[0](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs63 =     Deltatangentxis[0](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs64 =     Deltatangentxis[0](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs65 =     clhs15*(clhs0*clhs58 + clhs4*clhs51 + clhs59*clhs6 + clhs62*clhs7) + clhs22*(clhs16*clhs58 + clhs17*clhs59 + clhs4*clhs63 + clhs64*clhs7) + clhs23*(clhs18*clhs58 + clhs19*clhs59 - clhs20*clhs60 - clhs21*clhs61) + clhs26*(clhs10*clhs59 - clhs11*clhs60 - clhs13*clhs61 + clhs4 + clhs58*clhs9);
    const double clhs66 =     clhs2*clhs3*clhs65;
    const double clhs67 =     clhs1*(clhs0*clhs66 + clhs51*clhs53 + clhs54*clhs55 + clhs56*clhs57);
    const double clhs68 =     Deltatangentxis[1](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs69 =     DeltaJs[1]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs70 =     DeltaPhi[1][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs71 =     DeltaN1[1][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs72 =     DeltaN1[1][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs73 =     DeltaN2[1][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs74 =     DeltaN2[1][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs75 =     Deltatangentxis[1](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs76 =     Deltatangentxis[1](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs77 =     Deltatangentxis[1](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs78 =     clhs15*(clhs0*clhs71 + clhs4*clhs68 + clhs6*clhs72 + clhs7*clhs75) + clhs22*(clhs16*clhs71 + clhs17*clhs72 + clhs4*clhs76 + clhs7*clhs77) + clhs23*(clhs18*clhs71 + clhs19*clhs72 - clhs20*clhs73 - clhs21*clhs74 + clhs4) + clhs26*(clhs10*clhs72 - clhs11*clhs73 - clhs13*clhs74 + clhs71*clhs9);
    const double clhs79 =     clhs2*clhs3*clhs78;
    const double clhs80 =     clhs1*(clhs0*clhs79 + clhs53*clhs68 + clhs55*clhs69 + clhs57*clhs70);
    const double clhs81 =     Deltatangentxis[2](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs82 =     DeltaJs[2]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs83 =     DeltaPhi[2][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs84 =     DeltaN1[2][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs85 =     DeltaN1[2][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs86 =     DeltaN2[2][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs87 =     DeltaN2[2][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs88 =     Deltatangentxis[2](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs89 =     Deltatangentxis[2](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs90 =     Deltatangentxis[2](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs91 =     clhs15*(clhs0*clhs84 + clhs4*clhs81 + clhs6*clhs85 + clhs7*clhs88) + clhs22*(clhs16*clhs84 + clhs17*clhs85 + clhs4*clhs89 + clhs7*clhs90) + clhs23*(clhs18*clhs84 + clhs19*clhs85 - clhs20*clhs86 - clhs21*clhs87) + clhs26*(clhs10*clhs85 - clhs11*clhs86 - clhs13*clhs87 + clhs7 + clhs84*clhs9);
    const double clhs92 =     clhs2*clhs3*clhs91;
    const double clhs93 =     clhs1*(clhs0*clhs92 + clhs53*clhs81 + clhs55*clhs82 + clhs57*clhs83);
    const double clhs94 =     Deltatangentxis[3](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs95 =     DeltaJs[3]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs96 =     DeltaPhi[3][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs97 =     DeltaN1[3][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs98 =     DeltaN1[3][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs99 =     DeltaN2[3][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs100 =     DeltaN2[3][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs101 =     Deltatangentxis[3](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs102 =     Deltatangentxis[3](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs103 =     Deltatangentxis[3](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs104 =     clhs15*(clhs0*clhs97 + clhs101*clhs7 + clhs4*clhs94 + clhs6*clhs98) + clhs22*(clhs102*clhs4 + clhs103*clhs7 + clhs16*clhs97 + clhs17*clhs98) + clhs23*(-clhs100*clhs21 + clhs18*clhs97 + clhs19*clhs98 - clhs20*clhs99 + clhs7) + clhs26*(clhs10*clhs98 - clhs100*clhs13 - clhs11*clhs99 + clhs9*clhs97);
    const double clhs105 =     clhs104*clhs2*clhs3;
    const double clhs106 =     clhs1*(clhs0*clhs105 + clhs53*clhs94 + clhs55*clhs95 + clhs57*clhs96);
    const double clhs107 =     clhs16*clhs28;
    const double clhs108 =     clhs16*clhs35;
    const double clhs109 =     clhs16*clhs42;
    const double clhs110 =     clhs16*clhs49;
    const double clhs111 =     clhs16*clhs2*clhs52;
    const double clhs112 =     clhs16*clhs3*clhs52;
    const double clhs113 =     clhs1*(clhs111*clhs54 + clhs112*clhs56 + clhs16*clhs66 + clhs53*clhs63);
    const double clhs114 =     clhs1*(clhs111*clhs69 + clhs112*clhs70 + clhs16*clhs79 + clhs53*clhs76);
    const double clhs115 =     clhs1*(clhs111*clhs82 + clhs112*clhs83 + clhs16*clhs92 + clhs53*clhs89);
    const double clhs116 =     clhs1*(clhs102*clhs53 + clhs105*clhs16 + clhs111*clhs95 + clhs112*clhs96);
    const double clhs117 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs118 =     clhs1*clhs117*clhs27*clhs3;
    const double clhs119 =     clhs118*clhs6;
    const double clhs120 =     clhs1*clhs117*clhs3*clhs34;
    const double clhs121 =     clhs120*clhs6;
    const double clhs122 =     clhs1*clhs117*clhs3*clhs41;
    const double clhs123 =     clhs122*clhs6;
    const double clhs124 =     clhs1*clhs117*clhs3*clhs48;
    const double clhs125 =     clhs124*clhs6;
    const double clhs126 =     clhs117*clhs3*clhs52;
    const double clhs127 =     clhs117*clhs52*clhs6;
    const double clhs128 =     DeltaPhi[0][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs129 =     clhs3*clhs52*clhs6;
    const double clhs130 =     clhs117*clhs3*clhs65;
    const double clhs131 =     clhs1*(clhs126*clhs62 + clhs127*clhs54 + clhs128*clhs129 + clhs130*clhs6);
    const double clhs132 =     DeltaPhi[1][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs133 =     clhs117*clhs3*clhs78;
    const double clhs134 =     clhs1*(clhs126*clhs75 + clhs127*clhs69 + clhs129*clhs132 + clhs133*clhs6);
    const double clhs135 =     DeltaPhi[2][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs136 =     clhs117*clhs3*clhs91;
    const double clhs137 =     clhs1*(clhs126*clhs88 + clhs127*clhs82 + clhs129*clhs135 + clhs136*clhs6);
    const double clhs138 =     DeltaPhi[3][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs139 =     clhs104*clhs117*clhs3;
    const double clhs140 =     clhs1*(clhs101*clhs126 + clhs127*clhs95 + clhs129*clhs138 + clhs139*clhs6);
    const double clhs141 =     clhs118*clhs17;
    const double clhs142 =     clhs120*clhs17;
    const double clhs143 =     clhs122*clhs17;
    const double clhs144 =     clhs124*clhs17;
    const double clhs145 =     clhs117*clhs17*clhs52;
    const double clhs146 =     clhs17*clhs3*clhs52;
    const double clhs147 =     clhs1*(clhs126*clhs64 + clhs128*clhs146 + clhs130*clhs17 + clhs145*clhs54);
    const double clhs148 =     clhs1*(clhs126*clhs77 + clhs132*clhs146 + clhs133*clhs17 + clhs145*clhs69);
    const double clhs149 =     clhs1*(clhs126*clhs90 + clhs135*clhs146 + clhs136*clhs17 + clhs145*clhs82);
    const double clhs150 =     clhs1*(clhs103*clhs126 + clhs138*clhs146 + clhs139*clhs17 + clhs145*clhs95);

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
    lhs(8,1)=clhs36;
    lhs(8,2)=clhs43;
    lhs(8,3)=clhs50;
    lhs(8,4)=clhs67;
    lhs(8,5)=clhs80;
    lhs(8,6)=clhs93;
    lhs(8,7)=clhs106;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(8,12)=0;
    lhs(8,13)=0;
    lhs(8,14)=0;
    lhs(8,15)=0;
    lhs(9,0)=clhs107;
    lhs(9,1)=clhs108;
    lhs(9,2)=clhs109;
    lhs(9,3)=clhs110;
    lhs(9,4)=clhs113;
    lhs(9,5)=clhs114;
    lhs(9,6)=clhs115;
    lhs(9,7)=clhs116;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(9,12)=0;
    lhs(9,13)=0;
    lhs(9,14)=0;
    lhs(9,15)=0;
    lhs(10,0)=clhs119;
    lhs(10,1)=clhs121;
    lhs(10,2)=clhs123;
    lhs(10,3)=clhs125;
    lhs(10,4)=clhs131;
    lhs(10,5)=clhs134;
    lhs(10,6)=clhs137;
    lhs(10,7)=clhs140;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(10,12)=0;
    lhs(10,13)=0;
    lhs(10,14)=0;
    lhs(10,15)=0;
    lhs(11,0)=clhs141;
    lhs(11,1)=clhs142;
    lhs(11,2)=clhs143;
    lhs(11,3)=clhs144;
    lhs(11,4)=clhs147;
    lhs(11,5)=clhs148;
    lhs(11,6)=clhs149;
    lhs(11,7)=clhs150;
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;
    lhs(11,12)=0;
    lhs(11,13)=0;
    lhs(11,14)=0;
    lhs(11,15)=0;
    lhs(12,0)=clhs29;
    lhs(12,1)=clhs36;
    lhs(12,2)=clhs43;
    lhs(12,3)=clhs50;
    lhs(12,4)=clhs67;
    lhs(12,5)=clhs80;
    lhs(12,6)=clhs93;
    lhs(12,7)=clhs106;
    lhs(12,8)=0;
    lhs(12,9)=0;
    lhs(12,10)=0;
    lhs(12,11)=0;
    lhs(12,12)=0;
    lhs(12,13)=0;
    lhs(12,14)=0;
    lhs(12,15)=0;
    lhs(13,0)=clhs107;
    lhs(13,1)=clhs108;
    lhs(13,2)=clhs109;
    lhs(13,3)=clhs110;
    lhs(13,4)=clhs113;
    lhs(13,5)=clhs114;
    lhs(13,6)=clhs115;
    lhs(13,7)=clhs116;
    lhs(13,8)=0;
    lhs(13,9)=0;
    lhs(13,10)=0;
    lhs(13,11)=0;
    lhs(13,12)=0;
    lhs(13,13)=0;
    lhs(13,14)=0;
    lhs(13,15)=0;
    lhs(14,0)=clhs119;
    lhs(14,1)=clhs121;
    lhs(14,2)=clhs123;
    lhs(14,3)=clhs125;
    lhs(14,4)=clhs131;
    lhs(14,5)=clhs134;
    lhs(14,6)=clhs137;
    lhs(14,7)=clhs140;
    lhs(14,8)=0;
    lhs(14,9)=0;
    lhs(14,10)=0;
    lhs(14,11)=0;
    lhs(14,12)=0;
    lhs(14,13)=0;
    lhs(14,14)=0;
    lhs(14,15)=0;
    lhs(15,0)=clhs141;
    lhs(15,1)=clhs142;
    lhs(15,2)=clhs143;
    lhs(15,3)=clhs144;
    lhs(15,4)=clhs147;
    lhs(15,5)=clhs148;
    lhs(15,6)=clhs149;
    lhs(15,7)=clhs150;
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
//     const bounded_matrix<double, 2, 2> dlm             = rContactData.DoubleLagrangeMultipliers;
//     const double Dt              = rContactData.Dt;
//     const double epsilon  = rContactData.epsilon;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
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
    const bounded_matrix<double, 2, 2> dlm            = rContactData.DoubleLagrangeMultipliers;
//     const double Dt             = rContactData.Dt;
    const double epsilon = rContactData.epsilon;
    
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
    const double crhs19 =     -crhs18;
    const double crhs20 =     crhs16*normalslave(0,1); // CRHS16*NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs21 =     epsilon*(-crhs10 + crhs8 - crhs9);
    const double crhs22 =     -crhs21;
    const double crhs23 =     crhs1*crhs3;
    const double crhs24 =     crhs16*normalslave(1,0); // CRHS16*NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs25 =     crhs16*normalslave(1,1); // CRHS16*NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))

    rhs[0]=-crhs0*crhs7;
    rhs[1]=-crhs0*crhs11;
    rhs[2]=-crhs12*crhs7;
    rhs[3]=-crhs11*crhs12;
    rhs[4]=crhs13*crhs7;
    rhs[5]=crhs11*crhs13;
    rhs[6]=crhs14*crhs7;
    rhs[7]=crhs11*crhs14;
    rhs[8]=crhs15*(crhs17 + crhs19);
    rhs[9]=crhs15*(crhs20 + crhs22);
    rhs[10]=crhs23*(crhs19 + crhs24);
    rhs[11]=crhs23*(crhs22 + crhs25);
    rhs[12]=crhs15*(crhs17 + crhs18);
    rhs[13]=crhs15*(crhs20 + crhs21);
    rhs[14]=crhs23*(crhs18 + crhs24);
    rhs[15]=crhs23*(crhs21 + crhs25);

    
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
    const bounded_matrix<double, 2, 2> dlm             = rContactData.DoubleLagrangeMultipliers;
    const double Dt              = rContactData.Dt;
//     const double epsilon  = rContactData.epsilon;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
    
    const bounded_matrix<double, 2, 2> v1 = rContactData.v1;
    const bounded_matrix<double, 2, 2> v2 = rContactData.v2;
    
    const double crhs0 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs1 =     1.0/Dt;
    const double crhs2 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs3 =     N1[0]; // N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs4 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs5 =     N1[1]; // N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs6 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs7 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs8 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs9 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs10 =     (crhs0*crhs3 + crhs4*crhs5)*(crhs3*(Dt*v1(0,0)) + crhs5*(Dt*v1(1,0)) - crhs6*(Dt*v2(0,0)) - crhs7*(Dt*v2(1,0))) + (crhs3*crhs8 + crhs5*crhs9)*(crhs3*(Dt*v1(0,1)) + crhs5*(Dt*v1(1,1)) - crhs6*(Dt*v2(0,1)) - crhs7*(Dt*v2(1,1)));
    const double crhs11 =     crhs1*crhs10*crhs2*Phi[0]; // CRHS1*CRHS10*CRHS2*PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs12 =     -crhs0*crhs11;
    const double crhs13 =     -crhs11*crhs8;
    const double crhs14 =     crhs1*crhs10*crhs2*Phi[1]; // CRHS1*CRHS10*CRHS2*PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs15 =     -crhs14*crhs4;
    const double crhs16 =     -crhs14*crhs9;

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
    const bounded_matrix<double, 2, 2> dlm             = rContactData.DoubleLagrangeMultipliers;
    const double Dt              = rContactData.Dt;
//     const double epsilon  = rContactData.epsilon;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
    
//     const double sign_tangpress = boost::math::sign(augmented_tangent_lm);
    
    const bounded_matrix<double, 2, 2> v1 = rContactData.v1;
    const bounded_matrix<double, 2, 2> v2 = rContactData.v2;
    
    const double crhs0 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs1 =     1.0/Dt;
    const double crhs2 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs3 =     N1[0]; // N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs4 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs5 =     N1[1]; // N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs6 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs7 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs8 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs9 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs10 =     (crhs0*crhs3 + crhs4*crhs5)*(crhs3*(Dt*v1(0,0)) + crhs5*(Dt*v1(1,0)) - crhs6*(Dt*v2(0,0)) - crhs7*(Dt*v2(1,0))) + (crhs3*crhs8 + crhs5*crhs9)*(crhs3*(Dt*v1(0,1)) + crhs5*(Dt*v1(1,1)) - crhs6*(Dt*v2(0,1)) - crhs7*(Dt*v2(1,1)));
    const double crhs11 =     crhs1*crhs10*crhs2*Phi[0]; // CRHS1*CRHS10*CRHS2*PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs12 =     -crhs0*crhs11;
    const double crhs13 =     -crhs11*crhs8;
    const double crhs14 =     crhs1*crhs10*crhs2*Phi[1]; // CRHS1*CRHS10*CRHS2*PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs15 =     -crhs14*crhs4;
    const double crhs16 =     -crhs14*crhs9;

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
//     const bounded_matrix<double, 2, 2> dlm             = rContactData.DoubleLagrangeMultipliers;
//     const double Dt              = rContactData.Dt;
//     const double epsilon  = rContactData.epsilon;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
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
