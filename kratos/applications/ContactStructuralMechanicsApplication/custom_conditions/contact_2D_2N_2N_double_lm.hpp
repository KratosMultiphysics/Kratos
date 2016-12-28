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
    
    static inline bounded_matrix<double,16,16> ComputeGaussPointFrictionlessLHS(
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const ContactData<2,2>& rContactData
        )
{
    bounded_matrix<double,16,16> lhs;
    
    const double epsilon = rContactData.epsilon;
    const double integration_point_gap = inner_prod(rContactData.Gaps,N1);
    const bounded_matrix<double, 2, 2> normalslave = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave   = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm          = rContactData.LagrangeMultipliers;
    const bounded_matrix<double, 2, 2> dlm         = rContactData.DoubleLagrangeMultipliers;

    const std::vector<double> DeltaJs  = rContactData.DeltaJ_s;
    const std::vector<double> DeltaGap = rContactData.DeltaGap;
    const std::vector<array_1d<double,2>> DeltaPhi = rContactData.DeltaPhi;
    const std::vector<array_1d<double,2>> DeltaN1  = rContactData.DeltaN1;
    const std::vector<array_1d<double,2>> DeltaN2  = rContactData.DeltaN2;
    const std::vector<bounded_matrix<double, 2, 2>> DeltaNormals = rContactData.Delta_Normal_s;
    const std::vector<bounded_matrix<double, 2, 2>> Deltatangentxis = rContactData.Delta_Tangent_xi_s;
    
    const double Dt = rContactData.Dt;
    const bounded_matrix<double, 2, 2> v1 = rContactData.v1;
    const bounded_matrix<double, 2, 2> v2 = rContactData.v2;
    
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
    const double clhs125 =     clhs2*clhs3;
    const double clhs126 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs127 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs128 =     DeltaGap[4]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs129 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs130 =     1.0/Dt;
    const double clhs131 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs132 =     Dt*v1(0,0);
    const double clhs133 =     Dt*v1(1,0);
    const double clhs134 =     Dt*v2(0,0);
    const double clhs135 =     Dt*v2(1,0);
    const double clhs136 =     -clhs0*clhs134 + clhs109*clhs133 + clhs132*clhs93 - clhs135*clhs77;
    const double clhs137 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs138 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs139 =     Dt*v1(0,1);
    const double clhs140 =     Dt*v1(1,1);
    const double clhs141 =     Dt*v2(0,1);
    const double clhs142 =     Dt*v2(1,1);
    const double clhs143 =     -clhs0*clhs141 + clhs109*clhs140 + clhs139*clhs93 - clhs142*clhs77;
    const double clhs144 =     clhs109*clhs138 + clhs137*clhs93;
    const double clhs145 =     clhs109*clhs131 + clhs129*clhs93;
    const double clhs146 =     clhs130*(clhs136*(clhs110*clhs131 + clhs129*clhs94) + clhs143*(clhs110*clhs138 + clhs137*clhs94) + clhs144*(-clhs1*clhs141 + clhs110*clhs140 + clhs139*clhs94 - clhs142*clhs78) - clhs145*(clhs0 + clhs1*clhs134 - clhs110*clhs133 - clhs132*clhs94 + clhs135*clhs78));
    const double clhs147 =     clhs125*(-clhs126*clhs128 + clhs129*clhs146);
    const double clhs148 =     DeltaGap[5]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs149 =     clhs130*(clhs136*(clhs111*clhs131 + clhs129*clhs95) + clhs143*(clhs111*clhs138 + clhs137*clhs95) - clhs144*(clhs0 + clhs10*clhs141 - clhs111*clhs140 - clhs139*clhs95 + clhs142*clhs79) + clhs145*(-clhs10*clhs134 + clhs111*clhs133 + clhs132*clhs95 - clhs135*clhs79));
    const double clhs150 =     clhs125*(-clhs126*clhs148 + clhs129*clhs149);
    const double clhs151 =     DeltaGap[6]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs152 =     clhs130*(clhs136*(clhs112*clhs131 + clhs129*clhs96) + clhs143*(clhs112*clhs138 + clhs137*clhs96) + clhs144*(-clhs11*clhs141 + clhs112*clhs140 + clhs139*clhs96 - clhs142*clhs80) - clhs145*(clhs11*clhs134 - clhs112*clhs133 - clhs132*clhs96 + clhs135*clhs80 + clhs77));
    const double clhs153 =     clhs125*(-clhs126*clhs151 + clhs129*clhs152);
    const double clhs154 =     DeltaGap[7]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs155 =     clhs130*(clhs136*(clhs113*clhs131 + clhs129*clhs97) + clhs143*(clhs113*clhs138 + clhs137*clhs97) - clhs144*(-clhs113*clhs140 + clhs12*clhs141 - clhs139*clhs97 + clhs142*clhs81 + clhs77) + clhs145*(clhs113*clhs133 - clhs12*clhs134 + clhs132*clhs97 - clhs135*clhs81));
    const double clhs156 =     clhs125*(-clhs126*clhs154 + clhs129*clhs155);
    const double clhs157 =     clhs5 - clhs6 - clhs7;
    const double clhs158 =     clhs157*clhs3*epsilon;
    const double clhs159 =     clhs13*clhs158;
    const double clhs160 =     clhs157*clhs2*epsilon;
    const double clhs161 =     clhs160*clhs20;
    const double clhs162 =     clhs2*clhs3*epsilon;
    const double clhs163 =     clhs22 - clhs23 - clhs24;
    const double clhs164 =     clhs162*clhs163;
    const double clhs165 =     clhs127*clhs2*clhs3;
    const double clhs166 =     DeltaGap[0]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs167 =     clhs126*clhs2*clhs3;
    const double clhs168 =     clhs126*clhs127*clhs3;
    const double clhs169 =     clhs126*clhs127*clhs2;
    const double clhs170 =     Deltatangentxis[0](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs171 =     clhs136*clhs145 + clhs143*clhs144;
    const double clhs172 =     clhs130*clhs171*clhs2*clhs3;
    const double clhs173 =     clhs129*clhs130*clhs171*clhs3;
    const double clhs174 =     clhs129*clhs130*clhs171*clhs2;
    const double clhs175 =     Deltatangentxis[0](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs176 =     Deltatangentxis[0](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs177 =     Deltatangentxis[0](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs178 =     clhs136*(clhs109*clhs175 + clhs115*clhs131 + clhs129*clhs99 + clhs170*clhs93) + clhs143*(clhs109*clhs177 + clhs115*clhs138 + clhs137*clhs99 + clhs176*clhs93) + clhs144*(clhs115*clhs140 + clhs139*clhs99 - clhs141*clhs15 - clhs142*clhs83) + clhs145*(clhs115*clhs133 + clhs132*clhs99 - clhs134*clhs15 - clhs135*clhs83 + clhs93);
    const double clhs179 =     clhs130*clhs178*clhs2*clhs3;
    const double clhs180 =     -clhs13*clhs168 + clhs13*clhs173 - clhs165*DeltaNormals[0](0,0) - clhs166*clhs167 - clhs169*clhs20 + clhs170*clhs172 + clhs174*clhs20 + clhs179*tan1slave(0,0); // -CLHS13*CLHS168 + CLHS13*CLHS173 - CLHS165*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) - CLHS166*CLHS167 - CLHS169*CLHS20 + CLHS170*CLHS172 + CLHS174*CLHS20 + CLHS179*TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs181 =     clhs158*clhs26;
    const double clhs182 =     clhs160*clhs28;
    const double clhs183 =     clhs30 - clhs31 - clhs32;
    const double clhs184 =     clhs162*clhs183;
    const double clhs185 =     DeltaGap[1]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs186 =     Deltatangentxis[1](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs187 =     Deltatangentxis[1](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs188 =     Deltatangentxis[1](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs189 =     Deltatangentxis[1](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs190 =     clhs136*(clhs102*clhs129 + clhs109*clhs187 + clhs118*clhs131 + clhs186*clhs93) + clhs143*(clhs102*clhs137 + clhs109*clhs189 + clhs118*clhs138 + clhs188*clhs93) + clhs144*(clhs102*clhs139 + clhs118*clhs140 - clhs141*clhs27 - clhs142*clhs86 + clhs93) + clhs145*(clhs102*clhs132 + clhs118*clhs133 - clhs134*clhs27 - clhs135*clhs86);
    const double clhs191 =     clhs130*clhs190*clhs2*clhs3;
    const double clhs192 =     -clhs165*DeltaNormals[1](0,0) - clhs167*clhs185 - clhs168*clhs26 - clhs169*clhs28 + clhs172*clhs186 + clhs173*clhs26 + clhs174*clhs28 + clhs191*tan1slave(0,0); // -CLHS165*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) - CLHS167*CLHS185 - CLHS168*CLHS26 - CLHS169*CLHS28 + CLHS172*CLHS186 + CLHS173*CLHS26 + CLHS174*CLHS28 + CLHS191*TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs193 =     clhs158*clhs34;
    const double clhs194 =     clhs160*clhs36;
    const double clhs195 =     clhs38 - clhs39 - clhs40;
    const double clhs196 =     clhs162*clhs195;
    const double clhs197 =     DeltaGap[2]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs198 =     Deltatangentxis[2](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs199 =     Deltatangentxis[2](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs200 =     Deltatangentxis[2](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs201 =     Deltatangentxis[2](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs202 =     clhs136*(clhs103*clhs129 + clhs109*clhs199 + clhs119*clhs131 + clhs198*clhs93) + clhs143*(clhs103*clhs137 + clhs109*clhs201 + clhs119*clhs138 + clhs200*clhs93) + clhs144*(clhs103*clhs139 + clhs119*clhs140 - clhs141*clhs35 - clhs142*clhs87) + clhs145*(clhs103*clhs132 + clhs109 + clhs119*clhs133 - clhs134*clhs35 - clhs135*clhs87);
    const double clhs203 =     clhs130*clhs2*clhs202*clhs3;
    const double clhs204 =     -clhs165*DeltaNormals[2](0,0) - clhs167*clhs197 - clhs168*clhs34 - clhs169*clhs36 + clhs172*clhs198 + clhs173*clhs34 + clhs174*clhs36 + clhs203*tan1slave(0,0); // -CLHS165*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) - CLHS167*CLHS197 - CLHS168*CLHS34 - CLHS169*CLHS36 + CLHS172*CLHS198 + CLHS173*CLHS34 + CLHS174*CLHS36 + CLHS203*TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs205 =     clhs158*clhs42;
    const double clhs206 =     clhs160*clhs44;
    const double clhs207 =     clhs46 - clhs47 - clhs48;
    const double clhs208 =     clhs162*clhs207;
    const double clhs209 =     DeltaGap[3]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs210 =     Deltatangentxis[3](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs211 =     Deltatangentxis[3](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs212 =     Deltatangentxis[3](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs213 =     Deltatangentxis[3](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs214 =     clhs136*(clhs104*clhs129 + clhs109*clhs211 + clhs120*clhs131 + clhs210*clhs93) + clhs143*(clhs104*clhs137 + clhs109*clhs213 + clhs120*clhs138 + clhs212*clhs93) + clhs144*(clhs104*clhs139 + clhs109 + clhs120*clhs140 - clhs141*clhs43 - clhs142*clhs88) + clhs145*(clhs104*clhs132 + clhs120*clhs133 - clhs134*clhs43 - clhs135*clhs88);
    const double clhs215 =     clhs130*clhs2*clhs214*clhs3;
    const double clhs216 =     -clhs165*DeltaNormals[3](0,0) - clhs167*clhs209 - clhs168*clhs42 - clhs169*clhs44 + clhs172*clhs210 + clhs173*clhs42 + clhs174*clhs44 + clhs215*tan1slave(0,0); // -CLHS165*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) - CLHS167*CLHS209 - CLHS168*CLHS42 - CLHS169*CLHS44 + CLHS172*CLHS210 + CLHS173*CLHS42 + CLHS174*CLHS44 + CLHS215*TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs217 =     clhs2*epsilon;
    const double clhs218 =     clhs217*std::pow(clhs3, 2);
    const double clhs219 =     -clhs218;
    const double clhs220 =     clhs162*clhs4;
    const double clhs221 =     -clhs220;
    const double clhs222 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs223 =     clhs125*(-clhs128*clhs222 + clhs137*clhs146);
    const double clhs224 =     clhs125*(clhs137*clhs149 - clhs148*clhs222);
    const double clhs225 =     clhs125*(clhs137*clhs152 - clhs151*clhs222);
    const double clhs226 =     clhs125*(clhs137*clhs155 - clhs154*clhs222);
    const double clhs227 =     clhs52 - clhs53 - clhs54;
    const double clhs228 =     clhs227*clhs3*epsilon;
    const double clhs229 =     clhs13*clhs228;
    const double clhs230 =     clhs2*clhs227*epsilon;
    const double clhs231 =     clhs20*clhs230;
    const double clhs232 =     clhs61 - clhs62 - clhs63;
    const double clhs233 =     clhs162*clhs232;
    const double clhs234 =     clhs2*clhs222*clhs3;
    const double clhs235 =     clhs127*clhs222*clhs3;
    const double clhs236 =     clhs127*clhs2*clhs222;
    const double clhs237 =     clhs130*clhs137*clhs171*clhs3;
    const double clhs238 =     clhs130*clhs137*clhs171*clhs2;
    const double clhs239 =     -clhs13*clhs235 + clhs13*clhs237 - clhs165*DeltaNormals[0](0,1) - clhs166*clhs234 + clhs172*clhs176 + clhs179*tan1slave(0,1) - clhs20*clhs236 + clhs20*clhs238; // -CLHS13*CLHS235 + CLHS13*CLHS237 - CLHS165*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) - CLHS166*CLHS234 + CLHS172*CLHS176 + CLHS179*TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS20*CLHS236 + CLHS20*CLHS238
    const double clhs240 =     clhs228*clhs26;
    const double clhs241 =     clhs230*clhs28;
    const double clhs242 =     clhs65 - clhs66 - clhs67;
    const double clhs243 =     clhs162*clhs242;
    const double clhs244 =     -clhs165*DeltaNormals[1](0,1) + clhs172*clhs188 - clhs185*clhs234 + clhs191*tan1slave(0,1) - clhs235*clhs26 - clhs236*clhs28 + clhs237*clhs26 + clhs238*clhs28; // -CLHS165*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS172*CLHS188 - CLHS185*CLHS234 + CLHS191*TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS235*CLHS26 - CLHS236*CLHS28 + CLHS237*CLHS26 + CLHS238*CLHS28
    const double clhs245 =     clhs228*clhs34;
    const double clhs246 =     clhs230*clhs36;
    const double clhs247 =     clhs69 - clhs70 - clhs71;
    const double clhs248 =     clhs162*clhs247;
    const double clhs249 =     -clhs165*DeltaNormals[2](0,1) + clhs172*clhs200 - clhs197*clhs234 + clhs203*tan1slave(0,1) - clhs235*clhs34 - clhs236*clhs36 + clhs237*clhs34 + clhs238*clhs36; // -CLHS165*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS172*CLHS200 - CLHS197*CLHS234 + CLHS203*TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS235*CLHS34 - CLHS236*CLHS36 + CLHS237*CLHS34 + CLHS238*CLHS36
    const double clhs250 =     clhs228*clhs42;
    const double clhs251 =     clhs230*clhs44;
    const double clhs252 =     clhs73 - clhs74 - clhs75;
    const double clhs253 =     clhs162*clhs252;
    const double clhs254 =     -clhs165*DeltaNormals[3](0,1) + clhs172*clhs212 - clhs209*clhs234 + clhs215*tan1slave(0,1) - clhs235*clhs42 - clhs236*clhs44 + clhs237*clhs42 + clhs238*clhs44; // -CLHS165*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS172*CLHS212 - CLHS209*CLHS234 + CLHS215*TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS235*CLHS42 - CLHS236*CLHS44 + CLHS237*CLHS42 + CLHS238*CLHS44
    const double clhs255 =     clhs2*clhs4;
    const double clhs256 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs257 =     clhs255*(-clhs128*clhs256 + clhs131*clhs146);
    const double clhs258 =     clhs255*(clhs131*clhs149 - clhs148*clhs256);
    const double clhs259 =     clhs255*(clhs131*clhs152 - clhs151*clhs256);
    const double clhs260 =     clhs255*(clhs131*clhs155 - clhs154*clhs256);
    const double clhs261 =     clhs157*clhs4*epsilon;
    const double clhs262 =     clhs13*clhs261;
    const double clhs263 =     clhs160*clhs21;
    const double clhs264 =     clhs2*clhs4*epsilon;
    const double clhs265 =     clhs163*clhs264;
    const double clhs266 =     clhs127*clhs2*clhs4;
    const double clhs267 =     clhs2*clhs256*clhs4;
    const double clhs268 =     clhs127*clhs256*clhs4;
    const double clhs269 =     clhs127*clhs2*clhs256;
    const double clhs270 =     clhs130*clhs171*clhs2*clhs4;
    const double clhs271 =     clhs130*clhs131*clhs171*clhs4;
    const double clhs272 =     clhs130*clhs131*clhs171*clhs2;
    const double clhs273 =     clhs130*clhs178*clhs2*clhs4;
    const double clhs274 =     -clhs13*clhs268 + clhs13*clhs271 - clhs166*clhs267 + clhs175*clhs270 - clhs21*clhs269 + clhs21*clhs272 - clhs266*DeltaNormals[0](1,0) + clhs273*tan1slave(1,0); // -CLHS13*CLHS268 + CLHS13*CLHS271 - CLHS166*CLHS267 + CLHS175*CLHS270 - CLHS21*CLHS269 + CLHS21*CLHS272 - CLHS266*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS273*TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs275 =     clhs26*clhs261;
    const double clhs276 =     clhs160*clhs29;
    const double clhs277 =     clhs183*clhs264;
    const double clhs278 =     clhs130*clhs190*clhs2*clhs4;
    const double clhs279 =     -clhs185*clhs267 + clhs187*clhs270 - clhs26*clhs268 + clhs26*clhs271 - clhs266*DeltaNormals[1](1,0) - clhs269*clhs29 + clhs272*clhs29 + clhs278*tan1slave(1,0); // -CLHS185*CLHS267 + CLHS187*CLHS270 - CLHS26*CLHS268 + CLHS26*CLHS271 - CLHS266*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) - CLHS269*CLHS29 + CLHS272*CLHS29 + CLHS278*TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs280 =     clhs261*clhs34;
    const double clhs281 =     clhs160*clhs37;
    const double clhs282 =     clhs195*clhs264;
    const double clhs283 =     clhs130*clhs2*clhs202*clhs4;
    const double clhs284 =     -clhs197*clhs267 + clhs199*clhs270 - clhs266*DeltaNormals[2](1,0) - clhs268*clhs34 - clhs269*clhs37 + clhs271*clhs34 + clhs272*clhs37 + clhs283*tan1slave(1,0); // -CLHS197*CLHS267 + CLHS199*CLHS270 - CLHS266*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) - CLHS268*CLHS34 - CLHS269*CLHS37 + CLHS271*CLHS34 + CLHS272*CLHS37 + CLHS283*TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs285 =     clhs261*clhs42;
    const double clhs286 =     clhs160*clhs45;
    const double clhs287 =     clhs207*clhs264;
    const double clhs288 =     clhs130*clhs2*clhs214*clhs4;
    const double clhs289 =     -clhs209*clhs267 + clhs211*clhs270 - clhs266*DeltaNormals[3](1,0) - clhs268*clhs42 - clhs269*clhs45 + clhs271*clhs42 + clhs272*clhs45 + clhs288*tan1slave(1,0); // -CLHS209*CLHS267 + CLHS211*CLHS270 - CLHS266*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) - CLHS268*CLHS42 - CLHS269*CLHS45 + CLHS271*CLHS42 + CLHS272*CLHS45 + CLHS288*TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs290 =     clhs217*std::pow(clhs4, 2);
    const double clhs291 =     -clhs290;
    const double clhs292 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs293 =     clhs255*(-clhs128*clhs292 + clhs138*clhs146);
    const double clhs294 =     clhs255*(clhs138*clhs149 - clhs148*clhs292);
    const double clhs295 =     clhs255*(clhs138*clhs152 - clhs151*clhs292);
    const double clhs296 =     clhs255*(clhs138*clhs155 - clhs154*clhs292);
    const double clhs297 =     clhs227*clhs4*epsilon;
    const double clhs298 =     clhs13*clhs297;
    const double clhs299 =     clhs21*clhs230;
    const double clhs300 =     clhs232*clhs264;
    const double clhs301 =     clhs2*clhs292*clhs4;
    const double clhs302 =     clhs127*clhs292*clhs4;
    const double clhs303 =     clhs127*clhs2*clhs292;
    const double clhs304 =     clhs130*clhs138*clhs171*clhs4;
    const double clhs305 =     clhs130*clhs138*clhs171*clhs2;
    const double clhs306 =     -clhs13*clhs302 + clhs13*clhs304 - clhs166*clhs301 + clhs177*clhs270 - clhs21*clhs303 + clhs21*clhs305 - clhs266*DeltaNormals[0](1,1) + clhs273*tan1slave(1,1); // -CLHS13*CLHS302 + CLHS13*CLHS304 - CLHS166*CLHS301 + CLHS177*CLHS270 - CLHS21*CLHS303 + CLHS21*CLHS305 - CLHS266*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS273*TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs307 =     clhs26*clhs297;
    const double clhs308 =     clhs230*clhs29;
    const double clhs309 =     clhs242*clhs264;
    const double clhs310 =     -clhs185*clhs301 + clhs189*clhs270 - clhs26*clhs302 + clhs26*clhs304 - clhs266*DeltaNormals[1](1,1) + clhs278*tan1slave(1,1) - clhs29*clhs303 + clhs29*clhs305; // -CLHS185*CLHS301 + CLHS189*CLHS270 - CLHS26*CLHS302 + CLHS26*CLHS304 - CLHS266*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS278*TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS29*CLHS303 + CLHS29*CLHS305
    const double clhs311 =     clhs297*clhs34;
    const double clhs312 =     clhs230*clhs37;
    const double clhs313 =     clhs247*clhs264;
    const double clhs314 =     -clhs197*clhs301 + clhs201*clhs270 - clhs266*DeltaNormals[2](1,1) + clhs283*tan1slave(1,1) - clhs302*clhs34 - clhs303*clhs37 + clhs304*clhs34 + clhs305*clhs37; // -CLHS197*CLHS301 + CLHS201*CLHS270 - CLHS266*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS283*TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS302*CLHS34 - CLHS303*CLHS37 + CLHS304*CLHS34 + CLHS305*CLHS37
    const double clhs315 =     clhs297*clhs42;
    const double clhs316 =     clhs230*clhs45;
    const double clhs317 =     clhs252*clhs264;
    const double clhs318 =     -clhs209*clhs301 + clhs213*clhs270 - clhs266*DeltaNormals[3](1,1) + clhs288*tan1slave(1,1) - clhs302*clhs42 - clhs303*clhs45 + clhs304*clhs42 + clhs305*clhs45; // -CLHS209*CLHS301 + CLHS213*CLHS270 - CLHS266*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS288*TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS302*CLHS42 - CLHS303*CLHS45 + CLHS304*CLHS42 + CLHS305*CLHS45

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
    lhs(8,0)=clhs147;
    lhs(8,1)=clhs150;
    lhs(8,2)=clhs153;
    lhs(8,3)=clhs156;
    lhs(8,4)=clhs159 + clhs161 + clhs164 + clhs180;
    lhs(8,5)=clhs181 + clhs182 + clhs184 + clhs192;
    lhs(8,6)=clhs193 + clhs194 + clhs196 + clhs204;
    lhs(8,7)=clhs205 + clhs206 + clhs208 + clhs216;
    lhs(8,8)=clhs219;
    lhs(8,9)=0;
    lhs(8,10)=clhs221;
    lhs(8,11)=0;
    lhs(8,12)=clhs218;
    lhs(8,13)=0;
    lhs(8,14)=clhs220;
    lhs(8,15)=0;
    lhs(9,0)=clhs223;
    lhs(9,1)=clhs224;
    lhs(9,2)=clhs225;
    lhs(9,3)=clhs226;
    lhs(9,4)=clhs229 + clhs231 + clhs233 + clhs239;
    lhs(9,5)=clhs240 + clhs241 + clhs243 + clhs244;
    lhs(9,6)=clhs245 + clhs246 + clhs248 + clhs249;
    lhs(9,7)=clhs250 + clhs251 + clhs253 + clhs254;
    lhs(9,8)=0;
    lhs(9,9)=clhs219;
    lhs(9,10)=0;
    lhs(9,11)=clhs221;
    lhs(9,12)=0;
    lhs(9,13)=clhs218;
    lhs(9,14)=0;
    lhs(9,15)=clhs220;
    lhs(10,0)=clhs257;
    lhs(10,1)=clhs258;
    lhs(10,2)=clhs259;
    lhs(10,3)=clhs260;
    lhs(10,4)=clhs262 + clhs263 + clhs265 + clhs274;
    lhs(10,5)=clhs275 + clhs276 + clhs277 + clhs279;
    lhs(10,6)=clhs280 + clhs281 + clhs282 + clhs284;
    lhs(10,7)=clhs285 + clhs286 + clhs287 + clhs289;
    lhs(10,8)=clhs221;
    lhs(10,9)=0;
    lhs(10,10)=clhs291;
    lhs(10,11)=0;
    lhs(10,12)=clhs220;
    lhs(10,13)=0;
    lhs(10,14)=clhs290;
    lhs(10,15)=0;
    lhs(11,0)=clhs293;
    lhs(11,1)=clhs294;
    lhs(11,2)=clhs295;
    lhs(11,3)=clhs296;
    lhs(11,4)=clhs298 + clhs299 + clhs300 + clhs306;
    lhs(11,5)=clhs307 + clhs308 + clhs309 + clhs310;
    lhs(11,6)=clhs311 + clhs312 + clhs313 + clhs314;
    lhs(11,7)=clhs315 + clhs316 + clhs317 + clhs318;
    lhs(11,8)=0;
    lhs(11,9)=clhs221;
    lhs(11,10)=0;
    lhs(11,11)=clhs291;
    lhs(11,12)=0;
    lhs(11,13)=clhs220;
    lhs(11,14)=0;
    lhs(11,15)=clhs290;
    lhs(12,0)=clhs147;
    lhs(12,1)=clhs150;
    lhs(12,2)=clhs153;
    lhs(12,3)=clhs156;
    lhs(12,4)=-clhs159 - clhs161 - clhs164 + clhs180;
    lhs(12,5)=-clhs181 - clhs182 - clhs184 + clhs192;
    lhs(12,6)=-clhs193 - clhs194 - clhs196 + clhs204;
    lhs(12,7)=-clhs205 - clhs206 - clhs208 + clhs216;
    lhs(12,8)=clhs218;
    lhs(12,9)=0;
    lhs(12,10)=clhs220;
    lhs(12,11)=0;
    lhs(12,12)=clhs219;
    lhs(12,13)=0;
    lhs(12,14)=clhs221;
    lhs(12,15)=0;
    lhs(13,0)=clhs223;
    lhs(13,1)=clhs224;
    lhs(13,2)=clhs225;
    lhs(13,3)=clhs226;
    lhs(13,4)=-clhs229 - clhs231 - clhs233 + clhs239;
    lhs(13,5)=-clhs240 - clhs241 - clhs243 + clhs244;
    lhs(13,6)=-clhs245 - clhs246 - clhs248 + clhs249;
    lhs(13,7)=-clhs250 - clhs251 - clhs253 + clhs254;
    lhs(13,8)=0;
    lhs(13,9)=clhs218;
    lhs(13,10)=0;
    lhs(13,11)=clhs220;
    lhs(13,12)=0;
    lhs(13,13)=clhs219;
    lhs(13,14)=0;
    lhs(13,15)=clhs221;
    lhs(14,0)=clhs257;
    lhs(14,1)=clhs258;
    lhs(14,2)=clhs259;
    lhs(14,3)=clhs260;
    lhs(14,4)=-clhs262 - clhs263 - clhs265 + clhs274;
    lhs(14,5)=-clhs275 - clhs276 - clhs277 + clhs279;
    lhs(14,6)=-clhs280 - clhs281 - clhs282 + clhs284;
    lhs(14,7)=-clhs285 - clhs286 - clhs287 + clhs289;
    lhs(14,8)=clhs220;
    lhs(14,9)=0;
    lhs(14,10)=clhs290;
    lhs(14,11)=0;
    lhs(14,12)=clhs221;
    lhs(14,13)=0;
    lhs(14,14)=clhs291;
    lhs(14,15)=0;
    lhs(15,0)=clhs293;
    lhs(15,1)=clhs294;
    lhs(15,2)=clhs295;
    lhs(15,3)=clhs296;
    lhs(15,4)=-clhs298 - clhs299 - clhs300 + clhs306;
    lhs(15,5)=-clhs307 - clhs308 - clhs309 + clhs310;
    lhs(15,6)=-clhs311 - clhs312 - clhs313 + clhs314;
    lhs(15,7)=-clhs315 - clhs316 - clhs317 + clhs318;
    lhs(15,8)=0;
    lhs(15,9)=clhs220;
    lhs(15,10)=0;
    lhs(15,11)=clhs290;
    lhs(15,12)=0;
    lhs(15,13)=clhs221;
    lhs(15,14)=0;
    lhs(15,15)=clhs291;

    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,16,16> ComputeGaussPointLHS(
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const double mu, 
        const ContactData<2,2>& rContactData
        )
{
    bounded_matrix<double,16,16> lhs;
    
    const double epsilon = rContactData.epsilon;
    const double integration_point_gap = inner_prod(rContactData.Gaps,N1);
    const bounded_matrix<double, 2, 2> normalslave = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave   = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm          = rContactData.LagrangeMultipliers;
    const bounded_matrix<double, 2, 2> dlm         = rContactData.DoubleLagrangeMultipliers;
    
    const std::vector<double> DeltaJs = rContactData.DeltaJ_s;
    const std::vector<double> DeltaGap = rContactData.DeltaGap;
    const std::vector<array_1d<double,2>> DeltaPhi = rContactData.DeltaPhi;
    const std::vector<array_1d<double,2>> DeltaN1  = rContactData.DeltaN1;
    const std::vector<array_1d<double,2>> DeltaN2  = rContactData.DeltaN2;
    const std::vector<bounded_matrix<double, 2, 2>> DeltaNormals    = rContactData.Delta_Normal_s;
    const std::vector<bounded_matrix<double, 2, 2>> Deltatangentxis = rContactData.Delta_Tangent_xi_s;
    
    const array_1d<double, 1> Ctan = rContactData.Ctan;
    const std::vector<array_1d<double, 1>> DeltaCtan = rContactData.DeltaCtan;
    
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
    const double clhs19 =     DeltaPhi[0][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs20 =     DeltaPhi[0][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs21 =     clhs19*dlm(0,0) + clhs20*dlm(1,0);
    const double clhs22 =     clhs2*clhs21;
    const double clhs23 =     clhs19*lm(0,0);
    const double clhs24 =     clhs20*lm(1,0);
    const double clhs25 =     clhs2*(clhs23 + clhs24);
    const double clhs26 =     DeltaJs[1]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs27 =     DeltaN2[1][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs28 =     DeltaPhi[1][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs29 =     DeltaPhi[1][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs30 =     clhs28*dlm(0,0) + clhs29*dlm(1,0);
    const double clhs31 =     clhs2*clhs30;
    const double clhs32 =     clhs28*lm(0,0);
    const double clhs33 =     clhs29*lm(1,0);
    const double clhs34 =     clhs2*(clhs32 + clhs33);
    const double clhs35 =     DeltaJs[2]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs36 =     DeltaN2[2][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs37 =     DeltaPhi[2][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs38 =     DeltaPhi[2][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs39 =     clhs37*dlm(0,0) + clhs38*dlm(1,0);
    const double clhs40 =     clhs2*clhs39;
    const double clhs41 =     clhs37*lm(0,0);
    const double clhs42 =     clhs38*lm(1,0);
    const double clhs43 =     clhs2*(clhs41 + clhs42);
    const double clhs44 =     DeltaJs[3]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs45 =     DeltaN2[3][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs46 =     DeltaPhi[3][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs47 =     DeltaPhi[3][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs48 =     clhs46*dlm(0,0) + clhs47*dlm(1,0);
    const double clhs49 =     clhs2*clhs48;
    const double clhs50 =     clhs46*lm(0,0);
    const double clhs51 =     clhs47*lm(1,0);
    const double clhs52 =     clhs2*(clhs50 + clhs51);
    const double clhs53 =     clhs2*clhs3;
    const double clhs54 =     clhs0*clhs53;
    const double clhs55 =     clhs2*clhs4;
    const double clhs56 =     clhs0*clhs55;
    const double clhs57 =     clhs3*dlm(0,1) + clhs4*dlm(1,1);
    const double clhs58 =     clhs3*lm(0,1);
    const double clhs59 =     clhs4*lm(1,1);
    const double clhs60 =     clhs58 + clhs59;
    const double clhs61 =     clhs2*(clhs57 + clhs60);
    const double clhs62 =     clhs0*clhs57;
    const double clhs63 =     clhs2*clhs57;
    const double clhs64 =     clhs0*clhs60;
    const double clhs65 =     clhs2*clhs60;
    const double clhs66 =     clhs19*dlm(0,1) + clhs20*dlm(1,1);
    const double clhs67 =     clhs2*clhs66;
    const double clhs68 =     clhs19*lm(0,1);
    const double clhs69 =     clhs20*lm(1,1);
    const double clhs70 =     clhs2*(clhs68 + clhs69);
    const double clhs71 =     clhs28*dlm(0,1) + clhs29*dlm(1,1);
    const double clhs72 =     clhs2*clhs71;
    const double clhs73 =     clhs28*lm(0,1);
    const double clhs74 =     clhs29*lm(1,1);
    const double clhs75 =     clhs2*(clhs73 + clhs74);
    const double clhs76 =     clhs37*dlm(0,1) + clhs38*dlm(1,1);
    const double clhs77 =     clhs2*clhs76;
    const double clhs78 =     clhs37*lm(0,1);
    const double clhs79 =     clhs38*lm(1,1);
    const double clhs80 =     clhs2*(clhs78 + clhs79);
    const double clhs81 =     clhs46*dlm(0,1) + clhs47*dlm(1,1);
    const double clhs82 =     clhs2*clhs81;
    const double clhs83 =     clhs46*lm(0,1);
    const double clhs84 =     clhs47*lm(1,1);
    const double clhs85 =     clhs2*(clhs83 + clhs84);
    const double clhs86 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs87 =     DeltaN2[4][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs88 =     DeltaN2[5][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs89 =     DeltaN2[6][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs90 =     DeltaN2[7][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs91 =     clhs5*clhs86;
    const double clhs92 =     DeltaN2[0][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs93 =     clhs8*clhs86;
    const double clhs94 =     DeltaN2[1][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs95 =     DeltaN2[2][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs96 =     DeltaN2[3][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs97 =     clhs53*clhs86;
    const double clhs98 =     clhs55*clhs86;
    const double clhs99 =     clhs57*clhs86;
    const double clhs100 =     clhs60*clhs86;
    const double clhs101 =     N1[0]; // N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs102 =     DeltaN1[4][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs103 =     DeltaN1[5][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs104 =     DeltaN1[6][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs105 =     DeltaN1[7][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs106 =     clhs101*clhs5;
    const double clhs107 =     DeltaN1[0][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs108 =     clhs101*clhs8;
    const double clhs109 =     DeltaN1[1][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs110 =     DeltaN1[2][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs111 =     DeltaN1[3][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs112 =     -clhs101*clhs53;
    const double clhs113 =     -clhs101*clhs55;
    const double clhs114 =     clhs101*clhs57;
    const double clhs115 =     clhs101*clhs60;
    const double clhs116 =     N1[1]; // N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs117 =     DeltaN1[4][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs118 =     DeltaN1[5][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs119 =     DeltaN1[6][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs120 =     DeltaN1[7][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs121 =     clhs116*clhs5;
    const double clhs122 =     DeltaN1[0][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs123 =     clhs116*clhs8;
    const double clhs124 =     DeltaN1[1][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs125 =     DeltaN1[2][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs126 =     DeltaN1[3][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs127 =     -clhs116*clhs53;
    const double clhs128 =     -clhs116*clhs55;
    const double clhs129 =     clhs116*clhs57;
    const double clhs130 =     clhs116*clhs60;
    const double clhs131 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs132 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs133 =     DeltaGap[4]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs134 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs135 =     Ctan[0]; // CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1))
    const double clhs136 =     DeltaCtan[4][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U2(0,0))
    const double clhs137 =     clhs53*(-clhs131*clhs133 + clhs134*clhs136);
    const double clhs138 =     DeltaGap[5]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs139 =     DeltaCtan[5][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U2(0,1))
    const double clhs140 =     clhs53*(-clhs131*clhs138 + clhs134*clhs139);
    const double clhs141 =     DeltaGap[6]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs142 =     DeltaCtan[6][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U2(1,0))
    const double clhs143 =     clhs53*(-clhs131*clhs141 + clhs134*clhs142);
    const double clhs144 =     DeltaGap[7]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs145 =     DeltaCtan[7][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U2(1,1))
    const double clhs146 =     clhs53*(-clhs131*clhs144 + clhs134*clhs145);
    const double clhs147 =     clhs5 - clhs6 - clhs7;
    const double clhs148 =     clhs147*clhs3*epsilon;
    const double clhs149 =     clhs13*clhs148;
    const double clhs150 =     clhs147*clhs2*epsilon;
    const double clhs151 =     clhs150*clhs19;
    const double clhs152 =     clhs3*epsilon;
    const double clhs153 =     clhs2*(clhs21 - clhs23 - clhs24);
    const double clhs154 =     clhs152*clhs153;
    const double clhs155 =     clhs132*clhs2*clhs3;
    const double clhs156 =     DeltaGap[0]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs157 =     clhs131*clhs2*clhs3;
    const double clhs158 =     clhs131*clhs132*clhs3;
    const double clhs159 =     clhs131*clhs132*clhs2;
    const double clhs160 =     clhs135*clhs2*clhs3;
    const double clhs161 =     clhs134*clhs135*clhs3;
    const double clhs162 =     clhs134*clhs135*clhs2;
    const double clhs163 =     DeltaCtan[0][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U1(0,0))
    const double clhs164 =     clhs134*clhs2*clhs3;
    const double clhs165 =     -clhs13*clhs158 + clhs13*clhs161 - clhs155*DeltaNormals[0](0,0) - clhs156*clhs157 - clhs159*clhs19 + clhs160*Deltatangentxis[0](0,0) + clhs162*clhs19 + clhs163*clhs164; // -CLHS13*CLHS158 + CLHS13*CLHS161 - CLHS155*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) - CLHS156*CLHS157 - CLHS159*CLHS19 + CLHS160*DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS162*CLHS19 + CLHS163*CLHS164
    const double clhs166 =     clhs148*clhs26;
    const double clhs167 =     clhs150*clhs28;
    const double clhs168 =     clhs2*(clhs30 - clhs32 - clhs33);
    const double clhs169 =     clhs152*clhs168;
    const double clhs170 =     DeltaGap[1]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs171 =     DeltaCtan[1][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U1(0,1))
    const double clhs172 =     -clhs155*DeltaNormals[1](0,0) - clhs157*clhs170 - clhs158*clhs26 - clhs159*clhs28 + clhs160*Deltatangentxis[1](0,0) + clhs161*clhs26 + clhs162*clhs28 + clhs164*clhs171; // -CLHS155*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) - CLHS157*CLHS170 - CLHS158*CLHS26 - CLHS159*CLHS28 + CLHS160*DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS161*CLHS26 + CLHS162*CLHS28 + CLHS164*CLHS171
    const double clhs173 =     clhs148*clhs35;
    const double clhs174 =     clhs150*clhs37;
    const double clhs175 =     clhs2*(clhs39 - clhs41 - clhs42);
    const double clhs176 =     clhs152*clhs175;
    const double clhs177 =     DeltaGap[2]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs178 =     DeltaCtan[2][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U1(1,0))
    const double clhs179 =     -clhs155*DeltaNormals[2](0,0) - clhs157*clhs177 - clhs158*clhs35 - clhs159*clhs37 + clhs160*Deltatangentxis[2](0,0) + clhs161*clhs35 + clhs162*clhs37 + clhs164*clhs178; // -CLHS155*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) - CLHS157*CLHS177 - CLHS158*CLHS35 - CLHS159*CLHS37 + CLHS160*DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS161*CLHS35 + CLHS162*CLHS37 + CLHS164*CLHS178
    const double clhs180 =     clhs148*clhs44;
    const double clhs181 =     clhs150*clhs46;
    const double clhs182 =     clhs2*(clhs48 - clhs50 - clhs51);
    const double clhs183 =     clhs152*clhs182;
    const double clhs184 =     DeltaGap[3]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs185 =     DeltaCtan[3][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U1(1,1))
    const double clhs186 =     -clhs155*DeltaNormals[3](0,0) - clhs157*clhs184 - clhs158*clhs44 - clhs159*clhs46 + clhs160*Deltatangentxis[3](0,0) + clhs161*clhs44 + clhs162*clhs46 + clhs164*clhs185; // -CLHS155*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) - CLHS157*CLHS184 - CLHS158*CLHS44 - CLHS159*CLHS46 + CLHS160*DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS161*CLHS44 + CLHS162*CLHS46 + CLHS164*CLHS185
    const double clhs187 =     -clhs152;
    const double clhs188 =     DeltaCtan[8][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), LM(0,0))
    const double clhs189 =     clhs134*clhs188;
    const double clhs190 =     DeltaCtan[9][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), LM(0,1))
    const double clhs191 =     clhs164*clhs190;
    const double clhs192 =     clhs4*epsilon;
    const double clhs193 =     -clhs192;
    const double clhs194 =     DeltaCtan[10][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), LM(1,0))
    const double clhs195 =     clhs134*clhs194;
    const double clhs196 =     DeltaCtan[11][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), LM(1,1))
    const double clhs197 =     clhs164*clhs196;
    const double clhs198 =     clhs2*epsilon;
    const double clhs199 =     clhs198*std::pow(clhs3, 2);
    const double clhs200 =     clhs152*clhs55;
    const double clhs201 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs202 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs203 =     clhs53*(-clhs133*clhs201 + clhs136*clhs202);
    const double clhs204 =     clhs53*(-clhs138*clhs201 + clhs139*clhs202);
    const double clhs205 =     clhs53*(-clhs141*clhs201 + clhs142*clhs202);
    const double clhs206 =     clhs53*(-clhs144*clhs201 + clhs145*clhs202);
    const double clhs207 =     clhs57 - clhs58 - clhs59;
    const double clhs208 =     clhs207*clhs3*epsilon;
    const double clhs209 =     clhs13*clhs208;
    const double clhs210 =     clhs2*clhs207*epsilon;
    const double clhs211 =     clhs19*clhs210;
    const double clhs212 =     clhs2*(clhs66 - clhs68 - clhs69);
    const double clhs213 =     clhs152*clhs212;
    const double clhs214 =     clhs2*clhs201*clhs3;
    const double clhs215 =     clhs132*clhs201*clhs3;
    const double clhs216 =     clhs132*clhs2*clhs201;
    const double clhs217 =     clhs135*clhs202*clhs3;
    const double clhs218 =     clhs135*clhs2*clhs202;
    const double clhs219 =     clhs2*clhs202*clhs3;
    const double clhs220 =     -clhs13*clhs215 + clhs13*clhs217 - clhs155*DeltaNormals[0](0,1) - clhs156*clhs214 + clhs160*Deltatangentxis[0](0,1) + clhs163*clhs219 - clhs19*clhs216 + clhs19*clhs218; // -CLHS13*CLHS215 + CLHS13*CLHS217 - CLHS155*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) - CLHS156*CLHS214 + CLHS160*DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS163*CLHS219 - CLHS19*CLHS216 + CLHS19*CLHS218
    const double clhs221 =     clhs208*clhs26;
    const double clhs222 =     clhs210*clhs28;
    const double clhs223 =     clhs2*(clhs71 - clhs73 - clhs74);
    const double clhs224 =     clhs152*clhs223;
    const double clhs225 =     -clhs155*DeltaNormals[1](0,1) + clhs160*Deltatangentxis[1](0,1) - clhs170*clhs214 + clhs171*clhs219 - clhs215*clhs26 - clhs216*clhs28 + clhs217*clhs26 + clhs218*clhs28; // -CLHS155*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS160*DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) - CLHS170*CLHS214 + CLHS171*CLHS219 - CLHS215*CLHS26 - CLHS216*CLHS28 + CLHS217*CLHS26 + CLHS218*CLHS28
    const double clhs226 =     clhs208*clhs35;
    const double clhs227 =     clhs210*clhs37;
    const double clhs228 =     clhs2*(clhs76 - clhs78 - clhs79);
    const double clhs229 =     clhs152*clhs228;
    const double clhs230 =     -clhs155*DeltaNormals[2](0,1) + clhs160*Deltatangentxis[2](0,1) - clhs177*clhs214 + clhs178*clhs219 - clhs215*clhs35 - clhs216*clhs37 + clhs217*clhs35 + clhs218*clhs37; // -CLHS155*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS160*DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) - CLHS177*CLHS214 + CLHS178*CLHS219 - CLHS215*CLHS35 - CLHS216*CLHS37 + CLHS217*CLHS35 + CLHS218*CLHS37
    const double clhs231 =     clhs208*clhs44;
    const double clhs232 =     clhs210*clhs46;
    const double clhs233 =     clhs2*(clhs81 - clhs83 - clhs84);
    const double clhs234 =     clhs152*clhs233;
    const double clhs235 =     -clhs155*DeltaNormals[3](0,1) + clhs160*Deltatangentxis[3](0,1) - clhs184*clhs214 + clhs185*clhs219 - clhs215*clhs44 - clhs216*clhs46 + clhs217*clhs44 + clhs218*clhs46; // -CLHS155*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS160*DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) - CLHS184*CLHS214 + CLHS185*CLHS219 - CLHS215*CLHS44 - CLHS216*CLHS46 + CLHS217*CLHS44 + CLHS218*CLHS46
    const double clhs236 =     clhs188*clhs219;
    const double clhs237 =     clhs190*clhs202;
    const double clhs238 =     clhs194*clhs219;
    const double clhs239 =     clhs196*clhs202;
    const double clhs240 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs241 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs242 =     clhs55*(-clhs133*clhs240 + clhs136*clhs241);
    const double clhs243 =     clhs55*(-clhs138*clhs240 + clhs139*clhs241);
    const double clhs244 =     clhs55*(-clhs141*clhs240 + clhs142*clhs241);
    const double clhs245 =     clhs55*(-clhs144*clhs240 + clhs145*clhs241);
    const double clhs246 =     clhs147*clhs4*epsilon;
    const double clhs247 =     clhs13*clhs246;
    const double clhs248 =     clhs150*clhs20;
    const double clhs249 =     clhs153*clhs192;
    const double clhs250 =     clhs132*clhs2*clhs4;
    const double clhs251 =     clhs2*clhs240*clhs4;
    const double clhs252 =     clhs132*clhs240*clhs4;
    const double clhs253 =     clhs132*clhs2*clhs240;
    const double clhs254 =     clhs135*clhs2*clhs4;
    const double clhs255 =     clhs135*clhs241*clhs4;
    const double clhs256 =     clhs135*clhs2*clhs241;
    const double clhs257 =     clhs2*clhs241*clhs4;
    const double clhs258 =     -clhs13*clhs252 + clhs13*clhs255 - clhs156*clhs251 + clhs163*clhs257 - clhs20*clhs253 + clhs20*clhs256 - clhs250*DeltaNormals[0](1,0) + clhs254*Deltatangentxis[0](1,0); // -CLHS13*CLHS252 + CLHS13*CLHS255 - CLHS156*CLHS251 + CLHS163*CLHS257 - CLHS20*CLHS253 + CLHS20*CLHS256 - CLHS250*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS254*DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs259 =     clhs246*clhs26;
    const double clhs260 =     clhs150*clhs29;
    const double clhs261 =     clhs168*clhs192;
    const double clhs262 =     -clhs170*clhs251 + clhs171*clhs257 - clhs250*DeltaNormals[1](1,0) - clhs252*clhs26 - clhs253*clhs29 + clhs254*Deltatangentxis[1](1,0) + clhs255*clhs26 + clhs256*clhs29; // -CLHS170*CLHS251 + CLHS171*CLHS257 - CLHS250*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) - CLHS252*CLHS26 - CLHS253*CLHS29 + CLHS254*DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS255*CLHS26 + CLHS256*CLHS29
    const double clhs263 =     clhs246*clhs35;
    const double clhs264 =     clhs150*clhs38;
    const double clhs265 =     clhs175*clhs192;
    const double clhs266 =     -clhs177*clhs251 + clhs178*clhs257 - clhs250*DeltaNormals[2](1,0) - clhs252*clhs35 - clhs253*clhs38 + clhs254*Deltatangentxis[2](1,0) + clhs255*clhs35 + clhs256*clhs38; // -CLHS177*CLHS251 + CLHS178*CLHS257 - CLHS250*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) - CLHS252*CLHS35 - CLHS253*CLHS38 + CLHS254*DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS255*CLHS35 + CLHS256*CLHS38
    const double clhs267 =     clhs246*clhs44;
    const double clhs268 =     clhs150*clhs47;
    const double clhs269 =     clhs182*clhs192;
    const double clhs270 =     -clhs184*clhs251 + clhs185*clhs257 - clhs250*DeltaNormals[3](1,0) - clhs252*clhs44 - clhs253*clhs47 + clhs254*Deltatangentxis[3](1,0) + clhs255*clhs44 + clhs256*clhs47; // -CLHS184*CLHS251 + CLHS185*CLHS257 - CLHS250*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) - CLHS252*CLHS44 - CLHS253*CLHS47 + CLHS254*DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS255*CLHS44 + CLHS256*CLHS47
    const double clhs271 =     clhs188*clhs241;
    const double clhs272 =     clhs190*clhs257;
    const double clhs273 =     clhs194*clhs241;
    const double clhs274 =     clhs196*clhs257;
    const double clhs275 =     clhs198*std::pow(clhs4, 2);
    const double clhs276 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs277 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs278 =     clhs55*(-clhs133*clhs276 + clhs136*clhs277);
    const double clhs279 =     clhs55*(-clhs138*clhs276 + clhs139*clhs277);
    const double clhs280 =     clhs55*(-clhs141*clhs276 + clhs142*clhs277);
    const double clhs281 =     clhs55*(-clhs144*clhs276 + clhs145*clhs277);
    const double clhs282 =     clhs207*clhs4*epsilon;
    const double clhs283 =     clhs13*clhs282;
    const double clhs284 =     clhs20*clhs210;
    const double clhs285 =     clhs192*clhs212;
    const double clhs286 =     clhs2*clhs276*clhs4;
    const double clhs287 =     clhs132*clhs276*clhs4;
    const double clhs288 =     clhs132*clhs2*clhs276;
    const double clhs289 =     clhs135*clhs277*clhs4;
    const double clhs290 =     clhs135*clhs2*clhs277;
    const double clhs291 =     clhs2*clhs277*clhs4;
    const double clhs292 =     -clhs13*clhs287 + clhs13*clhs289 - clhs156*clhs286 + clhs163*clhs291 - clhs20*clhs288 + clhs20*clhs290 - clhs250*DeltaNormals[0](1,1) + clhs254*Deltatangentxis[0](1,1); // -CLHS13*CLHS287 + CLHS13*CLHS289 - CLHS156*CLHS286 + CLHS163*CLHS291 - CLHS20*CLHS288 + CLHS20*CLHS290 - CLHS250*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS254*DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs293 =     clhs26*clhs282;
    const double clhs294 =     clhs210*clhs29;
    const double clhs295 =     clhs192*clhs223;
    const double clhs296 =     -clhs170*clhs286 + clhs171*clhs291 - clhs250*DeltaNormals[1](1,1) + clhs254*Deltatangentxis[1](1,1) - clhs26*clhs287 + clhs26*clhs289 - clhs288*clhs29 + clhs29*clhs290; // -CLHS170*CLHS286 + CLHS171*CLHS291 - CLHS250*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS254*DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) - CLHS26*CLHS287 + CLHS26*CLHS289 - CLHS288*CLHS29 + CLHS29*CLHS290
    const double clhs297 =     clhs282*clhs35;
    const double clhs298 =     clhs210*clhs38;
    const double clhs299 =     clhs192*clhs228;
    const double clhs300 =     -clhs177*clhs286 + clhs178*clhs291 - clhs250*DeltaNormals[2](1,1) + clhs254*Deltatangentxis[2](1,1) - clhs287*clhs35 - clhs288*clhs38 + clhs289*clhs35 + clhs290*clhs38; // -CLHS177*CLHS286 + CLHS178*CLHS291 - CLHS250*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS254*DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) - CLHS287*CLHS35 - CLHS288*CLHS38 + CLHS289*CLHS35 + CLHS290*CLHS38
    const double clhs301 =     clhs282*clhs44;
    const double clhs302 =     clhs210*clhs47;
    const double clhs303 =     clhs192*clhs233;
    const double clhs304 =     -clhs184*clhs286 + clhs185*clhs291 - clhs250*DeltaNormals[3](1,1) + clhs254*Deltatangentxis[3](1,1) - clhs287*clhs44 - clhs288*clhs47 + clhs289*clhs44 + clhs290*clhs47; // -CLHS184*CLHS286 + CLHS185*CLHS291 - CLHS250*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS254*DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) - CLHS287*CLHS44 - CLHS288*CLHS47 + CLHS289*CLHS44 + CLHS290*CLHS47
    const double clhs305 =     clhs188*clhs291;
    const double clhs306 =     clhs190*clhs277;
    const double clhs307 =     clhs194*clhs291;
    const double clhs308 =     clhs196*clhs277;
    const double clhs309 =     -clhs199;
    const double clhs310 =     -clhs200;
    const double clhs311 =     -clhs275;

    lhs(0,0)=clhs1*clhs9;
    lhs(0,1)=clhs10*clhs9;
    lhs(0,2)=clhs11*clhs9;
    lhs(0,3)=clhs12*clhs9;
    lhs(0,4)=clhs0*clhs22 + clhs0*clhs25 + clhs13*clhs14 + clhs13*clhs17 + clhs15*clhs16 + clhs15*clhs18;
    lhs(0,5)=clhs0*clhs31 + clhs0*clhs34 + clhs14*clhs26 + clhs16*clhs27 + clhs17*clhs26 + clhs18*clhs27;
    lhs(0,6)=clhs0*clhs40 + clhs0*clhs43 + clhs14*clhs35 + clhs16*clhs36 + clhs17*clhs35 + clhs18*clhs36;
    lhs(0,7)=clhs0*clhs49 + clhs0*clhs52 + clhs14*clhs44 + clhs16*clhs45 + clhs17*clhs44 + clhs18*clhs45;
    lhs(0,8)=clhs54;
    lhs(0,9)=0;
    lhs(0,10)=clhs56;
    lhs(0,11)=0;
    lhs(0,12)=clhs54;
    lhs(0,13)=0;
    lhs(0,14)=clhs56;
    lhs(0,15)=0;
    lhs(1,0)=clhs1*clhs61;
    lhs(1,1)=clhs10*clhs61;
    lhs(1,2)=clhs11*clhs61;
    lhs(1,3)=clhs12*clhs61;
    lhs(1,4)=clhs0*clhs67 + clhs0*clhs70 + clhs13*clhs62 + clhs13*clhs64 + clhs15*clhs63 + clhs15*clhs65;
    lhs(1,5)=clhs0*clhs72 + clhs0*clhs75 + clhs26*clhs62 + clhs26*clhs64 + clhs27*clhs63 + clhs27*clhs65;
    lhs(1,6)=clhs0*clhs77 + clhs0*clhs80 + clhs35*clhs62 + clhs35*clhs64 + clhs36*clhs63 + clhs36*clhs65;
    lhs(1,7)=clhs0*clhs82 + clhs0*clhs85 + clhs44*clhs62 + clhs44*clhs64 + clhs45*clhs63 + clhs45*clhs65;
    lhs(1,8)=0;
    lhs(1,9)=clhs54;
    lhs(1,10)=0;
    lhs(1,11)=clhs56;
    lhs(1,12)=0;
    lhs(1,13)=clhs54;
    lhs(1,14)=0;
    lhs(1,15)=clhs56;
    lhs(2,0)=clhs87*clhs9;
    lhs(2,1)=clhs88*clhs9;
    lhs(2,2)=clhs89*clhs9;
    lhs(2,3)=clhs9*clhs90;
    lhs(2,4)=clhs13*clhs91 + clhs13*clhs93 + clhs16*clhs92 + clhs18*clhs92 + clhs22*clhs86 + clhs25*clhs86;
    lhs(2,5)=clhs16*clhs94 + clhs18*clhs94 + clhs26*clhs91 + clhs26*clhs93 + clhs31*clhs86 + clhs34*clhs86;
    lhs(2,6)=clhs16*clhs95 + clhs18*clhs95 + clhs35*clhs91 + clhs35*clhs93 + clhs40*clhs86 + clhs43*clhs86;
    lhs(2,7)=clhs16*clhs96 + clhs18*clhs96 + clhs44*clhs91 + clhs44*clhs93 + clhs49*clhs86 + clhs52*clhs86;
    lhs(2,8)=clhs97;
    lhs(2,9)=0;
    lhs(2,10)=clhs98;
    lhs(2,11)=0;
    lhs(2,12)=clhs97;
    lhs(2,13)=0;
    lhs(2,14)=clhs98;
    lhs(2,15)=0;
    lhs(3,0)=clhs61*clhs87;
    lhs(3,1)=clhs61*clhs88;
    lhs(3,2)=clhs61*clhs89;
    lhs(3,3)=clhs61*clhs90;
    lhs(3,4)=clhs100*clhs13 + clhs13*clhs99 + clhs63*clhs92 + clhs65*clhs92 + clhs67*clhs86 + clhs70*clhs86;
    lhs(3,5)=clhs100*clhs26 + clhs26*clhs99 + clhs63*clhs94 + clhs65*clhs94 + clhs72*clhs86 + clhs75*clhs86;
    lhs(3,6)=clhs100*clhs35 + clhs35*clhs99 + clhs63*clhs95 + clhs65*clhs95 + clhs77*clhs86 + clhs80*clhs86;
    lhs(3,7)=clhs100*clhs44 + clhs44*clhs99 + clhs63*clhs96 + clhs65*clhs96 + clhs82*clhs86 + clhs85*clhs86;
    lhs(3,8)=0;
    lhs(3,9)=clhs97;
    lhs(3,10)=0;
    lhs(3,11)=clhs98;
    lhs(3,12)=0;
    lhs(3,13)=clhs97;
    lhs(3,14)=0;
    lhs(3,15)=clhs98;
    lhs(4,0)=-clhs102*clhs9;
    lhs(4,1)=-clhs103*clhs9;
    lhs(4,2)=-clhs104*clhs9;
    lhs(4,3)=-clhs105*clhs9;
    lhs(4,4)=-clhs101*clhs22 - clhs101*clhs25 - clhs106*clhs13 - clhs107*clhs16 - clhs107*clhs18 - clhs108*clhs13;
    lhs(4,5)=-clhs101*clhs31 - clhs101*clhs34 - clhs106*clhs26 - clhs108*clhs26 - clhs109*clhs16 - clhs109*clhs18;
    lhs(4,6)=-clhs101*clhs40 - clhs101*clhs43 - clhs106*clhs35 - clhs108*clhs35 - clhs110*clhs16 - clhs110*clhs18;
    lhs(4,7)=-clhs101*clhs49 - clhs101*clhs52 - clhs106*clhs44 - clhs108*clhs44 - clhs111*clhs16 - clhs111*clhs18;
    lhs(4,8)=clhs112;
    lhs(4,9)=0;
    lhs(4,10)=clhs113;
    lhs(4,11)=0;
    lhs(4,12)=clhs112;
    lhs(4,13)=0;
    lhs(4,14)=clhs113;
    lhs(4,15)=0;
    lhs(5,0)=-clhs102*clhs61;
    lhs(5,1)=-clhs103*clhs61;
    lhs(5,2)=-clhs104*clhs61;
    lhs(5,3)=-clhs105*clhs61;
    lhs(5,4)=-clhs101*clhs67 - clhs101*clhs70 - clhs107*clhs63 - clhs107*clhs65 - clhs114*clhs13 - clhs115*clhs13;
    lhs(5,5)=-clhs101*clhs72 - clhs101*clhs75 - clhs109*clhs63 - clhs109*clhs65 - clhs114*clhs26 - clhs115*clhs26;
    lhs(5,6)=-clhs101*clhs77 - clhs101*clhs80 - clhs110*clhs63 - clhs110*clhs65 - clhs114*clhs35 - clhs115*clhs35;
    lhs(5,7)=-clhs101*clhs82 - clhs101*clhs85 - clhs111*clhs63 - clhs111*clhs65 - clhs114*clhs44 - clhs115*clhs44;
    lhs(5,8)=0;
    lhs(5,9)=clhs112;
    lhs(5,10)=0;
    lhs(5,11)=clhs113;
    lhs(5,12)=0;
    lhs(5,13)=clhs112;
    lhs(5,14)=0;
    lhs(5,15)=clhs113;
    lhs(6,0)=-clhs117*clhs9;
    lhs(6,1)=-clhs118*clhs9;
    lhs(6,2)=-clhs119*clhs9;
    lhs(6,3)=-clhs120*clhs9;
    lhs(6,4)=-clhs116*clhs22 - clhs116*clhs25 - clhs121*clhs13 - clhs122*clhs16 - clhs122*clhs18 - clhs123*clhs13;
    lhs(6,5)=-clhs116*clhs31 - clhs116*clhs34 - clhs121*clhs26 - clhs123*clhs26 - clhs124*clhs16 - clhs124*clhs18;
    lhs(6,6)=-clhs116*clhs40 - clhs116*clhs43 - clhs121*clhs35 - clhs123*clhs35 - clhs125*clhs16 - clhs125*clhs18;
    lhs(6,7)=-clhs116*clhs49 - clhs116*clhs52 - clhs121*clhs44 - clhs123*clhs44 - clhs126*clhs16 - clhs126*clhs18;
    lhs(6,8)=clhs127;
    lhs(6,9)=0;
    lhs(6,10)=clhs128;
    lhs(6,11)=0;
    lhs(6,12)=clhs127;
    lhs(6,13)=0;
    lhs(6,14)=clhs128;
    lhs(6,15)=0;
    lhs(7,0)=-clhs117*clhs61;
    lhs(7,1)=-clhs118*clhs61;
    lhs(7,2)=-clhs119*clhs61;
    lhs(7,3)=-clhs120*clhs61;
    lhs(7,4)=-clhs116*clhs67 - clhs116*clhs70 - clhs122*clhs63 - clhs122*clhs65 - clhs129*clhs13 - clhs13*clhs130;
    lhs(7,5)=-clhs116*clhs72 - clhs116*clhs75 - clhs124*clhs63 - clhs124*clhs65 - clhs129*clhs26 - clhs130*clhs26;
    lhs(7,6)=-clhs116*clhs77 - clhs116*clhs80 - clhs125*clhs63 - clhs125*clhs65 - clhs129*clhs35 - clhs130*clhs35;
    lhs(7,7)=-clhs116*clhs82 - clhs116*clhs85 - clhs126*clhs63 - clhs126*clhs65 - clhs129*clhs44 - clhs130*clhs44;
    lhs(7,8)=0;
    lhs(7,9)=clhs127;
    lhs(7,10)=0;
    lhs(7,11)=clhs128;
    lhs(7,12)=0;
    lhs(7,13)=clhs127;
    lhs(7,14)=0;
    lhs(7,15)=clhs128;
    lhs(8,0)=clhs137;
    lhs(8,1)=clhs140;
    lhs(8,2)=clhs143;
    lhs(8,3)=clhs146;
    lhs(8,4)=clhs149 + clhs151 + clhs154 + clhs165;
    lhs(8,5)=clhs166 + clhs167 + clhs169 + clhs172;
    lhs(8,6)=clhs173 + clhs174 + clhs176 + clhs179;
    lhs(8,7)=clhs180 + clhs181 + clhs183 + clhs186;
    lhs(8,8)=clhs53*(clhs187 + clhs189);
    lhs(8,9)=clhs191;
    lhs(8,10)=clhs53*(clhs193 + clhs195);
    lhs(8,11)=clhs197;
    lhs(8,12)=clhs199;
    lhs(8,13)=0;
    lhs(8,14)=clhs200;
    lhs(8,15)=0;
    lhs(9,0)=clhs203;
    lhs(9,1)=clhs204;
    lhs(9,2)=clhs205;
    lhs(9,3)=clhs206;
    lhs(9,4)=clhs209 + clhs211 + clhs213 + clhs220;
    lhs(9,5)=clhs221 + clhs222 + clhs224 + clhs225;
    lhs(9,6)=clhs226 + clhs227 + clhs229 + clhs230;
    lhs(9,7)=clhs231 + clhs232 + clhs234 + clhs235;
    lhs(9,8)=clhs236;
    lhs(9,9)=clhs53*(clhs187 + clhs237);
    lhs(9,10)=clhs238;
    lhs(9,11)=clhs53*(clhs193 + clhs239);
    lhs(9,12)=0;
    lhs(9,13)=clhs199;
    lhs(9,14)=0;
    lhs(9,15)=clhs200;
    lhs(10,0)=clhs242;
    lhs(10,1)=clhs243;
    lhs(10,2)=clhs244;
    lhs(10,3)=clhs245;
    lhs(10,4)=clhs247 + clhs248 + clhs249 + clhs258;
    lhs(10,5)=clhs259 + clhs260 + clhs261 + clhs262;
    lhs(10,6)=clhs263 + clhs264 + clhs265 + clhs266;
    lhs(10,7)=clhs267 + clhs268 + clhs269 + clhs270;
    lhs(10,8)=clhs55*(clhs187 + clhs271);
    lhs(10,9)=clhs272;
    lhs(10,10)=clhs55*(clhs193 + clhs273);
    lhs(10,11)=clhs274;
    lhs(10,12)=clhs200;
    lhs(10,13)=0;
    lhs(10,14)=clhs275;
    lhs(10,15)=0;
    lhs(11,0)=clhs278;
    lhs(11,1)=clhs279;
    lhs(11,2)=clhs280;
    lhs(11,3)=clhs281;
    lhs(11,4)=clhs283 + clhs284 + clhs285 + clhs292;
    lhs(11,5)=clhs293 + clhs294 + clhs295 + clhs296;
    lhs(11,6)=clhs297 + clhs298 + clhs299 + clhs300;
    lhs(11,7)=clhs301 + clhs302 + clhs303 + clhs304;
    lhs(11,8)=clhs305;
    lhs(11,9)=clhs55*(clhs187 + clhs306);
    lhs(11,10)=clhs307;
    lhs(11,11)=clhs55*(clhs193 + clhs308);
    lhs(11,12)=0;
    lhs(11,13)=clhs200;
    lhs(11,14)=0;
    lhs(11,15)=clhs275;
    lhs(12,0)=clhs137;
    lhs(12,1)=clhs140;
    lhs(12,2)=clhs143;
    lhs(12,3)=clhs146;
    lhs(12,4)=-clhs149 - clhs151 - clhs154 + clhs165;
    lhs(12,5)=-clhs166 - clhs167 - clhs169 + clhs172;
    lhs(12,6)=-clhs173 - clhs174 - clhs176 + clhs179;
    lhs(12,7)=-clhs180 - clhs181 - clhs183 + clhs186;
    lhs(12,8)=clhs53*(clhs152 + clhs189);
    lhs(12,9)=clhs191;
    lhs(12,10)=clhs53*(clhs192 + clhs195);
    lhs(12,11)=clhs197;
    lhs(12,12)=clhs309;
    lhs(12,13)=0;
    lhs(12,14)=clhs310;
    lhs(12,15)=0;
    lhs(13,0)=clhs203;
    lhs(13,1)=clhs204;
    lhs(13,2)=clhs205;
    lhs(13,3)=clhs206;
    lhs(13,4)=-clhs209 - clhs211 - clhs213 + clhs220;
    lhs(13,5)=-clhs221 - clhs222 - clhs224 + clhs225;
    lhs(13,6)=-clhs226 - clhs227 - clhs229 + clhs230;
    lhs(13,7)=-clhs231 - clhs232 - clhs234 + clhs235;
    lhs(13,8)=clhs236;
    lhs(13,9)=clhs53*(clhs152 + clhs237);
    lhs(13,10)=clhs238;
    lhs(13,11)=clhs53*(clhs192 + clhs239);
    lhs(13,12)=0;
    lhs(13,13)=clhs309;
    lhs(13,14)=0;
    lhs(13,15)=clhs310;
    lhs(14,0)=clhs242;
    lhs(14,1)=clhs243;
    lhs(14,2)=clhs244;
    lhs(14,3)=clhs245;
    lhs(14,4)=-clhs247 - clhs248 - clhs249 + clhs258;
    lhs(14,5)=-clhs259 - clhs260 - clhs261 + clhs262;
    lhs(14,6)=-clhs263 - clhs264 - clhs265 + clhs266;
    lhs(14,7)=-clhs267 - clhs268 - clhs269 + clhs270;
    lhs(14,8)=clhs55*(clhs152 + clhs271);
    lhs(14,9)=clhs272;
    lhs(14,10)=clhs55*(clhs192 + clhs273);
    lhs(14,11)=clhs274;
    lhs(14,12)=clhs310;
    lhs(14,13)=0;
    lhs(14,14)=clhs311;
    lhs(14,15)=0;
    lhs(15,0)=clhs278;
    lhs(15,1)=clhs279;
    lhs(15,2)=clhs280;
    lhs(15,3)=clhs281;
    lhs(15,4)=-clhs283 - clhs284 - clhs285 + clhs292;
    lhs(15,5)=-clhs293 - clhs294 - clhs295 + clhs296;
    lhs(15,6)=-clhs297 - clhs298 - clhs299 + clhs300;
    lhs(15,7)=-clhs301 - clhs302 - clhs303 + clhs304;
    lhs(15,8)=clhs305;
    lhs(15,9)=clhs55*(clhs152 + clhs306);
    lhs(15,10)=clhs307;
    lhs(15,11)=clhs55*(clhs192 + clhs308);
    lhs(15,12)=0;
    lhs(15,13)=clhs310;
    lhs(15,14)=0;
    lhs(15,15)=clhs311;

    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,16> ComputeGaussPointFrictionlessRHS(
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const ContactData<2,2>& rContactData
        )
{
    array_1d<double,16> rhs;
    
    const double epsilon = rContactData.epsilon;
    const double integration_point_gap = inner_prod(rContactData.Gaps,N1);
    const bounded_matrix<double, 2, 2> normalslave = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave   = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm          = rContactData.LagrangeMultipliers;
    const bounded_matrix<double, 2, 2> dlm         = rContactData.DoubleLagrangeMultipliers;
    
    const double Dt = rContactData.Dt;
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
    const double crhs16 =     epsilon*(crhs4 - crhs5 - crhs6);
    const double crhs17 =     -crhs16;
    const double crhs18 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs19 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs20 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs21 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs22 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs23 =     ((crhs13*crhs19 + crhs14*crhs20)*(-crhs0*(Dt*v2(0,0)) - crhs12*(Dt*v2(1,0)) + crhs13*(Dt*v1(0,0)) + crhs14*(Dt*v1(1,0))) + (crhs13*crhs21 + crhs14*crhs22)*(-crhs0*(Dt*v2(0,1)) - crhs12*(Dt*v2(1,1)) + crhs13*(Dt*v1(0,1)) + crhs14*(Dt*v1(1,1))))/Dt;
    const double crhs24 =     crhs18*normalslave(0,0) - crhs19*crhs23; // CRHS18*NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CRHS19*CRHS23
    const double crhs25 =     epsilon*(-crhs10 + crhs8 - crhs9);
    const double crhs26 =     -crhs25;
    const double crhs27 =     crhs18*normalslave(0,1) - crhs21*crhs23; // CRHS18*NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CRHS21*CRHS23
    const double crhs28 =     crhs1*crhs3;
    const double crhs29 =     crhs18*normalslave(1,0) - crhs20*crhs23; // CRHS18*NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CRHS20*CRHS23
    const double crhs30 =     crhs18*normalslave(1,1) - crhs22*crhs23; // CRHS18*NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CRHS22*CRHS23

    rhs[0]=-crhs0*crhs7;
    rhs[1]=-crhs0*crhs11;
    rhs[2]=-crhs12*crhs7;
    rhs[3]=-crhs11*crhs12;
    rhs[4]=crhs13*crhs7;
    rhs[5]=crhs11*crhs13;
    rhs[6]=crhs14*crhs7;
    rhs[7]=crhs11*crhs14;
    rhs[8]=crhs15*(crhs17 + crhs24);
    rhs[9]=crhs15*(crhs26 + crhs27);
    rhs[10]=crhs28*(crhs17 + crhs29);
    rhs[11]=crhs28*(crhs26 + crhs30);
    rhs[12]=crhs15*(crhs16 + crhs24);
    rhs[13]=crhs15*(crhs25 + crhs27);
    rhs[14]=crhs28*(crhs16 + crhs29);
    rhs[15]=crhs28*(crhs25 + crhs30);

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,16> ComputeGaussPointRHS(
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const double mu, 
        const ContactData<2,2>& rContactData
        )
{
    array_1d<double,16> rhs;
    
    const double epsilon = rContactData.epsilon;
    const double integration_point_gap = inner_prod(rContactData.Gaps,N1);
    const bounded_matrix<double, 2, 2> normalslave = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave   = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm          = rContactData.LagrangeMultipliers;    
    const bounded_matrix<double, 2, 2> dlm         = rContactData.DoubleLagrangeMultipliers;    
    
    const array_1d<double, 1> Ctan = rContactData.Ctan; 
    
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
    const double crhs16 =     epsilon*(crhs4 - crhs5 - crhs6);
    const double crhs17 =     -crhs16;
    const double crhs18 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs19 =     Ctan[0]; // CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1))
    const double crhs20 =     crhs18*normalslave(0,0) - crhs19*tan1slave(0,0); // CRHS18*NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CRHS19*TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs21 =     epsilon*(-crhs10 + crhs8 - crhs9);
    const double crhs22 =     -crhs21;
    const double crhs23 =     crhs18*normalslave(0,1) - crhs19*tan1slave(0,1); // CRHS18*NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CRHS19*TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs24 =     crhs1*crhs3;
    const double crhs25 =     crhs18*normalslave(1,0) - crhs19*tan1slave(1,0); // CRHS18*NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CRHS19*TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs26 =     crhs18*normalslave(1,1) - crhs19*tan1slave(1,1); // CRHS18*NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CRHS19*TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))

    rhs[0]=-crhs0*crhs7;
    rhs[1]=-crhs0*crhs11;
    rhs[2]=-crhs12*crhs7;
    rhs[3]=-crhs11*crhs12;
    rhs[4]=crhs13*crhs7;
    rhs[5]=crhs11*crhs13;
    rhs[6]=crhs14*crhs7;
    rhs[7]=crhs11*crhs14;
    rhs[8]=crhs15*(crhs17 + crhs20);
    rhs[9]=crhs15*(crhs22 + crhs23);
    rhs[10]=crhs24*(crhs17 + crhs25);
    rhs[11]=crhs24*(crhs22 + crhs26);
    rhs[12]=crhs15*(crhs16 + crhs20);
    rhs[13]=crhs15*(crhs21 + crhs23);
    rhs[14]=crhs24*(crhs16 + crhs25);
    rhs[15]=crhs24*(crhs21 + crhs26);

    
    return rhs;
}

private:
};// class Contact2D2N2NDLM
}
#endif /* KRATOS_CONTACT2D2N2NDLM defined */
