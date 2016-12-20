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

#if !defined(KRATOS_CONTACT2D2N2N)
#define KRATOS_CONTACT2D2N2N

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
        
class Contact2D2N2N
{
public:
    
    static inline bounded_matrix<double,12,12> ComputeGaussPointActiveLHS(
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
    bounded_matrix<double,12,12> lhs;
    
    const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
//     const double Dt              = rContactData.Dt;
//     const double epsilon_normal  = rContactData.epsilon_normal;
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
    const double clhs4 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs5 =     clhs3*lm(0,0) + clhs4*lm(1,0);
    const double clhs6 =     clhs2*clhs5;
    const double clhs7 =     DeltaN2[5][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs8 =     DeltaN2[6][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs9 =     DeltaN2[7][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs10 =     DeltaJs[0]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs11 =     clhs0*clhs5;
    const double clhs12 =     DeltaN2[0][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs13 =     clhs0*clhs2;
    const double clhs14 =     DeltaPhi[0][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs15 =     DeltaPhi[0][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs16 =     clhs14*lm(0,0) + clhs15*lm(1,0);
    const double clhs17 =     DeltaJs[1]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs18 =     DeltaN2[1][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs19 =     DeltaPhi[1][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs20 =     DeltaPhi[1][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs21 =     clhs19*lm(0,0) + clhs20*lm(1,0);
    const double clhs22 =     DeltaJs[2]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs23 =     DeltaN2[2][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs24 =     DeltaPhi[2][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs25 =     DeltaPhi[2][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs26 =     clhs24*lm(0,0) + clhs25*lm(1,0);
    const double clhs27 =     DeltaJs[3]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs28 =     DeltaN2[3][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs29 =     DeltaPhi[3][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs30 =     DeltaPhi[3][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs31 =     clhs29*lm(0,0) + clhs30*lm(1,0);
    const double clhs32 =     clhs13*clhs3;
    const double clhs33 =     clhs13*clhs4;
    const double clhs34 =     clhs3*lm(0,1) + clhs4*lm(1,1);
    const double clhs35 =     clhs2*clhs34;
    const double clhs36 =     clhs0*clhs34;
    const double clhs37 =     clhs14*lm(0,1) + clhs15*lm(1,1);
    const double clhs38 =     clhs19*lm(0,1) + clhs20*lm(1,1);
    const double clhs39 =     clhs24*lm(0,1) + clhs25*lm(1,1);
    const double clhs40 =     clhs29*lm(0,1) + clhs30*lm(1,1);
    const double clhs41 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs42 =     DeltaN2[4][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs43 =     DeltaN2[5][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs44 =     DeltaN2[6][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs45 =     DeltaN2[7][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs46 =     clhs41*clhs5;
    const double clhs47 =     DeltaN2[0][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs48 =     clhs2*clhs41;
    const double clhs49 =     DeltaN2[1][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs50 =     DeltaN2[2][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs51 =     DeltaN2[3][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs52 =     clhs3*clhs48;
    const double clhs53 =     clhs4*clhs48;
    const double clhs54 =     clhs34*clhs41;
    const double clhs55 =     N1[0]; // N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs56 =     DeltaN1[4][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs57 =     DeltaN1[5][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs58 =     DeltaN1[6][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs59 =     DeltaN1[7][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs60 =     clhs5*clhs55;
    const double clhs61 =     DeltaN1[0][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs62 =     clhs2*clhs55;
    const double clhs63 =     DeltaN1[1][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs64 =     DeltaN1[2][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs65 =     DeltaN1[3][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs66 =     -clhs3*clhs62;
    const double clhs67 =     -clhs4*clhs62;
    const double clhs68 =     clhs34*clhs55;
    const double clhs69 =     N1[1]; // N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs70 =     DeltaN1[4][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs71 =     DeltaN1[5][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs72 =     DeltaN1[6][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs73 =     DeltaN1[7][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs74 =     clhs5*clhs69;
    const double clhs75 =     DeltaN1[0][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs76 =     clhs2*clhs69;
    const double clhs77 =     DeltaN1[1][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs78 =     DeltaN1[2][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs79 =     DeltaN1[3][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs80 =     -clhs3*clhs76;
    const double clhs81 =     -clhs4*clhs76;
    const double clhs82 =     clhs34*clhs69;
    const double clhs83 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs84 =     DeltaGap[4]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs85 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs86 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs87 =     clhs2*clhs3*(GPnormal[0]*clhs85 + GPnormal[1]*clhs86);
    const double clhs88 =     DeltaGap[5]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs89 =     DeltaGap[6]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs90 =     DeltaGap[7]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs91 =     clhs2*clhs3*clhs83;
    const double clhs92 =     DeltaGap[0]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs93 =     clhs2*clhs3*clhs85;
    const double clhs94 =     clhs3*clhs83*clhs85;
    const double clhs95 =     clhs2*clhs83*clhs85;
    const double clhs96 =     clhs10*clhs94 + clhs14*clhs95 + clhs91*DeltaNormals[0](0,0) + clhs92*clhs93; // CLHS10*CLHS94 + CLHS14*CLHS95 + CLHS91*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS92*CLHS93
    const double clhs97 =     clhs2*clhs3*clhs86;
    const double clhs98 =     clhs3*clhs83*clhs86;
    const double clhs99 =     clhs2*clhs83*clhs86;
    const double clhs100 =     clhs10*clhs98 + clhs14*clhs99 + clhs91*DeltaNormals[0](0,1) + clhs92*clhs97; // CLHS10*CLHS98 + CLHS14*CLHS99 + CLHS91*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS92*CLHS97
    const double clhs101 =     DeltaGap[1]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs102 =     clhs101*clhs93 + clhs17*clhs94 + clhs19*clhs95 + clhs91*DeltaNormals[1](0,0); // CLHS101*CLHS93 + CLHS17*CLHS94 + CLHS19*CLHS95 + CLHS91*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs103 =     clhs101*clhs97 + clhs17*clhs98 + clhs19*clhs99 + clhs91*DeltaNormals[1](0,1); // CLHS101*CLHS97 + CLHS17*CLHS98 + CLHS19*CLHS99 + CLHS91*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs104 =     DeltaGap[2]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs105 =     clhs104*clhs93 + clhs22*clhs94 + clhs24*clhs95 + clhs91*DeltaNormals[2](0,0); // CLHS104*CLHS93 + CLHS22*CLHS94 + CLHS24*CLHS95 + CLHS91*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs106 =     clhs104*clhs97 + clhs22*clhs98 + clhs24*clhs99 + clhs91*DeltaNormals[2](0,1); // CLHS104*CLHS97 + CLHS22*CLHS98 + CLHS24*CLHS99 + CLHS91*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs107 =     DeltaGap[3]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs108 =     clhs107*clhs93 + clhs27*clhs94 + clhs29*clhs95 + clhs91*DeltaNormals[3](0,0); // CLHS107*CLHS93 + CLHS27*CLHS94 + CLHS29*CLHS95 + CLHS91*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs109 =     clhs107*clhs97 + clhs27*clhs98 + clhs29*clhs99 + clhs91*DeltaNormals[3](0,1); // CLHS107*CLHS97 + CLHS27*CLHS98 + CLHS29*CLHS99 + CLHS91*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs110 =     clhs2*clhs3*(GPtangent1[0]*clhs85 + GPtangent1[1]*clhs86);
    const double clhs111 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs112 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs113 =     clhs2*clhs4*(GPnormal[0]*clhs111 + GPnormal[1]*clhs112);
    const double clhs114 =     clhs2*clhs4*clhs83;
    const double clhs115 =     clhs111*clhs2*clhs4;
    const double clhs116 =     clhs111*clhs4*clhs83;
    const double clhs117 =     clhs111*clhs2*clhs83;
    const double clhs118 =     clhs10*clhs116 + clhs114*DeltaNormals[0](1,0) + clhs115*clhs92 + clhs117*clhs15; // CLHS10*CLHS116 + CLHS114*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS115*CLHS92 + CLHS117*CLHS15
    const double clhs119 =     clhs112*clhs2*clhs4;
    const double clhs120 =     clhs112*clhs4*clhs83;
    const double clhs121 =     clhs112*clhs2*clhs83;
    const double clhs122 =     clhs10*clhs120 + clhs114*DeltaNormals[0](1,1) + clhs119*clhs92 + clhs121*clhs15; // CLHS10*CLHS120 + CLHS114*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS119*CLHS92 + CLHS121*CLHS15
    const double clhs123 =     clhs101*clhs115 + clhs114*DeltaNormals[1](1,0) + clhs116*clhs17 + clhs117*clhs20; // CLHS101*CLHS115 + CLHS114*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS116*CLHS17 + CLHS117*CLHS20
    const double clhs124 =     clhs101*clhs119 + clhs114*DeltaNormals[1](1,1) + clhs120*clhs17 + clhs121*clhs20; // CLHS101*CLHS119 + CLHS114*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS120*CLHS17 + CLHS121*CLHS20
    const double clhs125 =     clhs104*clhs115 + clhs114*DeltaNormals[2](1,0) + clhs116*clhs22 + clhs117*clhs25; // CLHS104*CLHS115 + CLHS114*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS116*CLHS22 + CLHS117*CLHS25
    const double clhs126 =     clhs104*clhs119 + clhs114*DeltaNormals[2](1,1) + clhs120*clhs22 + clhs121*clhs25; // CLHS104*CLHS119 + CLHS114*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS120*CLHS22 + CLHS121*CLHS25
    const double clhs127 =     clhs107*clhs115 + clhs114*DeltaNormals[3](1,0) + clhs116*clhs27 + clhs117*clhs30; // CLHS107*CLHS115 + CLHS114*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS116*CLHS27 + CLHS117*CLHS30
    const double clhs128 =     clhs107*clhs119 + clhs114*DeltaNormals[3](1,1) + clhs120*clhs27 + clhs121*clhs30; // CLHS107*CLHS119 + CLHS114*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS120*CLHS27 + CLHS121*CLHS30
    const double clhs129 =     clhs2*clhs4*(GPtangent1[0]*clhs111 + GPtangent1[1]*clhs112);

    lhs(0,0)=clhs1*clhs6;
    lhs(0,1)=clhs6*clhs7;
    lhs(0,2)=clhs6*clhs8;
    lhs(0,3)=clhs6*clhs9;
    lhs(0,4)=clhs10*clhs11 + clhs12*clhs6 + clhs13*clhs16;
    lhs(0,5)=clhs11*clhs17 + clhs13*clhs21 + clhs18*clhs6;
    lhs(0,6)=clhs11*clhs22 + clhs13*clhs26 + clhs23*clhs6;
    lhs(0,7)=clhs11*clhs27 + clhs13*clhs31 + clhs28*clhs6;
    lhs(0,8)=clhs32;
    lhs(0,9)=0;
    lhs(0,10)=clhs33;
    lhs(0,11)=0;
    lhs(1,0)=clhs1*clhs35;
    lhs(1,1)=clhs35*clhs7;
    lhs(1,2)=clhs35*clhs8;
    lhs(1,3)=clhs35*clhs9;
    lhs(1,4)=clhs10*clhs36 + clhs12*clhs35 + clhs13*clhs37;
    lhs(1,5)=clhs13*clhs38 + clhs17*clhs36 + clhs18*clhs35;
    lhs(1,6)=clhs13*clhs39 + clhs22*clhs36 + clhs23*clhs35;
    lhs(1,7)=clhs13*clhs40 + clhs27*clhs36 + clhs28*clhs35;
    lhs(1,8)=0;
    lhs(1,9)=clhs32;
    lhs(1,10)=0;
    lhs(1,11)=clhs33;
    lhs(2,0)=clhs42*clhs6;
    lhs(2,1)=clhs43*clhs6;
    lhs(2,2)=clhs44*clhs6;
    lhs(2,3)=clhs45*clhs6;
    lhs(2,4)=clhs10*clhs46 + clhs16*clhs48 + clhs47*clhs6;
    lhs(2,5)=clhs17*clhs46 + clhs21*clhs48 + clhs49*clhs6;
    lhs(2,6)=clhs22*clhs46 + clhs26*clhs48 + clhs50*clhs6;
    lhs(2,7)=clhs27*clhs46 + clhs31*clhs48 + clhs51*clhs6;
    lhs(2,8)=clhs52;
    lhs(2,9)=0;
    lhs(2,10)=clhs53;
    lhs(2,11)=0;
    lhs(3,0)=clhs35*clhs42;
    lhs(3,1)=clhs35*clhs43;
    lhs(3,2)=clhs35*clhs44;
    lhs(3,3)=clhs35*clhs45;
    lhs(3,4)=clhs10*clhs54 + clhs35*clhs47 + clhs37*clhs48;
    lhs(3,5)=clhs17*clhs54 + clhs35*clhs49 + clhs38*clhs48;
    lhs(3,6)=clhs22*clhs54 + clhs35*clhs50 + clhs39*clhs48;
    lhs(3,7)=clhs27*clhs54 + clhs35*clhs51 + clhs40*clhs48;
    lhs(3,8)=0;
    lhs(3,9)=clhs52;
    lhs(3,10)=0;
    lhs(3,11)=clhs53;
    lhs(4,0)=-clhs56*clhs6;
    lhs(4,1)=-clhs57*clhs6;
    lhs(4,2)=-clhs58*clhs6;
    lhs(4,3)=-clhs59*clhs6;
    lhs(4,4)=-clhs10*clhs60 - clhs16*clhs62 - clhs6*clhs61;
    lhs(4,5)=-clhs17*clhs60 - clhs21*clhs62 - clhs6*clhs63;
    lhs(4,6)=-clhs22*clhs60 - clhs26*clhs62 - clhs6*clhs64;
    lhs(4,7)=-clhs27*clhs60 - clhs31*clhs62 - clhs6*clhs65;
    lhs(4,8)=clhs66;
    lhs(4,9)=0;
    lhs(4,10)=clhs67;
    lhs(4,11)=0;
    lhs(5,0)=-clhs35*clhs56;
    lhs(5,1)=-clhs35*clhs57;
    lhs(5,2)=-clhs35*clhs58;
    lhs(5,3)=-clhs35*clhs59;
    lhs(5,4)=-clhs10*clhs68 - clhs35*clhs61 - clhs37*clhs62;
    lhs(5,5)=-clhs17*clhs68 - clhs35*clhs63 - clhs38*clhs62;
    lhs(5,6)=-clhs22*clhs68 - clhs35*clhs64 - clhs39*clhs62;
    lhs(5,7)=-clhs27*clhs68 - clhs35*clhs65 - clhs40*clhs62;
    lhs(5,8)=0;
    lhs(5,9)=clhs66;
    lhs(5,10)=0;
    lhs(5,11)=clhs67;
    lhs(6,0)=-clhs6*clhs70;
    lhs(6,1)=-clhs6*clhs71;
    lhs(6,2)=-clhs6*clhs72;
    lhs(6,3)=-clhs6*clhs73;
    lhs(6,4)=-clhs10*clhs74 - clhs16*clhs76 - clhs6*clhs75;
    lhs(6,5)=-clhs17*clhs74 - clhs21*clhs76 - clhs6*clhs77;
    lhs(6,6)=-clhs22*clhs74 - clhs26*clhs76 - clhs6*clhs78;
    lhs(6,7)=-clhs27*clhs74 - clhs31*clhs76 - clhs6*clhs79;
    lhs(6,8)=clhs80;
    lhs(6,9)=0;
    lhs(6,10)=clhs81;
    lhs(6,11)=0;
    lhs(7,0)=-clhs35*clhs70;
    lhs(7,1)=-clhs35*clhs71;
    lhs(7,2)=-clhs35*clhs72;
    lhs(7,3)=-clhs35*clhs73;
    lhs(7,4)=-clhs10*clhs82 - clhs35*clhs75 - clhs37*clhs76;
    lhs(7,5)=-clhs17*clhs82 - clhs35*clhs77 - clhs38*clhs76;
    lhs(7,6)=-clhs22*clhs82 - clhs35*clhs78 - clhs39*clhs76;
    lhs(7,7)=-clhs27*clhs82 - clhs35*clhs79 - clhs40*clhs76;
    lhs(7,8)=0;
    lhs(7,9)=clhs80;
    lhs(7,10)=0;
    lhs(7,11)=clhs81;
    lhs(8,0)=-clhs84*clhs87;
    lhs(8,1)=-clhs87*clhs88;
    lhs(8,2)=-clhs87*clhs89;
    lhs(8,3)=-clhs87*clhs90;
    lhs(8,4)=-GPnormal[0]*clhs96 - GPnormal[1]*clhs100;
    lhs(8,5)=-GPnormal[0]*clhs102 - GPnormal[1]*clhs103;
    lhs(8,6)=-GPnormal[0]*clhs105 - GPnormal[1]*clhs106;
    lhs(8,7)=-GPnormal[0]*clhs108 - GPnormal[1]*clhs109;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=-clhs110*clhs84;
    lhs(9,1)=-clhs110*clhs88;
    lhs(9,2)=-clhs110*clhs89;
    lhs(9,3)=-clhs110*clhs90;
    lhs(9,4)=-GPtangent1[0]*clhs96 - GPtangent1[1]*clhs100;
    lhs(9,5)=-GPtangent1[0]*clhs102 - GPtangent1[1]*clhs103;
    lhs(9,6)=-GPtangent1[0]*clhs105 - GPtangent1[1]*clhs106;
    lhs(9,7)=-GPtangent1[0]*clhs108 - GPtangent1[1]*clhs109;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=-clhs113*clhs84;
    lhs(10,1)=-clhs113*clhs88;
    lhs(10,2)=-clhs113*clhs89;
    lhs(10,3)=-clhs113*clhs90;
    lhs(10,4)=-GPnormal[0]*clhs118 - GPnormal[1]*clhs122;
    lhs(10,5)=-GPnormal[0]*clhs123 - GPnormal[1]*clhs124;
    lhs(10,6)=-GPnormal[0]*clhs125 - GPnormal[1]*clhs126;
    lhs(10,7)=-GPnormal[0]*clhs127 - GPnormal[1]*clhs128;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=-clhs129*clhs84;
    lhs(11,1)=-clhs129*clhs88;
    lhs(11,2)=-clhs129*clhs89;
    lhs(11,3)=-clhs129*clhs90;
    lhs(11,4)=-GPtangent1[0]*clhs118 - GPtangent1[1]*clhs122;
    lhs(11,5)=-GPtangent1[0]*clhs123 - GPtangent1[1]*clhs124;
    lhs(11,6)=-GPtangent1[0]*clhs125 - GPtangent1[1]*clhs126;
    lhs(11,7)=-GPtangent1[0]*clhs127 - GPtangent1[1]*clhs128;
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;

    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,12,12> ComputeGaussPointStickLHS(
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
    bounded_matrix<double,12,12> lhs;
    
    const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
    const double Dt              = rContactData.Dt;
//     const double epsilon_normal  = rContactData.epsilon_normal;
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
    const double clhs29 =     DeltaN1[5][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs30 =     DeltaN1[5][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs31 =     DeltaN2[5][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs32 =     DeltaN2[5][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs33 =     clhs17*(clhs29*clhs3 + clhs30*clhs8) + clhs23*(clhs18*clhs30 + clhs29*clhs4) - clhs24*(clhs14 - clhs19*clhs29 - clhs20*clhs30 + clhs21*clhs31 + clhs22*clhs32) + clhs27*(clhs11*clhs29 + clhs12*clhs30 - clhs13*clhs31 - clhs15*clhs32);
    const double clhs34 =     DeltaN1[6][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs35 =     DeltaN1[6][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs36 =     DeltaN2[6][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs37 =     DeltaN2[6][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs38 =     clhs17*(clhs3*clhs34 + clhs35*clhs8) + clhs23*(clhs18*clhs35 + clhs34*clhs4) + clhs24*(clhs19*clhs34 + clhs20*clhs35 - clhs21*clhs36 - clhs22*clhs37) - clhs27*(-clhs11*clhs34 - clhs12*clhs35 + clhs13*clhs36 + clhs15*clhs37 + clhs16);
    const double clhs39 =     DeltaN1[7][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs40 =     DeltaN1[7][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs41 =     DeltaN2[7][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs42 =     DeltaN2[7][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs43 =     clhs17*(clhs3*clhs39 + clhs40*clhs8) + clhs23*(clhs18*clhs40 + clhs39*clhs4) - clhs24*(clhs16 - clhs19*clhs39 - clhs20*clhs40 + clhs21*clhs41 + clhs22*clhs42) + clhs27*(clhs11*clhs39 + clhs12*clhs40 - clhs13*clhs41 - clhs15*clhs42);
    const double clhs44 =     Deltatangentxis[0](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs45 =     clhs17*clhs27 + clhs23*clhs24;
    const double clhs46 =     clhs1*clhs2*clhs45;
    const double clhs47 =     DeltaJs[0]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs48 =     clhs1*clhs3*clhs45;
    const double clhs49 =     DeltaPhi[0][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs50 =     clhs2*clhs3*clhs45;
    const double clhs51 =     DeltaN1[0][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs52 =     DeltaN1[0][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs53 =     DeltaN2[0][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs54 =     DeltaN2[0][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs55 =     Deltatangentxis[0](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs56 =     Deltatangentxis[0](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs57 =     Deltatangentxis[0](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs58 =     clhs17*(clhs3*clhs51 + clhs44*clhs6 + clhs52*clhs8 + clhs55*clhs9) + clhs23*(clhs18*clhs52 + clhs4*clhs51 + clhs56*clhs6 + clhs57*clhs9) + clhs24*(clhs19*clhs51 + clhs20*clhs52 - clhs21*clhs53 - clhs22*clhs54) + clhs27*(clhs11*clhs51 + clhs12*clhs52 - clhs13*clhs53 - clhs15*clhs54 + clhs6);
    const double clhs59 =     clhs1*clhs2*clhs58;
    const double clhs60 =     clhs3*clhs59 + clhs44*clhs46 + clhs47*clhs48 + clhs49*clhs50;
    const double clhs61 =     clhs1*clhs4*clhs45;
    const double clhs62 =     clhs2*clhs4*clhs45;
    const double clhs63 =     clhs4*clhs59 + clhs46*clhs56 + clhs47*clhs61 + clhs49*clhs62;
    const double clhs64 =     Deltatangentxis[1](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs65 =     DeltaJs[1]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs66 =     DeltaPhi[1][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs67 =     DeltaN1[1][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs68 =     DeltaN1[1][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs69 =     DeltaN2[1][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs70 =     DeltaN2[1][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs71 =     Deltatangentxis[1](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs72 =     Deltatangentxis[1](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs73 =     Deltatangentxis[1](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs74 =     clhs17*(clhs3*clhs67 + clhs6*clhs64 + clhs68*clhs8 + clhs71*clhs9) + clhs23*(clhs18*clhs68 + clhs4*clhs67 + clhs6*clhs72 + clhs73*clhs9) + clhs24*(clhs19*clhs67 + clhs20*clhs68 - clhs21*clhs69 - clhs22*clhs70 + clhs6) + clhs27*(clhs11*clhs67 + clhs12*clhs68 - clhs13*clhs69 - clhs15*clhs70);
    const double clhs75 =     clhs1*clhs2*clhs74;
    const double clhs76 =     clhs3*clhs75 + clhs46*clhs64 + clhs48*clhs65 + clhs50*clhs66;
    const double clhs77 =     clhs4*clhs75 + clhs46*clhs72 + clhs61*clhs65 + clhs62*clhs66;
    const double clhs78 =     Deltatangentxis[2](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs79 =     DeltaJs[2]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs80 =     DeltaPhi[2][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs81 =     DeltaN1[2][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs82 =     DeltaN1[2][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs83 =     DeltaN2[2][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs84 =     DeltaN2[2][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs85 =     Deltatangentxis[2](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs86 =     Deltatangentxis[2](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs87 =     Deltatangentxis[2](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs88 =     clhs17*(clhs3*clhs81 + clhs6*clhs78 + clhs8*clhs82 + clhs85*clhs9) + clhs23*(clhs18*clhs82 + clhs4*clhs81 + clhs6*clhs86 + clhs87*clhs9) + clhs24*(clhs19*clhs81 + clhs20*clhs82 - clhs21*clhs83 - clhs22*clhs84) + clhs27*(clhs11*clhs81 + clhs12*clhs82 - clhs13*clhs83 - clhs15*clhs84 + clhs9);
    const double clhs89 =     clhs1*clhs2*clhs88;
    const double clhs90 =     clhs3*clhs89 + clhs46*clhs78 + clhs48*clhs79 + clhs50*clhs80;
    const double clhs91 =     clhs4*clhs89 + clhs46*clhs86 + clhs61*clhs79 + clhs62*clhs80;
    const double clhs92 =     Deltatangentxis[3](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs93 =     DeltaJs[3]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs94 =     DeltaPhi[3][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs95 =     DeltaN1[3][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs96 =     DeltaN1[3][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs97 =     DeltaN2[3][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs98 =     DeltaN2[3][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs99 =     Deltatangentxis[3](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs100 =     Deltatangentxis[3](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs101 =     Deltatangentxis[3](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs102 =     clhs17*(clhs3*clhs95 + clhs6*clhs92 + clhs8*clhs96 + clhs9*clhs99) + clhs23*(clhs100*clhs6 + clhs101*clhs9 + clhs18*clhs96 + clhs4*clhs95) + clhs24*(clhs19*clhs95 + clhs20*clhs96 - clhs21*clhs97 - clhs22*clhs98 + clhs9) + clhs27*(clhs11*clhs95 + clhs12*clhs96 - clhs13*clhs97 - clhs15*clhs98);
    const double clhs103 =     clhs1*clhs102*clhs2;
    const double clhs104 =     clhs103*clhs3 + clhs46*clhs92 + clhs48*clhs93 + clhs50*clhs94;
    const double clhs105 =     clhs100*clhs46 + clhs103*clhs4 + clhs61*clhs93 + clhs62*clhs94;
    const double clhs106 =     clhs0*clhs1*clhs2*(GPtangent1[0]*clhs3 + GPtangent1[1]*clhs4);
    const double clhs107 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs108 =     clhs0*clhs107*clhs2*(GPnormal[0]*clhs8 + GPnormal[1]*clhs18);
    const double clhs109 =     clhs107*clhs2*clhs45;
    const double clhs110 =     clhs107*clhs45*clhs8;
    const double clhs111 =     DeltaPhi[0][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs112 =     clhs2*clhs45*clhs8;
    const double clhs113 =     clhs107*clhs2*clhs58;
    const double clhs114 =     clhs109*clhs55 + clhs110*clhs47 + clhs111*clhs112 + clhs113*clhs8;
    const double clhs115 =     clhs107*clhs18*clhs45;
    const double clhs116 =     clhs18*clhs2*clhs45;
    const double clhs117 =     clhs109*clhs57 + clhs111*clhs116 + clhs113*clhs18 + clhs115*clhs47;
    const double clhs118 =     DeltaPhi[1][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs119 =     clhs107*clhs2*clhs74;
    const double clhs120 =     clhs109*clhs71 + clhs110*clhs65 + clhs112*clhs118 + clhs119*clhs8;
    const double clhs121 =     clhs109*clhs73 + clhs115*clhs65 + clhs116*clhs118 + clhs119*clhs18;
    const double clhs122 =     DeltaPhi[2][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs123 =     clhs107*clhs2*clhs88;
    const double clhs124 =     clhs109*clhs85 + clhs110*clhs79 + clhs112*clhs122 + clhs123*clhs8;
    const double clhs125 =     clhs109*clhs87 + clhs115*clhs79 + clhs116*clhs122 + clhs123*clhs18;
    const double clhs126 =     DeltaPhi[3][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs127 =     clhs102*clhs107*clhs2;
    const double clhs128 =     clhs109*clhs99 + clhs110*clhs93 + clhs112*clhs126 + clhs127*clhs8;
    const double clhs129 =     clhs101*clhs109 + clhs115*clhs93 + clhs116*clhs126 + clhs127*clhs18;
    const double clhs130 =     clhs0*clhs107*clhs2*(GPtangent1[0]*clhs8 + GPtangent1[1]*clhs18);

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
    lhs(8,0)=clhs28*clhs5;
    lhs(8,1)=clhs33*clhs5;
    lhs(8,2)=clhs38*clhs5;
    lhs(8,3)=clhs43*clhs5;
    lhs(8,4)=clhs0*(GPnormal[0]*clhs60 + GPnormal[1]*clhs63);
    lhs(8,5)=clhs0*(GPnormal[0]*clhs76 + GPnormal[1]*clhs77);
    lhs(8,6)=clhs0*(GPnormal[0]*clhs90 + GPnormal[1]*clhs91);
    lhs(8,7)=clhs0*(GPnormal[0]*clhs104 + GPnormal[1]*clhs105);
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=clhs106*clhs28;
    lhs(9,1)=clhs106*clhs33;
    lhs(9,2)=clhs106*clhs38;
    lhs(9,3)=clhs106*clhs43;
    lhs(9,4)=clhs0*(GPtangent1[0]*clhs60 + GPtangent1[1]*clhs63);
    lhs(9,5)=clhs0*(GPtangent1[0]*clhs76 + GPtangent1[1]*clhs77);
    lhs(9,6)=clhs0*(GPtangent1[0]*clhs90 + GPtangent1[1]*clhs91);
    lhs(9,7)=clhs0*(GPtangent1[0]*clhs104 + GPtangent1[1]*clhs105);
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=clhs108*clhs28;
    lhs(10,1)=clhs108*clhs33;
    lhs(10,2)=clhs108*clhs38;
    lhs(10,3)=clhs108*clhs43;
    lhs(10,4)=clhs0*(GPnormal[0]*clhs114 + GPnormal[1]*clhs117);
    lhs(10,5)=clhs0*(GPnormal[0]*clhs120 + GPnormal[1]*clhs121);
    lhs(10,6)=clhs0*(GPnormal[0]*clhs124 + GPnormal[1]*clhs125);
    lhs(10,7)=clhs0*(GPnormal[0]*clhs128 + GPnormal[1]*clhs129);
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=clhs130*clhs28;
    lhs(11,1)=clhs130*clhs33;
    lhs(11,2)=clhs130*clhs38;
    lhs(11,3)=clhs130*clhs43;
    lhs(11,4)=clhs0*(GPtangent1[0]*clhs114 + GPtangent1[1]*clhs117);
    lhs(11,5)=clhs0*(GPtangent1[0]*clhs120 + GPtangent1[1]*clhs121);
    lhs(11,6)=clhs0*(GPtangent1[0]*clhs124 + GPtangent1[1]*clhs125);
    lhs(11,7)=clhs0*(GPtangent1[0]*clhs128 + GPtangent1[1]*clhs129);
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;

    
    return lhs;
}
    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,12,12> ComputeGaussPointSlipLHS(
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
    bounded_matrix<double,12,12> lhs;
    
    const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
    const double Dt              = rContactData.Dt;
//     const double epsilon_normal  = rContactData.epsilon_normal;
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
    const double clhs29 =     DeltaN1[5][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs30 =     DeltaN1[5][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs31 =     DeltaN2[5][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs32 =     DeltaN2[5][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs33 =     clhs17*(clhs29*clhs3 + clhs30*clhs8) + clhs23*(clhs18*clhs30 + clhs29*clhs4) - clhs24*(clhs14 - clhs19*clhs29 - clhs20*clhs30 + clhs21*clhs31 + clhs22*clhs32) + clhs27*(clhs11*clhs29 + clhs12*clhs30 - clhs13*clhs31 - clhs15*clhs32);
    const double clhs34 =     DeltaN1[6][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs35 =     DeltaN1[6][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs36 =     DeltaN2[6][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs37 =     DeltaN2[6][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs38 =     clhs17*(clhs3*clhs34 + clhs35*clhs8) + clhs23*(clhs18*clhs35 + clhs34*clhs4) + clhs24*(clhs19*clhs34 + clhs20*clhs35 - clhs21*clhs36 - clhs22*clhs37) - clhs27*(-clhs11*clhs34 - clhs12*clhs35 + clhs13*clhs36 + clhs15*clhs37 + clhs16);
    const double clhs39 =     DeltaN1[7][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs40 =     DeltaN1[7][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs41 =     DeltaN2[7][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs42 =     DeltaN2[7][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs43 =     clhs17*(clhs3*clhs39 + clhs40*clhs8) + clhs23*(clhs18*clhs40 + clhs39*clhs4) - clhs24*(clhs16 - clhs19*clhs39 - clhs20*clhs40 + clhs21*clhs41 + clhs22*clhs42) + clhs27*(clhs11*clhs39 + clhs12*clhs40 - clhs13*clhs41 - clhs15*clhs42);
    const double clhs44 =     Deltatangentxis[0](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs45 =     clhs17*clhs27 + clhs23*clhs24;
    const double clhs46 =     clhs1*clhs2*clhs45;
    const double clhs47 =     DeltaJs[0]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs48 =     clhs1*clhs3*clhs45;
    const double clhs49 =     DeltaPhi[0][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs50 =     clhs2*clhs3*clhs45;
    const double clhs51 =     DeltaN1[0][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs52 =     DeltaN1[0][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs53 =     DeltaN2[0][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs54 =     DeltaN2[0][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs55 =     Deltatangentxis[0](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs56 =     Deltatangentxis[0](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs57 =     Deltatangentxis[0](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs58 =     clhs17*(clhs3*clhs51 + clhs44*clhs6 + clhs52*clhs8 + clhs55*clhs9) + clhs23*(clhs18*clhs52 + clhs4*clhs51 + clhs56*clhs6 + clhs57*clhs9) + clhs24*(clhs19*clhs51 + clhs20*clhs52 - clhs21*clhs53 - clhs22*clhs54) + clhs27*(clhs11*clhs51 + clhs12*clhs52 - clhs13*clhs53 - clhs15*clhs54 + clhs6);
    const double clhs59 =     clhs1*clhs2*clhs58;
    const double clhs60 =     clhs3*clhs59 + clhs44*clhs46 + clhs47*clhs48 + clhs49*clhs50;
    const double clhs61 =     clhs1*clhs4*clhs45;
    const double clhs62 =     clhs2*clhs4*clhs45;
    const double clhs63 =     clhs4*clhs59 + clhs46*clhs56 + clhs47*clhs61 + clhs49*clhs62;
    const double clhs64 =     Deltatangentxis[1](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs65 =     DeltaJs[1]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs66 =     DeltaPhi[1][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs67 =     DeltaN1[1][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs68 =     DeltaN1[1][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs69 =     DeltaN2[1][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs70 =     DeltaN2[1][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs71 =     Deltatangentxis[1](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs72 =     Deltatangentxis[1](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs73 =     Deltatangentxis[1](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs74 =     clhs17*(clhs3*clhs67 + clhs6*clhs64 + clhs68*clhs8 + clhs71*clhs9) + clhs23*(clhs18*clhs68 + clhs4*clhs67 + clhs6*clhs72 + clhs73*clhs9) + clhs24*(clhs19*clhs67 + clhs20*clhs68 - clhs21*clhs69 - clhs22*clhs70 + clhs6) + clhs27*(clhs11*clhs67 + clhs12*clhs68 - clhs13*clhs69 - clhs15*clhs70);
    const double clhs75 =     clhs1*clhs2*clhs74;
    const double clhs76 =     clhs3*clhs75 + clhs46*clhs64 + clhs48*clhs65 + clhs50*clhs66;
    const double clhs77 =     clhs4*clhs75 + clhs46*clhs72 + clhs61*clhs65 + clhs62*clhs66;
    const double clhs78 =     Deltatangentxis[2](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs79 =     DeltaJs[2]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs80 =     DeltaPhi[2][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs81 =     DeltaN1[2][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs82 =     DeltaN1[2][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs83 =     DeltaN2[2][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs84 =     DeltaN2[2][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs85 =     Deltatangentxis[2](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs86 =     Deltatangentxis[2](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs87 =     Deltatangentxis[2](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs88 =     clhs17*(clhs3*clhs81 + clhs6*clhs78 + clhs8*clhs82 + clhs85*clhs9) + clhs23*(clhs18*clhs82 + clhs4*clhs81 + clhs6*clhs86 + clhs87*clhs9) + clhs24*(clhs19*clhs81 + clhs20*clhs82 - clhs21*clhs83 - clhs22*clhs84) + clhs27*(clhs11*clhs81 + clhs12*clhs82 - clhs13*clhs83 - clhs15*clhs84 + clhs9);
    const double clhs89 =     clhs1*clhs2*clhs88;
    const double clhs90 =     clhs3*clhs89 + clhs46*clhs78 + clhs48*clhs79 + clhs50*clhs80;
    const double clhs91 =     clhs4*clhs89 + clhs46*clhs86 + clhs61*clhs79 + clhs62*clhs80;
    const double clhs92 =     Deltatangentxis[3](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs93 =     DeltaJs[3]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs94 =     DeltaPhi[3][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs95 =     DeltaN1[3][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs96 =     DeltaN1[3][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs97 =     DeltaN2[3][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs98 =     DeltaN2[3][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs99 =     Deltatangentxis[3](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs100 =     Deltatangentxis[3](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs101 =     Deltatangentxis[3](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs102 =     clhs17*(clhs3*clhs95 + clhs6*clhs92 + clhs8*clhs96 + clhs9*clhs99) + clhs23*(clhs100*clhs6 + clhs101*clhs9 + clhs18*clhs96 + clhs4*clhs95) + clhs24*(clhs19*clhs95 + clhs20*clhs96 - clhs21*clhs97 - clhs22*clhs98 + clhs9) + clhs27*(clhs11*clhs95 + clhs12*clhs96 - clhs13*clhs97 - clhs15*clhs98);
    const double clhs103 =     clhs1*clhs102*clhs2;
    const double clhs104 =     clhs103*clhs3 + clhs46*clhs92 + clhs48*clhs93 + clhs50*clhs94;
    const double clhs105 =     clhs100*clhs46 + clhs103*clhs4 + clhs61*clhs93 + clhs62*clhs94;
    const double clhs106 =     clhs0*clhs1*clhs2*(GPtangent1[0]*clhs3 + GPtangent1[1]*clhs4);
    const double clhs107 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs108 =     clhs0*clhs107*clhs2*(GPnormal[0]*clhs8 + GPnormal[1]*clhs18);
    const double clhs109 =     clhs107*clhs2*clhs45;
    const double clhs110 =     clhs107*clhs45*clhs8;
    const double clhs111 =     DeltaPhi[0][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs112 =     clhs2*clhs45*clhs8;
    const double clhs113 =     clhs107*clhs2*clhs58;
    const double clhs114 =     clhs109*clhs55 + clhs110*clhs47 + clhs111*clhs112 + clhs113*clhs8;
    const double clhs115 =     clhs107*clhs18*clhs45;
    const double clhs116 =     clhs18*clhs2*clhs45;
    const double clhs117 =     clhs109*clhs57 + clhs111*clhs116 + clhs113*clhs18 + clhs115*clhs47;
    const double clhs118 =     DeltaPhi[1][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs119 =     clhs107*clhs2*clhs74;
    const double clhs120 =     clhs109*clhs71 + clhs110*clhs65 + clhs112*clhs118 + clhs119*clhs8;
    const double clhs121 =     clhs109*clhs73 + clhs115*clhs65 + clhs116*clhs118 + clhs119*clhs18;
    const double clhs122 =     DeltaPhi[2][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs123 =     clhs107*clhs2*clhs88;
    const double clhs124 =     clhs109*clhs85 + clhs110*clhs79 + clhs112*clhs122 + clhs123*clhs8;
    const double clhs125 =     clhs109*clhs87 + clhs115*clhs79 + clhs116*clhs122 + clhs123*clhs18;
    const double clhs126 =     DeltaPhi[3][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs127 =     clhs102*clhs107*clhs2;
    const double clhs128 =     clhs109*clhs99 + clhs110*clhs93 + clhs112*clhs126 + clhs127*clhs8;
    const double clhs129 =     clhs101*clhs109 + clhs115*clhs93 + clhs116*clhs126 + clhs127*clhs18;
    const double clhs130 =     clhs0*clhs107*clhs2*(GPtangent1[0]*clhs8 + GPtangent1[1]*clhs18);

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
    lhs(8,0)=clhs28*clhs5;
    lhs(8,1)=clhs33*clhs5;
    lhs(8,2)=clhs38*clhs5;
    lhs(8,3)=clhs43*clhs5;
    lhs(8,4)=clhs0*(GPnormal[0]*clhs60 + GPnormal[1]*clhs63);
    lhs(8,5)=clhs0*(GPnormal[0]*clhs76 + GPnormal[1]*clhs77);
    lhs(8,6)=clhs0*(GPnormal[0]*clhs90 + GPnormal[1]*clhs91);
    lhs(8,7)=clhs0*(GPnormal[0]*clhs104 + GPnormal[1]*clhs105);
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=clhs106*clhs28;
    lhs(9,1)=clhs106*clhs33;
    lhs(9,2)=clhs106*clhs38;
    lhs(9,3)=clhs106*clhs43;
    lhs(9,4)=clhs0*(GPtangent1[0]*clhs60 + GPtangent1[1]*clhs63);
    lhs(9,5)=clhs0*(GPtangent1[0]*clhs76 + GPtangent1[1]*clhs77);
    lhs(9,6)=clhs0*(GPtangent1[0]*clhs90 + GPtangent1[1]*clhs91);
    lhs(9,7)=clhs0*(GPtangent1[0]*clhs104 + GPtangent1[1]*clhs105);
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=clhs108*clhs28;
    lhs(10,1)=clhs108*clhs33;
    lhs(10,2)=clhs108*clhs38;
    lhs(10,3)=clhs108*clhs43;
    lhs(10,4)=clhs0*(GPnormal[0]*clhs114 + GPnormal[1]*clhs117);
    lhs(10,5)=clhs0*(GPnormal[0]*clhs120 + GPnormal[1]*clhs121);
    lhs(10,6)=clhs0*(GPnormal[0]*clhs124 + GPnormal[1]*clhs125);
    lhs(10,7)=clhs0*(GPnormal[0]*clhs128 + GPnormal[1]*clhs129);
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=clhs130*clhs28;
    lhs(11,1)=clhs130*clhs33;
    lhs(11,2)=clhs130*clhs38;
    lhs(11,3)=clhs130*clhs43;
    lhs(11,4)=clhs0*(GPtangent1[0]*clhs114 + GPtangent1[1]*clhs117);
    lhs(11,5)=clhs0*(GPtangent1[0]*clhs120 + GPtangent1[1]*clhs121);
    lhs(11,6)=clhs0*(GPtangent1[0]*clhs124 + GPtangent1[1]*clhs125);
    lhs(11,7)=clhs0*(GPtangent1[0]*clhs128 + GPtangent1[1]*clhs129);
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;

    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,12,12> ComputeGaussPointInactiveLHS(
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
    bounded_matrix<double,12,12> lhs;
    
//     const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
//     const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
//     const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
//     const double Dt              = rContactData.Dt;
//     const double epsilon_normal  = rContactData.epsilon_normal;
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
    
    static inline array_1d<double,12> ComputeGaussPointActiveRHS(
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
    array_1d<double,12> rhs;
    
    const bounded_matrix<double, 2, 2> normalslave    = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave      = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm             = rContactData.LagrangeMultipliers;
//     const double Dt             = rContactData.Dt;
//     const double epsilon_normal = rContactData.epsilon_normal;
    
    const array_1d<double,2> GPnormal     = prod(trans(normalslave), N1);
    const array_1d<double,2> GPtangent1   = prod(trans(tan1slave), N1);
    
    const bounded_matrix<double, 2, 2> v1 = rContactData.v1;
    const bounded_matrix<double, 2, 2> v2 = rContactData.v2;
    
    const double crhs0 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs1 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs2 =     Phi[0]; // PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs3 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs4 =     crhs1*(crhs2*lm(0,0) + crhs3*lm(1,0));
    const double crhs5 =     crhs1*(crhs2*lm(0,1) + crhs3*lm(1,1));
    const double crhs6 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs7 =     N1[0]; // N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs8 =     N1[1]; // N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs9 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs10 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs11 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs12 =     crhs1*crhs11*crhs2;
    const double crhs13 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs14 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs15 =     crhs1*crhs11*crhs3;

    rhs[0]=-crhs0*crhs4;
    rhs[1]=-crhs0*crhs5;
    rhs[2]=-crhs4*crhs6;
    rhs[3]=-crhs5*crhs6;
    rhs[4]=crhs4*crhs7;
    rhs[5]=crhs5*crhs7;
    rhs[6]=crhs4*crhs8;
    rhs[7]=crhs5*crhs8;
    rhs[8]=crhs12*(GPnormal[0]*crhs9 + GPnormal[1]*crhs10);
    rhs[9]=crhs12*(GPtangent1[0]*crhs9 + GPtangent1[1]*crhs10);
    rhs[10]=crhs15*(GPnormal[0]*crhs13 + GPnormal[1]*crhs14);
    rhs[11]=crhs15*(GPtangent1[0]*crhs13 + GPtangent1[1]*crhs14);

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,12> ComputeGaussPointStickRHS(
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
    array_1d<double,12> rhs;
    
    const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
    const double Dt              = rContactData.Dt;
//     const double epsilon_normal  = rContactData.epsilon_normal;
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
    const double crhs12 =     crhs10*crhs2*crhs3*Phi[1]; // CRHS10*CRHS2*CRHS3*PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))

    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=-crhs11*(GPnormal[0]*crhs0 + GPnormal[1]*crhs1);
    rhs[9]=-crhs11*(GPtangent1[0]*crhs0 + GPtangent1[1]*crhs1);
    rhs[10]=-crhs12*(GPnormal[0]*crhs5 + GPnormal[1]*crhs9);
    rhs[11]=-crhs12*(GPtangent1[0]*crhs5 + GPtangent1[1]*crhs9);

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,12> ComputeGaussPointSlipRHS(
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
    array_1d<double,12> rhs;
    
    const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
    const double Dt              = rContactData.Dt;
//     const double epsilon_normal  = rContactData.epsilon_normal;
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
    const double crhs12 =     crhs10*crhs2*crhs3*Phi[1]; // CRHS10*CRHS2*CRHS3*PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))

    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=-crhs11*(GPnormal[0]*crhs0 + GPnormal[1]*crhs1);
    rhs[9]=-crhs11*(GPtangent1[0]*crhs0 + GPtangent1[1]*crhs1);
    rhs[10]=-crhs12*(GPnormal[0]*crhs5 + GPnormal[1]*crhs9);
    rhs[11]=-crhs12*(GPtangent1[0]*crhs5 + GPtangent1[1]*crhs9);

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,12> ComputeGaussPointInactiveRHS(
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
    array_1d<double,12> rhs;
    
//     const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
//     const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
//     const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
//     const double Dt              = rContactData.Dt;
//     const double epsilon_normal  = rContactData.epsilon_normal;
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
};// class Contact2D2N2N
}
#endif /* KRATOS_CONTACT2D2N2N defined */
