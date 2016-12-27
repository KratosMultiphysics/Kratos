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
        const ContactData<2,2>& rContactData
        )
{
    bounded_matrix<double,12,12> lhs;
    
    const double integration_point_gap = inner_prod(rContactData.Gaps,N1);
    const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
    
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
    const double clhs86 =     clhs2*clhs3*clhs85;
    const double clhs87 =     DeltaGap[5]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs88 =     DeltaGap[6]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs89 =     DeltaGap[7]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs90 =     clhs2*clhs3*clhs83;
    const double clhs91 =     DeltaGap[0]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs92 =     clhs3*clhs83*clhs85;
    const double clhs93 =     clhs2*clhs83*clhs85;
    const double clhs94 =     DeltaGap[1]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs95 =     DeltaGap[2]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs96 =     DeltaGap[3]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs97 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs98 =     clhs2*clhs3*clhs97;
    const double clhs99 =     clhs3*clhs83*clhs97;
    const double clhs100 =     clhs2*clhs83*clhs97;
    const double clhs101 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs102 =     clhs101*clhs2*clhs4;
    const double clhs103 =     clhs2*clhs4*clhs83;
    const double clhs104 =     clhs101*clhs4*clhs83;
    const double clhs105 =     clhs101*clhs2*clhs83;
    const double clhs106 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs107 =     clhs106*clhs2*clhs4;
    const double clhs108 =     clhs106*clhs4*clhs83;
    const double clhs109 =     clhs106*clhs2*clhs83;

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
    lhs(8,0)=-clhs84*clhs86;
    lhs(8,1)=-clhs86*clhs87;
    lhs(8,2)=-clhs86*clhs88;
    lhs(8,3)=-clhs86*clhs89;
    lhs(8,4)=-clhs10*clhs92 - clhs14*clhs93 - clhs86*clhs91 - clhs90*DeltaNormals[0](0,0);
    lhs(8,5)=-clhs17*clhs92 - clhs19*clhs93 - clhs86*clhs94 - clhs90*DeltaNormals[1](0,0);
    lhs(8,6)=-clhs22*clhs92 - clhs24*clhs93 - clhs86*clhs95 - clhs90*DeltaNormals[2](0,0);
    lhs(8,7)=-clhs27*clhs92 - clhs29*clhs93 - clhs86*clhs96 - clhs90*DeltaNormals[3](0,0);
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=-clhs84*clhs98;
    lhs(9,1)=-clhs87*clhs98;
    lhs(9,2)=-clhs88*clhs98;
    lhs(9,3)=-clhs89*clhs98;
    lhs(9,4)=-clhs10*clhs99 - clhs100*clhs14 - clhs90*DeltaNormals[0](0,1) - clhs91*clhs98;
    lhs(9,5)=-clhs100*clhs19 - clhs17*clhs99 - clhs90*DeltaNormals[1](0,1) - clhs94*clhs98;
    lhs(9,6)=-clhs100*clhs24 - clhs22*clhs99 - clhs90*DeltaNormals[2](0,1) - clhs95*clhs98;
    lhs(9,7)=-clhs100*clhs29 - clhs27*clhs99 - clhs90*DeltaNormals[3](0,1) - clhs96*clhs98;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=-clhs102*clhs84;
    lhs(10,1)=-clhs102*clhs87;
    lhs(10,2)=-clhs102*clhs88;
    lhs(10,3)=-clhs102*clhs89;
    lhs(10,4)=-clhs10*clhs104 - clhs102*clhs91 - clhs103*DeltaNormals[0](1,0) - clhs105*clhs15;
    lhs(10,5)=-clhs102*clhs94 - clhs103*DeltaNormals[1](1,0) - clhs104*clhs17 - clhs105*clhs20;
    lhs(10,6)=-clhs102*clhs95 - clhs103*DeltaNormals[2](1,0) - clhs104*clhs22 - clhs105*clhs25;
    lhs(10,7)=-clhs102*clhs96 - clhs103*DeltaNormals[3](1,0) - clhs104*clhs27 - clhs105*clhs30;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=-clhs107*clhs84;
    lhs(11,1)=-clhs107*clhs87;
    lhs(11,2)=-clhs107*clhs88;
    lhs(11,3)=-clhs107*clhs89;
    lhs(11,4)=-clhs10*clhs108 - clhs103*DeltaNormals[0](1,1) - clhs107*clhs91 - clhs109*clhs15;
    lhs(11,5)=-clhs103*DeltaNormals[1](1,1) - clhs107*clhs94 - clhs108*clhs17 - clhs109*clhs20;
    lhs(11,6)=-clhs103*DeltaNormals[2](1,1) - clhs107*clhs95 - clhs108*clhs22 - clhs109*clhs25;
    lhs(11,7)=-clhs103*DeltaNormals[3](1,1) - clhs107*clhs96 - clhs108*clhs27 - clhs109*clhs30;
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
        const ContactData<2,2>& rContactData
        )
{
    bounded_matrix<double,12,12> lhs;
    
    const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
    const double Dt              = rContactData.Dt;
    
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
    const double clhs29 =     DeltaN1[5][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs30 =     DeltaN1[5][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs31 =     DeltaN2[5][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs32 =     DeltaN2[5][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs33 =     clhs15*(clhs0*clhs29 + clhs30*clhs6) + clhs22*(clhs16*clhs29 + clhs17*clhs30) - clhs23*(clhs12 - clhs18*clhs29 - clhs19*clhs30 + clhs20*clhs31 + clhs21*clhs32) + clhs26*(clhs10*clhs30 - clhs11*clhs31 - clhs13*clhs32 + clhs29*clhs9);
    const double clhs34 =     clhs1*clhs2*clhs3*clhs33;
    const double clhs35 =     DeltaN1[6][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs36 =     DeltaN1[6][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs37 =     DeltaN2[6][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs38 =     DeltaN2[6][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs39 =     clhs15*(clhs0*clhs35 + clhs36*clhs6) + clhs22*(clhs16*clhs35 + clhs17*clhs36) + clhs23*(clhs18*clhs35 + clhs19*clhs36 - clhs20*clhs37 - clhs21*clhs38) - clhs26*(-clhs10*clhs36 + clhs11*clhs37 + clhs13*clhs38 + clhs14 - clhs35*clhs9);
    const double clhs40 =     clhs1*clhs2*clhs3*clhs39;
    const double clhs41 =     DeltaN1[7][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs42 =     DeltaN1[7][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs43 =     DeltaN2[7][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs44 =     DeltaN2[7][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs45 =     clhs15*(clhs0*clhs41 + clhs42*clhs6) + clhs22*(clhs16*clhs41 + clhs17*clhs42) - clhs23*(clhs14 - clhs18*clhs41 - clhs19*clhs42 + clhs20*clhs43 + clhs21*clhs44) + clhs26*(clhs10*clhs42 - clhs11*clhs43 - clhs13*clhs44 + clhs41*clhs9);
    const double clhs46 =     clhs1*clhs2*clhs3*clhs45;
    const double clhs47 =     Deltatangentxis[0](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs48 =     clhs15*clhs26 + clhs22*clhs23;
    const double clhs49 =     clhs2*clhs3*clhs48;
    const double clhs50 =     DeltaJs[0]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs51 =     clhs0*clhs2*clhs48;
    const double clhs52 =     DeltaPhi[0][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs53 =     clhs0*clhs3*clhs48;
    const double clhs54 =     DeltaN1[0][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs55 =     DeltaN1[0][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs56 =     DeltaN2[0][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs57 =     DeltaN2[0][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs58 =     Deltatangentxis[0](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs59 =     Deltatangentxis[0](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs60 =     Deltatangentxis[0](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs61 =     clhs15*(clhs0*clhs54 + clhs4*clhs47 + clhs55*clhs6 + clhs58*clhs7) + clhs22*(clhs16*clhs54 + clhs17*clhs55 + clhs4*clhs59 + clhs60*clhs7) + clhs23*(clhs18*clhs54 + clhs19*clhs55 - clhs20*clhs56 - clhs21*clhs57) + clhs26*(clhs10*clhs55 - clhs11*clhs56 - clhs13*clhs57 + clhs4 + clhs54*clhs9);
    const double clhs62 =     clhs2*clhs3*clhs61;
    const double clhs63 =     Deltatangentxis[1](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs64 =     DeltaJs[1]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs65 =     DeltaPhi[1][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs66 =     DeltaN1[1][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs67 =     DeltaN1[1][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs68 =     DeltaN2[1][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs69 =     DeltaN2[1][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs70 =     Deltatangentxis[1](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs71 =     Deltatangentxis[1](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs72 =     Deltatangentxis[1](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs73 =     clhs15*(clhs0*clhs66 + clhs4*clhs63 + clhs6*clhs67 + clhs7*clhs70) + clhs22*(clhs16*clhs66 + clhs17*clhs67 + clhs4*clhs71 + clhs7*clhs72) + clhs23*(clhs18*clhs66 + clhs19*clhs67 - clhs20*clhs68 - clhs21*clhs69 + clhs4) + clhs26*(clhs10*clhs67 - clhs11*clhs68 - clhs13*clhs69 + clhs66*clhs9);
    const double clhs74 =     clhs2*clhs3*clhs73;
    const double clhs75 =     Deltatangentxis[2](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs76 =     DeltaJs[2]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs77 =     DeltaPhi[2][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs78 =     DeltaN1[2][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs79 =     DeltaN1[2][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs80 =     DeltaN2[2][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs81 =     DeltaN2[2][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs82 =     Deltatangentxis[2](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs83 =     Deltatangentxis[2](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs84 =     Deltatangentxis[2](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs85 =     clhs15*(clhs0*clhs78 + clhs4*clhs75 + clhs6*clhs79 + clhs7*clhs82) + clhs22*(clhs16*clhs78 + clhs17*clhs79 + clhs4*clhs83 + clhs7*clhs84) + clhs23*(clhs18*clhs78 + clhs19*clhs79 - clhs20*clhs80 - clhs21*clhs81) + clhs26*(clhs10*clhs79 - clhs11*clhs80 - clhs13*clhs81 + clhs7 + clhs78*clhs9);
    const double clhs86 =     clhs2*clhs3*clhs85;
    const double clhs87 =     Deltatangentxis[3](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs88 =     DeltaJs[3]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs89 =     DeltaPhi[3][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs90 =     DeltaN1[3][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs91 =     DeltaN1[3][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs92 =     DeltaN2[3][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs93 =     DeltaN2[3][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs94 =     Deltatangentxis[3](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs95 =     Deltatangentxis[3](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs96 =     Deltatangentxis[3](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs97 =     clhs15*(clhs0*clhs90 + clhs4*clhs87 + clhs6*clhs91 + clhs7*clhs94) + clhs22*(clhs16*clhs90 + clhs17*clhs91 + clhs4*clhs95 + clhs7*clhs96) + clhs23*(clhs18*clhs90 + clhs19*clhs91 - clhs20*clhs92 - clhs21*clhs93 + clhs7) + clhs26*(clhs10*clhs91 - clhs11*clhs92 - clhs13*clhs93 + clhs9*clhs90);
    const double clhs98 =     clhs2*clhs3*clhs97;
    const double clhs99 =     clhs16*clhs2*clhs48;
    const double clhs100 =     clhs16*clhs3*clhs48;
    const double clhs101 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs102 =     clhs1*clhs101*clhs27*clhs3;
    const double clhs103 =     clhs1*clhs101*clhs3*clhs33;
    const double clhs104 =     clhs1*clhs101*clhs3*clhs39;
    const double clhs105 =     clhs1*clhs101*clhs3*clhs45;
    const double clhs106 =     clhs101*clhs3*clhs48;
    const double clhs107 =     clhs101*clhs48*clhs6;
    const double clhs108 =     DeltaPhi[0][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs109 =     clhs3*clhs48*clhs6;
    const double clhs110 =     clhs101*clhs3*clhs61;
    const double clhs111 =     DeltaPhi[1][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs112 =     clhs101*clhs3*clhs73;
    const double clhs113 =     DeltaPhi[2][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs114 =     clhs101*clhs3*clhs85;
    const double clhs115 =     DeltaPhi[3][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs116 =     clhs101*clhs3*clhs97;
    const double clhs117 =     clhs101*clhs17*clhs48;
    const double clhs118 =     clhs17*clhs3*clhs48;

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
    lhs(8,0)=clhs0*clhs28;
    lhs(8,1)=clhs0*clhs34;
    lhs(8,2)=clhs0*clhs40;
    lhs(8,3)=clhs0*clhs46;
    lhs(8,4)=clhs1*(clhs0*clhs62 + clhs47*clhs49 + clhs50*clhs51 + clhs52*clhs53);
    lhs(8,5)=clhs1*(clhs0*clhs74 + clhs49*clhs63 + clhs51*clhs64 + clhs53*clhs65);
    lhs(8,6)=clhs1*(clhs0*clhs86 + clhs49*clhs75 + clhs51*clhs76 + clhs53*clhs77);
    lhs(8,7)=clhs1*(clhs0*clhs98 + clhs49*clhs87 + clhs51*clhs88 + clhs53*clhs89);
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=clhs16*clhs28;
    lhs(9,1)=clhs16*clhs34;
    lhs(9,2)=clhs16*clhs40;
    lhs(9,3)=clhs16*clhs46;
    lhs(9,4)=clhs1*(clhs100*clhs52 + clhs16*clhs62 + clhs49*clhs59 + clhs50*clhs99);
    lhs(9,5)=clhs1*(clhs100*clhs65 + clhs16*clhs74 + clhs49*clhs71 + clhs64*clhs99);
    lhs(9,6)=clhs1*(clhs100*clhs77 + clhs16*clhs86 + clhs49*clhs83 + clhs76*clhs99);
    lhs(9,7)=clhs1*(clhs100*clhs89 + clhs16*clhs98 + clhs49*clhs95 + clhs88*clhs99);
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=clhs102*clhs6;
    lhs(10,1)=clhs103*clhs6;
    lhs(10,2)=clhs104*clhs6;
    lhs(10,3)=clhs105*clhs6;
    lhs(10,4)=clhs1*(clhs106*clhs58 + clhs107*clhs50 + clhs108*clhs109 + clhs110*clhs6);
    lhs(10,5)=clhs1*(clhs106*clhs70 + clhs107*clhs64 + clhs109*clhs111 + clhs112*clhs6);
    lhs(10,6)=clhs1*(clhs106*clhs82 + clhs107*clhs76 + clhs109*clhs113 + clhs114*clhs6);
    lhs(10,7)=clhs1*(clhs106*clhs94 + clhs107*clhs88 + clhs109*clhs115 + clhs116*clhs6);
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=clhs102*clhs17;
    lhs(11,1)=clhs103*clhs17;
    lhs(11,2)=clhs104*clhs17;
    lhs(11,3)=clhs105*clhs17;
    lhs(11,4)=clhs1*(clhs106*clhs60 + clhs108*clhs118 + clhs110*clhs17 + clhs117*clhs50);
    lhs(11,5)=clhs1*(clhs106*clhs72 + clhs111*clhs118 + clhs112*clhs17 + clhs117*clhs64);
    lhs(11,6)=clhs1*(clhs106*clhs84 + clhs113*clhs118 + clhs114*clhs17 + clhs117*clhs76);
    lhs(11,7)=clhs1*(clhs106*clhs96 + clhs115*clhs118 + clhs116*clhs17 + clhs117*clhs88);
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
        const ContactData<2,2>& rContactData
        )
{
    bounded_matrix<double,12,12> lhs;
    
    const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
    const double Dt              = rContactData.Dt;

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
    const double clhs29 =     DeltaN1[5][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs30 =     DeltaN1[5][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs31 =     DeltaN2[5][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs32 =     DeltaN2[5][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs33 =     clhs15*(clhs0*clhs29 + clhs30*clhs6) + clhs22*(clhs16*clhs29 + clhs17*clhs30) - clhs23*(clhs12 - clhs18*clhs29 - clhs19*clhs30 + clhs20*clhs31 + clhs21*clhs32) + clhs26*(clhs10*clhs30 - clhs11*clhs31 - clhs13*clhs32 + clhs29*clhs9);
    const double clhs34 =     clhs1*clhs2*clhs3*clhs33;
    const double clhs35 =     DeltaN1[6][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs36 =     DeltaN1[6][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs37 =     DeltaN2[6][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs38 =     DeltaN2[6][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs39 =     clhs15*(clhs0*clhs35 + clhs36*clhs6) + clhs22*(clhs16*clhs35 + clhs17*clhs36) + clhs23*(clhs18*clhs35 + clhs19*clhs36 - clhs20*clhs37 - clhs21*clhs38) - clhs26*(-clhs10*clhs36 + clhs11*clhs37 + clhs13*clhs38 + clhs14 - clhs35*clhs9);
    const double clhs40 =     clhs1*clhs2*clhs3*clhs39;
    const double clhs41 =     DeltaN1[7][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs42 =     DeltaN1[7][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs43 =     DeltaN2[7][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs44 =     DeltaN2[7][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs45 =     clhs15*(clhs0*clhs41 + clhs42*clhs6) + clhs22*(clhs16*clhs41 + clhs17*clhs42) - clhs23*(clhs14 - clhs18*clhs41 - clhs19*clhs42 + clhs20*clhs43 + clhs21*clhs44) + clhs26*(clhs10*clhs42 - clhs11*clhs43 - clhs13*clhs44 + clhs41*clhs9);
    const double clhs46 =     clhs1*clhs2*clhs3*clhs45;
    const double clhs47 =     Deltatangentxis[0](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs48 =     clhs15*clhs26 + clhs22*clhs23;
    const double clhs49 =     clhs2*clhs3*clhs48;
    const double clhs50 =     DeltaJs[0]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs51 =     clhs0*clhs2*clhs48;
    const double clhs52 =     DeltaPhi[0][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs53 =     clhs0*clhs3*clhs48;
    const double clhs54 =     DeltaN1[0][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs55 =     DeltaN1[0][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs56 =     DeltaN2[0][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs57 =     DeltaN2[0][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs58 =     Deltatangentxis[0](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs59 =     Deltatangentxis[0](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs60 =     Deltatangentxis[0](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs61 =     clhs15*(clhs0*clhs54 + clhs4*clhs47 + clhs55*clhs6 + clhs58*clhs7) + clhs22*(clhs16*clhs54 + clhs17*clhs55 + clhs4*clhs59 + clhs60*clhs7) + clhs23*(clhs18*clhs54 + clhs19*clhs55 - clhs20*clhs56 - clhs21*clhs57) + clhs26*(clhs10*clhs55 - clhs11*clhs56 - clhs13*clhs57 + clhs4 + clhs54*clhs9);
    const double clhs62 =     clhs2*clhs3*clhs61;
    const double clhs63 =     Deltatangentxis[1](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs64 =     DeltaJs[1]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs65 =     DeltaPhi[1][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs66 =     DeltaN1[1][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs67 =     DeltaN1[1][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs68 =     DeltaN2[1][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs69 =     DeltaN2[1][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs70 =     Deltatangentxis[1](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs71 =     Deltatangentxis[1](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs72 =     Deltatangentxis[1](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs73 =     clhs15*(clhs0*clhs66 + clhs4*clhs63 + clhs6*clhs67 + clhs7*clhs70) + clhs22*(clhs16*clhs66 + clhs17*clhs67 + clhs4*clhs71 + clhs7*clhs72) + clhs23*(clhs18*clhs66 + clhs19*clhs67 - clhs20*clhs68 - clhs21*clhs69 + clhs4) + clhs26*(clhs10*clhs67 - clhs11*clhs68 - clhs13*clhs69 + clhs66*clhs9);
    const double clhs74 =     clhs2*clhs3*clhs73;
    const double clhs75 =     Deltatangentxis[2](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs76 =     DeltaJs[2]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs77 =     DeltaPhi[2][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs78 =     DeltaN1[2][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs79 =     DeltaN1[2][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs80 =     DeltaN2[2][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs81 =     DeltaN2[2][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs82 =     Deltatangentxis[2](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs83 =     Deltatangentxis[2](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs84 =     Deltatangentxis[2](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs85 =     clhs15*(clhs0*clhs78 + clhs4*clhs75 + clhs6*clhs79 + clhs7*clhs82) + clhs22*(clhs16*clhs78 + clhs17*clhs79 + clhs4*clhs83 + clhs7*clhs84) + clhs23*(clhs18*clhs78 + clhs19*clhs79 - clhs20*clhs80 - clhs21*clhs81) + clhs26*(clhs10*clhs79 - clhs11*clhs80 - clhs13*clhs81 + clhs7 + clhs78*clhs9);
    const double clhs86 =     clhs2*clhs3*clhs85;
    const double clhs87 =     Deltatangentxis[3](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs88 =     DeltaJs[3]; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs89 =     DeltaPhi[3][0]; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs90 =     DeltaN1[3][0]; // DERIVATIVE(N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs91 =     DeltaN1[3][1]; // DERIVATIVE(N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs92 =     DeltaN2[3][0]; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs93 =     DeltaN2[3][1]; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs94 =     Deltatangentxis[3](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs95 =     Deltatangentxis[3](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs96 =     Deltatangentxis[3](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs97 =     clhs15*(clhs0*clhs90 + clhs4*clhs87 + clhs6*clhs91 + clhs7*clhs94) + clhs22*(clhs16*clhs90 + clhs17*clhs91 + clhs4*clhs95 + clhs7*clhs96) + clhs23*(clhs18*clhs90 + clhs19*clhs91 - clhs20*clhs92 - clhs21*clhs93 + clhs7) + clhs26*(clhs10*clhs91 - clhs11*clhs92 - clhs13*clhs93 + clhs9*clhs90);
    const double clhs98 =     clhs2*clhs3*clhs97;
    const double clhs99 =     clhs16*clhs2*clhs48;
    const double clhs100 =     clhs16*clhs3*clhs48;
    const double clhs101 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs102 =     clhs1*clhs101*clhs27*clhs3;
    const double clhs103 =     clhs1*clhs101*clhs3*clhs33;
    const double clhs104 =     clhs1*clhs101*clhs3*clhs39;
    const double clhs105 =     clhs1*clhs101*clhs3*clhs45;
    const double clhs106 =     clhs101*clhs3*clhs48;
    const double clhs107 =     clhs101*clhs48*clhs6;
    const double clhs108 =     DeltaPhi[0][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs109 =     clhs3*clhs48*clhs6;
    const double clhs110 =     clhs101*clhs3*clhs61;
    const double clhs111 =     DeltaPhi[1][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs112 =     clhs101*clhs3*clhs73;
    const double clhs113 =     DeltaPhi[2][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs114 =     clhs101*clhs3*clhs85;
    const double clhs115 =     DeltaPhi[3][1]; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs116 =     clhs101*clhs3*clhs97;
    const double clhs117 =     clhs101*clhs17*clhs48;
    const double clhs118 =     clhs17*clhs3*clhs48;

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
    lhs(8,0)=clhs0*clhs28;
    lhs(8,1)=clhs0*clhs34;
    lhs(8,2)=clhs0*clhs40;
    lhs(8,3)=clhs0*clhs46;
    lhs(8,4)=clhs1*(clhs0*clhs62 + clhs47*clhs49 + clhs50*clhs51 + clhs52*clhs53);
    lhs(8,5)=clhs1*(clhs0*clhs74 + clhs49*clhs63 + clhs51*clhs64 + clhs53*clhs65);
    lhs(8,6)=clhs1*(clhs0*clhs86 + clhs49*clhs75 + clhs51*clhs76 + clhs53*clhs77);
    lhs(8,7)=clhs1*(clhs0*clhs98 + clhs49*clhs87 + clhs51*clhs88 + clhs53*clhs89);
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=clhs16*clhs28;
    lhs(9,1)=clhs16*clhs34;
    lhs(9,2)=clhs16*clhs40;
    lhs(9,3)=clhs16*clhs46;
    lhs(9,4)=clhs1*(clhs100*clhs52 + clhs16*clhs62 + clhs49*clhs59 + clhs50*clhs99);
    lhs(9,5)=clhs1*(clhs100*clhs65 + clhs16*clhs74 + clhs49*clhs71 + clhs64*clhs99);
    lhs(9,6)=clhs1*(clhs100*clhs77 + clhs16*clhs86 + clhs49*clhs83 + clhs76*clhs99);
    lhs(9,7)=clhs1*(clhs100*clhs89 + clhs16*clhs98 + clhs49*clhs95 + clhs88*clhs99);
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=clhs102*clhs6;
    lhs(10,1)=clhs103*clhs6;
    lhs(10,2)=clhs104*clhs6;
    lhs(10,3)=clhs105*clhs6;
    lhs(10,4)=clhs1*(clhs106*clhs58 + clhs107*clhs50 + clhs108*clhs109 + clhs110*clhs6);
    lhs(10,5)=clhs1*(clhs106*clhs70 + clhs107*clhs64 + clhs109*clhs111 + clhs112*clhs6);
    lhs(10,6)=clhs1*(clhs106*clhs82 + clhs107*clhs76 + clhs109*clhs113 + clhs114*clhs6);
    lhs(10,7)=clhs1*(clhs106*clhs94 + clhs107*clhs88 + clhs109*clhs115 + clhs116*clhs6);
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=clhs102*clhs17;
    lhs(11,1)=clhs103*clhs17;
    lhs(11,2)=clhs104*clhs17;
    lhs(11,3)=clhs105*clhs17;
    lhs(11,4)=clhs1*(clhs106*clhs60 + clhs108*clhs118 + clhs110*clhs17 + clhs117*clhs50);
    lhs(11,5)=clhs1*(clhs106*clhs72 + clhs111*clhs118 + clhs112*clhs17 + clhs117*clhs64);
    lhs(11,6)=clhs1*(clhs106*clhs84 + clhs113*clhs118 + clhs114*clhs17 + clhs117*clhs76);
    lhs(11,7)=clhs1*(clhs106*clhs96 + clhs115*clhs118 + clhs116*clhs17 + clhs117*clhs88);
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;

    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,12> ComputeGaussPointActiveRHS(
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const ContactData<2,2>& rContactData
        )
{
    array_1d<double,12> rhs;
    
    const double integration_point_gap = inner_prod(rContactData.Gaps,N1);
    const bounded_matrix<double, 2, 2> normalslave    = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave      = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm             = rContactData.LagrangeMultipliers;
    
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
    const double crhs9 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs10 =     crhs1*crhs2*crhs9;
    const double crhs11 =     crhs1*crhs3*crhs9;

    rhs[0]=-crhs0*crhs4;
    rhs[1]=-crhs0*crhs5;
    rhs[2]=-crhs4*crhs6;
    rhs[3]=-crhs5*crhs6;
    rhs[4]=crhs4*crhs7;
    rhs[5]=crhs5*crhs7;
    rhs[6]=crhs4*crhs8;
    rhs[7]=crhs5*crhs8;
    rhs[8]=crhs10*normalslave(0,0);
    rhs[9]=crhs10*normalslave(0,1);
    rhs[10]=crhs11*normalslave(1,0);
    rhs[11]=crhs11*normalslave(1,1);

    
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
        const ContactData<2,2>& rContactData
        )
{
    array_1d<double,12> rhs;
    
    const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
    const double Dt              = rContactData.Dt;

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
    const double crhs12 =     crhs1*crhs10*crhs2*Phi[1]; // CRHS1*CRHS10*CRHS2*PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))

    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=-crhs0*crhs11;
    rhs[9]=-crhs11*crhs8;
    rhs[10]=-crhs12*crhs4;
    rhs[11]=-crhs12*crhs9;

    
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
        const ContactData<2,2>& rContactData
        )
{
    array_1d<double,12> rhs;
    
    const bounded_matrix<double, 2, 2> normalslave     = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave       = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm              = rContactData.LagrangeMultipliers;
    const double Dt              = rContactData.Dt;
    
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
    const double crhs12 =     crhs1*crhs10*crhs2*Phi[1]; // CRHS1*CRHS10*CRHS2*PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))

    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=-crhs0*crhs11;
    rhs[9]=-crhs11*crhs8;
    rhs[10]=-crhs12*crhs4;
    rhs[11]=-crhs12*crhs9;

    
    return rhs;
}
private:
};// class Contact2D2N2N
}
#endif /* KRATOS_CONTACT2D2N2N defined */
