// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
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
    
    static inline bounded_matrix<double,12,12> ComputeGaussPointFrictionlessLHS(
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const ContactData<2,2>& rContactData
        )
{
    bounded_matrix<double,12,12> lhs;
    
    const double integration_point_gap = inner_prod(rContactData.Gaps,N1);
    const bounded_matrix<double, 2, 2> normalslave = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave   = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm          = rContactData.LagrangeMultipliers;

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
    const double clhs83 =     clhs2*clhs3;
    const double clhs84 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs85 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs86 =     DeltaGap[4]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs87 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs88 =     1.0/Dt;
    const double clhs89 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs90 =     Dt*v1(0,0);
    const double clhs91 =     Dt*v1(1,0);
    const double clhs92 =     Dt*v2(0,0);
    const double clhs93 =     Dt*v2(1,0);
    const double clhs94 =     -clhs0*clhs92 - clhs41*clhs93 + clhs55*clhs90 + clhs69*clhs91;
    const double clhs95 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs96 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs97 =     Dt*v1(0,1);
    const double clhs98 =     Dt*v1(1,1);
    const double clhs99 =     Dt*v2(0,1);
    const double clhs100 =     Dt*v2(1,1);
    const double clhs101 =     -clhs0*clhs99 - clhs100*clhs41 + clhs55*clhs97 + clhs69*clhs98;
    const double clhs102 =     clhs55*clhs95 + clhs69*clhs96;
    const double clhs103 =     clhs55*clhs87 + clhs69*clhs89;
    const double clhs104 =     clhs88*(clhs101*(clhs56*clhs95 + clhs70*clhs96) + clhs102*(-clhs1*clhs99 - clhs100*clhs42 + clhs56*clhs97 + clhs70*clhs98) - clhs103*(clhs0 + clhs1*clhs92 + clhs42*clhs93 - clhs56*clhs90 - clhs70*clhs91) + clhs94*(clhs56*clhs87 + clhs70*clhs89));
    const double clhs105 =     DeltaGap[5]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs106 =     clhs88*(clhs101*(clhs57*clhs95 + clhs71*clhs96) - clhs102*(clhs0 + clhs100*clhs43 - clhs57*clhs97 + clhs7*clhs99 - clhs71*clhs98) + clhs103*(-clhs43*clhs93 + clhs57*clhs90 - clhs7*clhs92 + clhs71*clhs91) + clhs94*(clhs57*clhs87 + clhs71*clhs89));
    const double clhs107 =     DeltaGap[6]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs108 =     clhs88*(clhs101*(clhs58*clhs95 + clhs72*clhs96) + clhs102*(-clhs100*clhs44 + clhs58*clhs97 + clhs72*clhs98 - clhs8*clhs99) - clhs103*(clhs41 + clhs44*clhs93 - clhs58*clhs90 - clhs72*clhs91 + clhs8*clhs92) + clhs94*(clhs58*clhs87 + clhs72*clhs89));
    const double clhs109 =     DeltaGap[7]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs110 =     clhs88*(clhs101*(clhs59*clhs95 + clhs73*clhs96) - clhs102*(clhs100*clhs45 + clhs41 - clhs59*clhs97 - clhs73*clhs98 + clhs9*clhs99) + clhs103*(-clhs45*clhs93 + clhs59*clhs90 + clhs73*clhs91 - clhs9*clhs92) + clhs94*(clhs59*clhs87 + clhs73*clhs89));
    const double clhs111 =     clhs2*clhs3*clhs85;
    const double clhs112 =     DeltaGap[0]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs113 =     clhs2*clhs3*clhs84;
    const double clhs114 =     clhs3*clhs84*clhs85;
    const double clhs115 =     clhs2*clhs84*clhs85;
    const double clhs116 =     Deltatangentxis[0](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs117 =     clhs101*clhs102 + clhs103*clhs94;
    const double clhs118 =     clhs117*clhs2*clhs3*clhs88;
    const double clhs119 =     clhs117*clhs3*clhs87*clhs88;
    const double clhs120 =     clhs117*clhs2*clhs87*clhs88;
    const double clhs121 =     Deltatangentxis[0](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs122 =     Deltatangentxis[0](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs123 =     Deltatangentxis[0](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs124 =     clhs101*(clhs122*clhs55 + clhs123*clhs69 + clhs61*clhs95 + clhs75*clhs96) + clhs102*(-clhs100*clhs47 - clhs12*clhs99 + clhs61*clhs97 + clhs75*clhs98) + clhs103*(-clhs12*clhs92 - clhs47*clhs93 + clhs55 + clhs61*clhs90 + clhs75*clhs91) + clhs94*(clhs116*clhs55 + clhs121*clhs69 + clhs61*clhs87 + clhs75*clhs89);
    const double clhs125 =     clhs124*clhs2*clhs3*clhs88;
    const double clhs126 =     DeltaGap[1]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs127 =     Deltatangentxis[1](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs128 =     Deltatangentxis[1](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs129 =     Deltatangentxis[1](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs130 =     Deltatangentxis[1](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs131 =     clhs101*(clhs129*clhs55 + clhs130*clhs69 + clhs63*clhs95 + clhs77*clhs96) + clhs102*(-clhs100*clhs49 - clhs18*clhs99 + clhs55 + clhs63*clhs97 + clhs77*clhs98) + clhs103*(-clhs18*clhs92 - clhs49*clhs93 + clhs63*clhs90 + clhs77*clhs91) + clhs94*(clhs127*clhs55 + clhs128*clhs69 + clhs63*clhs87 + clhs77*clhs89);
    const double clhs132 =     clhs131*clhs2*clhs3*clhs88;
    const double clhs133 =     DeltaGap[2]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs134 =     Deltatangentxis[2](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs135 =     Deltatangentxis[2](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs136 =     Deltatangentxis[2](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs137 =     Deltatangentxis[2](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs138 =     clhs101*(clhs136*clhs55 + clhs137*clhs69 + clhs64*clhs95 + clhs78*clhs96) + clhs102*(-clhs100*clhs50 - clhs23*clhs99 + clhs64*clhs97 + clhs78*clhs98) + clhs103*(-clhs23*clhs92 - clhs50*clhs93 + clhs64*clhs90 + clhs69 + clhs78*clhs91) + clhs94*(clhs134*clhs55 + clhs135*clhs69 + clhs64*clhs87 + clhs78*clhs89);
    const double clhs139 =     clhs138*clhs2*clhs3*clhs88;
    const double clhs140 =     DeltaGap[3]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs141 =     Deltatangentxis[3](0,0); // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs142 =     Deltatangentxis[3](1,0); // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs143 =     Deltatangentxis[3](0,1); // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs144 =     Deltatangentxis[3](1,1); // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs145 =     clhs101*(clhs143*clhs55 + clhs144*clhs69 + clhs65*clhs95 + clhs79*clhs96) + clhs102*(-clhs100*clhs51 - clhs28*clhs99 + clhs65*clhs97 + clhs69 + clhs79*clhs98) + clhs103*(-clhs28*clhs92 - clhs51*clhs93 + clhs65*clhs90 + clhs79*clhs91) + clhs94*(clhs141*clhs55 + clhs142*clhs69 + clhs65*clhs87 + clhs79*clhs89);
    const double clhs146 =     clhs145*clhs2*clhs3*clhs88;
    const double clhs147 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs148 =     clhs147*clhs2*clhs3;
    const double clhs149 =     clhs147*clhs3*clhs85;
    const double clhs150 =     clhs147*clhs2*clhs85;
    const double clhs151 =     clhs117*clhs3*clhs88*clhs95;
    const double clhs152 =     clhs117*clhs2*clhs88*clhs95;
    const double clhs153 =     clhs2*clhs4;
    const double clhs154 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs155 =     clhs2*clhs4*clhs85;
    const double clhs156 =     clhs154*clhs2*clhs4;
    const double clhs157 =     clhs154*clhs4*clhs85;
    const double clhs158 =     clhs154*clhs2*clhs85;
    const double clhs159 =     clhs117*clhs2*clhs4*clhs88;
    const double clhs160 =     clhs117*clhs4*clhs88*clhs89;
    const double clhs161 =     clhs117*clhs2*clhs88*clhs89;
    const double clhs162 =     clhs124*clhs2*clhs4*clhs88;
    const double clhs163 =     clhs131*clhs2*clhs4*clhs88;
    const double clhs164 =     clhs138*clhs2*clhs4*clhs88;
    const double clhs165 =     clhs145*clhs2*clhs4*clhs88;
    const double clhs166 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs167 =     clhs166*clhs2*clhs4;
    const double clhs168 =     clhs166*clhs4*clhs85;
    const double clhs169 =     clhs166*clhs2*clhs85;
    const double clhs170 =     clhs117*clhs4*clhs88*clhs96;
    const double clhs171 =     clhs117*clhs2*clhs88*clhs96;

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
    lhs(8,0)=clhs83*(clhs104*clhs87 - clhs84*clhs86);
    lhs(8,1)=clhs83*(-clhs105*clhs84 + clhs106*clhs87);
    lhs(8,2)=clhs83*(-clhs107*clhs84 + clhs108*clhs87);
    lhs(8,3)=clhs83*(-clhs109*clhs84 + clhs110*clhs87);
    lhs(8,4)=-clhs10*clhs114 + clhs10*clhs119 - clhs111*DeltaNormals[0](0,0) - clhs112*clhs113 - clhs115*clhs14 + clhs116*clhs118 + clhs120*clhs14 + clhs125*tan1slave(0,0);
    lhs(8,5)=-clhs111*DeltaNormals[1](0,0) - clhs113*clhs126 - clhs114*clhs17 - clhs115*clhs19 + clhs118*clhs127 + clhs119*clhs17 + clhs120*clhs19 + clhs132*tan1slave(0,0);
    lhs(8,6)=-clhs111*DeltaNormals[2](0,0) - clhs113*clhs133 - clhs114*clhs22 - clhs115*clhs24 + clhs118*clhs134 + clhs119*clhs22 + clhs120*clhs24 + clhs139*tan1slave(0,0);
    lhs(8,7)=-clhs111*DeltaNormals[3](0,0) - clhs113*clhs140 - clhs114*clhs27 - clhs115*clhs29 + clhs118*clhs141 + clhs119*clhs27 + clhs120*clhs29 + clhs146*tan1slave(0,0);
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=clhs83*(clhs104*clhs95 - clhs147*clhs86);
    lhs(9,1)=clhs83*(-clhs105*clhs147 + clhs106*clhs95);
    lhs(9,2)=clhs83*(-clhs107*clhs147 + clhs108*clhs95);
    lhs(9,3)=clhs83*(-clhs109*clhs147 + clhs110*clhs95);
    lhs(9,4)=-clhs10*clhs149 + clhs10*clhs151 - clhs111*DeltaNormals[0](0,1) - clhs112*clhs148 + clhs118*clhs122 + clhs125*tan1slave(0,1) - clhs14*clhs150 + clhs14*clhs152;
    lhs(9,5)=-clhs111*DeltaNormals[1](0,1) + clhs118*clhs129 - clhs126*clhs148 + clhs132*tan1slave(0,1) - clhs149*clhs17 - clhs150*clhs19 + clhs151*clhs17 + clhs152*clhs19;
    lhs(9,6)=-clhs111*DeltaNormals[2](0,1) + clhs118*clhs136 - clhs133*clhs148 + clhs139*tan1slave(0,1) - clhs149*clhs22 - clhs150*clhs24 + clhs151*clhs22 + clhs152*clhs24;
    lhs(9,7)=-clhs111*DeltaNormals[3](0,1) + clhs118*clhs143 - clhs140*clhs148 + clhs146*tan1slave(0,1) - clhs149*clhs27 - clhs150*clhs29 + clhs151*clhs27 + clhs152*clhs29;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=clhs153*(clhs104*clhs89 - clhs154*clhs86);
    lhs(10,1)=clhs153*(-clhs105*clhs154 + clhs106*clhs89);
    lhs(10,2)=clhs153*(-clhs107*clhs154 + clhs108*clhs89);
    lhs(10,3)=clhs153*(-clhs109*clhs154 + clhs110*clhs89);
    lhs(10,4)=-clhs10*clhs157 + clhs10*clhs160 - clhs112*clhs156 + clhs121*clhs159 - clhs15*clhs158 + clhs15*clhs161 - clhs155*DeltaNormals[0](1,0) + clhs162*tan1slave(1,0);
    lhs(10,5)=-clhs126*clhs156 + clhs128*clhs159 - clhs155*DeltaNormals[1](1,0) - clhs157*clhs17 - clhs158*clhs20 + clhs160*clhs17 + clhs161*clhs20 + clhs163*tan1slave(1,0);
    lhs(10,6)=-clhs133*clhs156 + clhs135*clhs159 - clhs155*DeltaNormals[2](1,0) - clhs157*clhs22 - clhs158*clhs25 + clhs160*clhs22 + clhs161*clhs25 + clhs164*tan1slave(1,0);
    lhs(10,7)=-clhs140*clhs156 + clhs142*clhs159 - clhs155*DeltaNormals[3](1,0) - clhs157*clhs27 - clhs158*clhs30 + clhs160*clhs27 + clhs161*clhs30 + clhs165*tan1slave(1,0);
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=clhs153*(clhs104*clhs96 - clhs166*clhs86);
    lhs(11,1)=clhs153*(-clhs105*clhs166 + clhs106*clhs96);
    lhs(11,2)=clhs153*(-clhs107*clhs166 + clhs108*clhs96);
    lhs(11,3)=clhs153*(-clhs109*clhs166 + clhs110*clhs96);
    lhs(11,4)=-clhs10*clhs168 + clhs10*clhs170 - clhs112*clhs167 + clhs123*clhs159 - clhs15*clhs169 + clhs15*clhs171 - clhs155*DeltaNormals[0](1,1) + clhs162*tan1slave(1,1);
    lhs(11,5)=-clhs126*clhs167 + clhs130*clhs159 - clhs155*DeltaNormals[1](1,1) + clhs163*tan1slave(1,1) - clhs168*clhs17 - clhs169*clhs20 + clhs17*clhs170 + clhs171*clhs20;
    lhs(11,6)=-clhs133*clhs167 + clhs137*clhs159 - clhs155*DeltaNormals[2](1,1) + clhs164*tan1slave(1,1) - clhs168*clhs22 - clhs169*clhs25 + clhs170*clhs22 + clhs171*clhs25;
    lhs(11,7)=-clhs140*clhs167 + clhs144*clhs159 - clhs155*DeltaNormals[3](1,1) + clhs165*tan1slave(1,1) - clhs168*clhs27 - clhs169*clhs30 + clhs170*clhs27 + clhs171*clhs30;
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;

    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,12,12> ComputeGaussPointLHS(
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const double mu, 
        const ContactData<2,2>& rContactData
        )
{
    bounded_matrix<double,12,12> lhs;
    
    const double integration_point_gap = inner_prod(rContactData.Gaps,N1);
    const bounded_matrix<double, 2, 2> normalslave = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave   = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm          = rContactData.LagrangeMultipliers;
    
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
    const double clhs83 =     clhs2*clhs3;
    const double clhs84 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs85 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs86 =     DeltaGap[4]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs87 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs88 =     Ctan[0]; // CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1))
    const double clhs89 =     DeltaCtan[4][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U2(0,0))
    const double clhs90 =     DeltaGap[5]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs91 =     DeltaCtan[5][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U2(0,1))
    const double clhs92 =     DeltaGap[6]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs93 =     DeltaCtan[6][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U2(1,0))
    const double clhs94 =     DeltaGap[7]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs95 =     DeltaCtan[7][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U2(1,1))
    const double clhs96 =     clhs2*clhs3*clhs85;
    const double clhs97 =     DeltaGap[0]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs98 =     clhs2*clhs3*clhs84;
    const double clhs99 =     clhs3*clhs84*clhs85;
    const double clhs100 =     clhs2*clhs84*clhs85;
    const double clhs101 =     clhs2*clhs3*clhs88;
    const double clhs102 =     clhs3*clhs87*clhs88;
    const double clhs103 =     clhs2*clhs87*clhs88;
    const double clhs104 =     DeltaCtan[0][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U1(0,0))
    const double clhs105 =     clhs2*clhs3*clhs87;
    const double clhs106 =     DeltaGap[1]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs107 =     DeltaCtan[1][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U1(0,1))
    const double clhs108 =     DeltaGap[2]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs109 =     DeltaCtan[2][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U1(1,0))
    const double clhs110 =     DeltaGap[3]; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs111 =     DeltaCtan[3][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), U1(1,1))
    const double clhs112 =     DeltaCtan[8][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), LM(0,0))
    const double clhs113 =     DeltaCtan[9][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), LM(0,1))
    const double clhs114 =     DeltaCtan[10][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), LM(1,0))
    const double clhs115 =     DeltaCtan[11][0]; // DERIVATIVE(CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1)), LM(1,1))
    const double clhs116 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs117 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs118 =     clhs116*clhs2*clhs3;
    const double clhs119 =     clhs116*clhs3*clhs85;
    const double clhs120 =     clhs116*clhs2*clhs85;
    const double clhs121 =     clhs117*clhs3*clhs88;
    const double clhs122 =     clhs117*clhs2*clhs88;
    const double clhs123 =     clhs117*clhs2*clhs3;
    const double clhs124 =     clhs2*clhs4;
    const double clhs125 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs126 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs127 =     clhs2*clhs4*clhs85;
    const double clhs128 =     clhs125*clhs2*clhs4;
    const double clhs129 =     clhs125*clhs4*clhs85;
    const double clhs130 =     clhs125*clhs2*clhs85;
    const double clhs131 =     clhs2*clhs4*clhs88;
    const double clhs132 =     clhs126*clhs4*clhs88;
    const double clhs133 =     clhs126*clhs2*clhs88;
    const double clhs134 =     clhs126*clhs2*clhs4;
    const double clhs135 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs136 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs137 =     clhs135*clhs2*clhs4;
    const double clhs138 =     clhs135*clhs4*clhs85;
    const double clhs139 =     clhs135*clhs2*clhs85;
    const double clhs140 =     clhs136*clhs4*clhs88;
    const double clhs141 =     clhs136*clhs2*clhs88;
    const double clhs142 =     clhs136*clhs2*clhs4;

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
    lhs(8,0)=clhs83*(-clhs84*clhs86 + clhs87*clhs89);
    lhs(8,1)=clhs83*(-clhs84*clhs90 + clhs87*clhs91);
    lhs(8,2)=clhs83*(-clhs84*clhs92 + clhs87*clhs93);
    lhs(8,3)=clhs83*(-clhs84*clhs94 + clhs87*clhs95);
    lhs(8,4)=clhs10*clhs102 - clhs10*clhs99 - clhs100*clhs14 + clhs101*Deltatangentxis[0](0,0) + clhs103*clhs14 + clhs104*clhs105 - clhs96*DeltaNormals[0](0,0) - clhs97*clhs98;
    lhs(8,5)=-clhs100*clhs19 + clhs101*Deltatangentxis[1](0,0) + clhs102*clhs17 + clhs103*clhs19 + clhs105*clhs107 - clhs106*clhs98 - clhs17*clhs99 - clhs96*DeltaNormals[1](0,0);
    lhs(8,6)=-clhs100*clhs24 + clhs101*Deltatangentxis[2](0,0) + clhs102*clhs22 + clhs103*clhs24 + clhs105*clhs109 - clhs108*clhs98 - clhs22*clhs99 - clhs96*DeltaNormals[2](0,0);
    lhs(8,7)=-clhs100*clhs29 + clhs101*Deltatangentxis[3](0,0) + clhs102*clhs27 + clhs103*clhs29 + clhs105*clhs111 - clhs110*clhs98 - clhs27*clhs99 - clhs96*DeltaNormals[3](0,0);
    lhs(8,8)=clhs105*clhs112;
    lhs(8,9)=clhs105*clhs113;
    lhs(8,10)=clhs105*clhs114;
    lhs(8,11)=clhs105*clhs115;
    lhs(9,0)=clhs83*(-clhs116*clhs86 + clhs117*clhs89);
    lhs(9,1)=clhs83*(-clhs116*clhs90 + clhs117*clhs91);
    lhs(9,2)=clhs83*(-clhs116*clhs92 + clhs117*clhs93);
    lhs(9,3)=clhs83*(-clhs116*clhs94 + clhs117*clhs95);
    lhs(9,4)=-clhs10*clhs119 + clhs10*clhs121 + clhs101*Deltatangentxis[0](0,1) + clhs104*clhs123 - clhs118*clhs97 - clhs120*clhs14 + clhs122*clhs14 - clhs96*DeltaNormals[0](0,1);
    lhs(9,5)=clhs101*Deltatangentxis[1](0,1) - clhs106*clhs118 + clhs107*clhs123 - clhs119*clhs17 - clhs120*clhs19 + clhs121*clhs17 + clhs122*clhs19 - clhs96*DeltaNormals[1](0,1);
    lhs(9,6)=clhs101*Deltatangentxis[2](0,1) - clhs108*clhs118 + clhs109*clhs123 - clhs119*clhs22 - clhs120*clhs24 + clhs121*clhs22 + clhs122*clhs24 - clhs96*DeltaNormals[2](0,1);
    lhs(9,7)=clhs101*Deltatangentxis[3](0,1) - clhs110*clhs118 + clhs111*clhs123 - clhs119*clhs27 - clhs120*clhs29 + clhs121*clhs27 + clhs122*clhs29 - clhs96*DeltaNormals[3](0,1);
    lhs(9,8)=clhs112*clhs123;
    lhs(9,9)=clhs113*clhs123;
    lhs(9,10)=clhs114*clhs123;
    lhs(9,11)=clhs115*clhs123;
    lhs(10,0)=clhs124*(-clhs125*clhs86 + clhs126*clhs89);
    lhs(10,1)=clhs124*(-clhs125*clhs90 + clhs126*clhs91);
    lhs(10,2)=clhs124*(-clhs125*clhs92 + clhs126*clhs93);
    lhs(10,3)=clhs124*(-clhs125*clhs94 + clhs126*clhs95);
    lhs(10,4)=-clhs10*clhs129 + clhs10*clhs132 + clhs104*clhs134 - clhs127*DeltaNormals[0](1,0) - clhs128*clhs97 - clhs130*clhs15 + clhs131*Deltatangentxis[0](1,0) + clhs133*clhs15;
    lhs(10,5)=-clhs106*clhs128 + clhs107*clhs134 - clhs127*DeltaNormals[1](1,0) - clhs129*clhs17 - clhs130*clhs20 + clhs131*Deltatangentxis[1](1,0) + clhs132*clhs17 + clhs133*clhs20;
    lhs(10,6)=-clhs108*clhs128 + clhs109*clhs134 - clhs127*DeltaNormals[2](1,0) - clhs129*clhs22 - clhs130*clhs25 + clhs131*Deltatangentxis[2](1,0) + clhs132*clhs22 + clhs133*clhs25;
    lhs(10,7)=-clhs110*clhs128 + clhs111*clhs134 - clhs127*DeltaNormals[3](1,0) - clhs129*clhs27 - clhs130*clhs30 + clhs131*Deltatangentxis[3](1,0) + clhs132*clhs27 + clhs133*clhs30;
    lhs(10,8)=clhs112*clhs134;
    lhs(10,9)=clhs113*clhs134;
    lhs(10,10)=clhs114*clhs134;
    lhs(10,11)=clhs115*clhs134;
    lhs(11,0)=clhs124*(-clhs135*clhs86 + clhs136*clhs89);
    lhs(11,1)=clhs124*(-clhs135*clhs90 + clhs136*clhs91);
    lhs(11,2)=clhs124*(-clhs135*clhs92 + clhs136*clhs93);
    lhs(11,3)=clhs124*(-clhs135*clhs94 + clhs136*clhs95);
    lhs(11,4)=-clhs10*clhs138 + clhs10*clhs140 + clhs104*clhs142 - clhs127*DeltaNormals[0](1,1) + clhs131*Deltatangentxis[0](1,1) - clhs137*clhs97 - clhs139*clhs15 + clhs141*clhs15;
    lhs(11,5)=-clhs106*clhs137 + clhs107*clhs142 - clhs127*DeltaNormals[1](1,1) + clhs131*Deltatangentxis[1](1,1) - clhs138*clhs17 - clhs139*clhs20 + clhs140*clhs17 + clhs141*clhs20;
    lhs(11,6)=-clhs108*clhs137 + clhs109*clhs142 - clhs127*DeltaNormals[2](1,1) + clhs131*Deltatangentxis[2](1,1) - clhs138*clhs22 - clhs139*clhs25 + clhs140*clhs22 + clhs141*clhs25;
    lhs(11,7)=-clhs110*clhs137 + clhs111*clhs142 - clhs127*DeltaNormals[3](1,1) + clhs131*Deltatangentxis[3](1,1) - clhs138*clhs27 - clhs139*clhs30 + clhs140*clhs27 + clhs141*clhs30;
    lhs(11,8)=clhs112*clhs142;
    lhs(11,9)=clhs113*clhs142;
    lhs(11,10)=clhs114*clhs142;
    lhs(11,11)=clhs115*clhs142;

    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,12> ComputeGaussPointFrictionlessRHS(
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const ContactData<2,2>& rContactData
        )
{
    array_1d<double,12> rhs;
    
    const double integration_point_gap = inner_prod(rContactData.Gaps,N1);
    const bounded_matrix<double, 2, 2> normalslave = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave   = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm          = rContactData.LagrangeMultipliers;
    
    const double Dt = rContactData.Dt;
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
    const double crhs9 =     crhs1*crhs2;
    const double crhs10 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs11 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs12 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs13 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs14 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs15 =     ((crhs11*crhs7 + crhs12*crhs8)*(-crhs0*(Dt*v2(0,0)) - crhs6*(Dt*v2(1,0)) + crhs7*(Dt*v1(0,0)) + crhs8*(Dt*v1(1,0))) + (crhs13*crhs7 + crhs14*crhs8)*(-crhs0*(Dt*v2(0,1)) - crhs6*(Dt*v2(1,1)) + crhs7*(Dt*v1(0,1)) + crhs8*(Dt*v1(1,1))))/Dt;
    const double crhs16 =     crhs1*crhs3;

    rhs[0]=-crhs0*crhs4;
    rhs[1]=-crhs0*crhs5;
    rhs[2]=-crhs4*crhs6;
    rhs[3]=-crhs5*crhs6;
    rhs[4]=crhs4*crhs7;
    rhs[5]=crhs5*crhs7;
    rhs[6]=crhs4*crhs8;
    rhs[7]=crhs5*crhs8;
    rhs[8]=crhs9*(crhs10*normalslave(0,0) - crhs11*crhs15);
    rhs[9]=crhs9*(crhs10*normalslave(0,1) - crhs13*crhs15);
    rhs[10]=crhs16*(crhs10*normalslave(1,0) - crhs12*crhs15);
    rhs[11]=crhs16*(crhs10*normalslave(1,1) - crhs14*crhs15);

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,12> ComputeGaussPointRHS(
        const array_1d<double,2> N1, 
        const array_1d<double,2> N2, 
        const array_1d<double,2> Phi, 
        const double detJ, 
        const double mu, 
        const ContactData<2,2>& rContactData
        )
{
    array_1d<double,12> rhs;
    
    const double integration_point_gap = inner_prod(rContactData.Gaps,N1);
    const bounded_matrix<double, 2, 2> normalslave = rContactData.Normal_s;
    const bounded_matrix<double, 2, 2> tan1slave   = rContactData.Tangent_xi_s;
    const bounded_matrix<double, 2, 2> lm          = rContactData.LagrangeMultipliers;    
    
    const array_1d<double, 1> Ctan = rContactData.Ctan; 
    
    const double crhs0 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs1 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs2 =     Phi[0]; // PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs3 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs4 =     crhs1*(crhs2*lm(0,0) + crhs3*lm(1,0));
    const double crhs5 =     crhs1*(crhs2*lm(0,1) + crhs3*lm(1,1));
    const double crhs6 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs7 =     N1[0]; // N1[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs8 =     N1[1]; // N1[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs9 =     crhs1*crhs2;
    const double crhs10 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs11 =     Ctan[0]; // CTANG[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1), LM(0,0), LM(0,1), LM(1,0), LM(1,1))
    const double crhs12 =     crhs1*crhs3;

    rhs[0]=-crhs0*crhs4;
    rhs[1]=-crhs0*crhs5;
    rhs[2]=-crhs4*crhs6;
    rhs[3]=-crhs5*crhs6;
    rhs[4]=crhs4*crhs7;
    rhs[5]=crhs5*crhs7;
    rhs[6]=crhs4*crhs8;
    rhs[7]=crhs5*crhs8;
    rhs[8]=// tan1slave(0,0)
crhs9*(crhs10*normalslave(0,0) - crhs11*tan1slave(0,0));
    rhs[9]=// tan1slave(0,1)
crhs9*(crhs10*normalslave(0,1) - crhs11*tan1slave(0,1));
    rhs[10]=// tan1slave(1,0)
crhs12*(crhs10*normalslave(1,0) - crhs11*tan1slave(1,0));
    rhs[11]=// tan1slave(1,1)
crhs12*(crhs10*normalslave(1,1) - crhs11*tan1slave(1,1));

    
    return rhs;
}

private:
};// class Contact2D2N2N
}
#endif /* KRATOS_CONTACT2D2N2N defined */
