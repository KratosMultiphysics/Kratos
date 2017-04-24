// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrándiz
//

// System includes

// External includes

// Project includes
/* Mortar includes */
#include "custom_conditions/ALM_frictional_mortar_contact_condition.h"

/* Utilities */
#include "custom_utilities/contact_utilities.h"

namespace Kratos 
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesPointerType pProperties ) const
{
    return boost::make_shared< AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes> >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return boost::make_shared< AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes> >( NewId, pGeom, pProperties );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes>::~AugmentedLagrangianMethodFrictionalMortarContactCondition( )
{
}

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/


/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 12, 12> AugmentedLagrangianMethodFrictionalMortarContactCondition<2,2>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const double& rPenaltyFactor,
        const double& rScaleFactor,
        const unsigned int& rActiveInactive
        )
{
    bounded_matrix<double,12,12> lhs = ZeroMatrix(12, 12);
    
    // Master segment info
    GeometryType& CurrentMasterElement = mThisMasterElements[rMasterElementIndex]->GetGeometry();

    // Initialize values
    const bounded_matrix<double, 2, 2> u1 = ContactUtilities::GetVariableMatrix<2,2>(this->GetGeometry(), DISPLACEMENT, 0);
    const bounded_matrix<double, 2, 2> u1old = ContactUtilities::GetVariableMatrix<2,2>(this->GetGeometry(), DISPLACEMENT, 1);
    const bounded_matrix<double, 2, 2> u2 = ContactUtilities::GetVariableMatrix<2,2>(CurrentMasterElement, DISPLACEMENT, 0);
    const bounded_matrix<double, 2, 2> u2old = ContactUtilities::GetVariableMatrix<2,2>(CurrentMasterElement, DISPLACEMENT, 1);
    const bounded_matrix<double, 2, 2> X1 = ContactUtilities::GetCoordinates<2,2>(this->GetGeometry(), false);
    const bounded_matrix<double, 2, 2> X2 = ContactUtilities::GetCoordinates<2,2>(CurrentMasterElement, false);
    
    const bounded_matrix<double, 2, 2> lm = ContactUtilities::GetVariableMatrix<2,2>(this->GetGeometry(), VECTOR_LAGRANGE_MULTIPLIER, 0); 
    
    const bounded_matrix<double, 2, 2> normalslave = ContactUtilities::GetVariableMatrix<2,2>(this->GetGeometry(),  NORMAL);
    const bounded_matrix<double, 2, 2> tangentxislave = ContactUtilities::GetVariableMatrix<2,2>(this->GetGeometry(),  TANGENT_XI);
    const bounded_matrix<double, 2, 2> tangentetaslave = ContactUtilities::GetVariableMatrix<2,2>(this->GetGeometry(),  TANGENT_ETA);
    
    // Mortar operators
    const bounded_matrix<double, 2, 2> MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 2, 2> DOperator = rMortarConditionMatrices.DOperator;
    // Mortar operators derivatives
    const array_1d<bounded_matrix<double, 2, 2>, 8> DeltaMOperator = rMortarConditionMatrices.DeltaMOperator;
    const array_1d<bounded_matrix<double, 2, 2>, 8> DeltaDOperator = rMortarConditionMatrices.DeltaDOperator;

    // We get the friction coefficient
    const array_1d<double, 2> mu = GetFrictionCoefficient();

    // NODE 0
    if (this->GetGeometry()[0].Is(ACTIVE) == false ) // INACTIVE
    {
        const double clhs0 =     0.5*std::pow(rScaleFactor, 2.0)/rPenaltyFactor;
        const double clhs1 =     clhs0*(normalslave(0,0)*normalslave(0,1) + tangentetaslave(0,0)*tangentetaslave(0,1) + tangentxislave(0,0)*tangentxislave(0,1));
    
        lhs(8,8)+=clhs0*(std::pow(normalslave(0,0), 2) + std::pow(tangentetaslave(0,0), 2) + std::pow(tangentxislave(0,0), 2));
        lhs(8,9)+=clhs1;
        lhs(9,8)+=clhs1;
        lhs(9,9)+=clhs0*(std::pow(normalslave(0,1), 2) + std::pow(tangentetaslave(0,1), 2) + std::pow(tangentxislave(0,1), 2));
    }
    else if (this->GetGeometry()[0].Is(SLIP) == true ) // ACTIVE-SLIP
    {
        const double clhs0 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs1 =     DeltaMOperator[4](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs2 =     rScaleFactor*(lm(0,0)*normalslave(0,0) + lm(0,1)*normalslave(0,1));
        const double clhs3 =     X1(0,0) + u1(0,0);
        const double clhs4 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs5 =     X1(1,0) + u1(1,0);
        const double clhs6 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs7 =     X2(0,0) + u2(0,0);
        const double clhs8 =     X2(1,0) + u2(1,0);
        const double clhs9 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs10 =     X1(0,1) + u1(0,1);
        const double clhs11 =     X1(1,1) + u1(1,1);
        const double clhs12 =     X2(0,1) + u2(0,1);
        const double clhs13 =     X2(1,1) + u2(1,1);
        const double clhs14 =     rPenaltyFactor*(normalslave(0,0)*(-clhs0*clhs7 + clhs3*clhs4 + clhs5*clhs6 - clhs8*clhs9) + normalslave(0,1)*(-clhs0*clhs12 + clhs10*clhs4 + clhs11*clhs6 - clhs13*clhs9));
        const double clhs15 =     clhs14 - clhs2;
        const double clhs16 =     clhs15*normalslave(0,0);
        const double clhs17 =     DeltaDOperator[4](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs18 =     DeltaDOperator[4](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs19 =     DeltaMOperator[4](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs20 =     normalslave(0,1)*(-clhs1*clhs12 + clhs10*clhs17 + clhs11*clhs18 - clhs13*clhs19);
        const double clhs21 =     clhs1*clhs7;
        const double clhs22 =     clhs19*clhs8;
        const double clhs23 =     clhs17*clhs3;
        const double clhs24 =     clhs18*clhs5;
        const double clhs25 =     clhs20 - normalslave(0,0)*(clhs0 + clhs21 + clhs22 - clhs23 - clhs24);
        const double clhs26 =     clhs25*normalslave(0,0)*rPenaltyFactor;
        const double clhs27 =     lm(0,0)*tangentetaslave(0,0) + lm(0,1)*tangentetaslave(0,1);
        const double clhs28 =     lm(0,0)*tangentxislave(0,0) + lm(0,1)*tangentxislave(0,1);
        const double clhs29 =     rScaleFactor*(clhs27*tangentetaslave(0,0) + clhs28*tangentxislave(0,0));
        const double clhs30 =     X1(0,0) + u1old(0,0);
        const double clhs31 =     X1(1,0) + u1old(1,0);
        const double clhs32 =     X2(0,0) + u2old(0,0);
        const double clhs33 =     X2(1,0) + u2old(1,0);
        const double clhs34 =     -clhs0*clhs32 + clhs30*clhs4 + clhs31*clhs6 - clhs33*clhs9;
        const double clhs35 =     X1(0,1) + u1old(0,1);
        const double clhs36 =     X1(1,1) + u1old(1,1);
        const double clhs37 =     X2(0,1) + u2old(0,1);
        const double clhs38 =     X2(1,1) + u2old(1,1);
        const double clhs39 =     -clhs0*clhs37 + clhs35*clhs4 + clhs36*clhs6 - clhs38*clhs9;
        const double clhs40 =     rPenaltyFactor*(clhs34*tangentetaslave(0,0) + clhs39*tangentetaslave(0,1));
        const double clhs41 =     rPenaltyFactor*(clhs34*tangentxislave(0,0) + clhs39*tangentxislave(0,1));
        const double clhs42 =     clhs40 + clhs41;
        const double clhs43 =     clhs29 + clhs42;
        const double clhs44 =     rScaleFactor*(clhs27*tangentetaslave(0,1) + clhs28*tangentxislave(0,1));
        const double clhs45 =     clhs42 + clhs44;
        const double clhs46 =     std::pow(clhs43, 2) + std::pow(clhs45, 2);
        const double clhs47 =     std::pow(clhs46, -1.0L/2.0L);
        const double clhs48 =     lm(1,0)*tangentetaslave(1,0) + lm(1,1)*tangentetaslave(1,1);
        const double clhs49 =     lm(1,0)*tangentxislave(1,0) + lm(1,1)*tangentxislave(1,1);
        const double clhs50 =     rScaleFactor*(clhs48*tangentetaslave(1,0) + clhs49*tangentxislave(1,0));
        const double clhs51 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs52 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs53 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs54 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs55 =     clhs30*clhs51 + clhs31*clhs52 - clhs32*clhs53 - clhs33*clhs54;
        const double clhs56 =     clhs35*clhs51 + clhs36*clhs52 - clhs37*clhs53 - clhs38*clhs54;
        const double clhs57 =     rPenaltyFactor*(clhs55*tangentetaslave(1,0) + clhs56*tangentetaslave(1,1));
        const double clhs58 =     rPenaltyFactor*(clhs55*tangentxislave(1,0) + clhs56*tangentxislave(1,1));
        const double clhs59 =     clhs57 + clhs58;
        const double clhs60 =     clhs50 + clhs59;
        const double clhs61 =     rScaleFactor*(clhs48*tangentetaslave(1,1) + clhs49*tangentxislave(1,1));
        const double clhs62 =     clhs59 + clhs61;
        const double clhs63 =     std::pow(clhs60, 2) + std::pow(clhs62, 2);
        const double clhs64 =     std::pow(clhs63, -1.0L/2.0L);
        const double clhs65 =     clhs15*clhs43*clhs47*clhs64*mu[0];
        const double clhs66 =     -clhs1*clhs32 + clhs17*clhs30 + clhs18*clhs31 - clhs19*clhs33;
        const double clhs67 =     -clhs1*clhs37 + clhs17*clhs35 + clhs18*clhs36 - clhs19*clhs38;
        const double clhs68 =     clhs66*tangentetaslave(0,0) + clhs66*tangentxislave(0,0) + clhs67*tangentetaslave(0,1) + clhs67*tangentxislave(0,1);
        const double clhs69 =     clhs15*clhs47*clhs64*clhs68*mu[0]*rPenaltyFactor;
        const double clhs70 =     clhs0*clhs69;
        const double clhs71 =     clhs25*clhs43*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs72 =     clhs29 + 2*clhs40 + 2*clhs41 + clhs44;
        const double clhs73 =     std::pow(clhs46, -3.0L/2.0L);
        const double clhs74 =     clhs15*clhs43*clhs64*clhs68*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs75 =     DeltaDOperator[4](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs76 =     DeltaDOperator[4](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs77 =     DeltaMOperator[4](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs78 =     DeltaMOperator[4](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs79 =     clhs30*clhs75 + clhs31*clhs76 - clhs32*clhs77 - clhs33*clhs78;
        const double clhs80 =     clhs35*clhs75 + clhs36*clhs76 - clhs37*clhs77 - clhs38*clhs78;
        const double clhs81 =     clhs79*tangentetaslave(1,0) + clhs79*tangentxislave(1,0) + clhs80*tangentetaslave(1,1) + clhs80*tangentxislave(1,1);
        const double clhs82 =     clhs50 + 2*clhs57 + 2*clhs58 + clhs61;
        const double clhs83 =     std::pow(clhs63, -3.0L/2.0L);
        const double clhs84 =     clhs15*clhs43*clhs47*clhs81*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs85 =     DeltaMOperator[5](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs86 =     DeltaDOperator[5](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs87 =     DeltaDOperator[5](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs88 =     DeltaMOperator[5](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs89 =     normalslave(0,0)*(clhs3*clhs86 + clhs5*clhs87 - clhs7*clhs85 - clhs8*clhs88);
        const double clhs90 =     clhs12*clhs85;
        const double clhs91 =     clhs13*clhs88;
        const double clhs92 =     clhs10*clhs86;
        const double clhs93 =     clhs11*clhs87;
        const double clhs94 =     clhs89 - normalslave(0,1)*(clhs0 + clhs90 + clhs91 - clhs92 - clhs93);
        const double clhs95 =     clhs94*normalslave(0,0)*rPenaltyFactor;
        const double clhs96 =     clhs30*clhs86 + clhs31*clhs87 - clhs32*clhs85 - clhs33*clhs88;
        const double clhs97 =     clhs35*clhs86 + clhs36*clhs87 - clhs37*clhs85 - clhs38*clhs88;
        const double clhs98 =     clhs96*tangentetaslave(0,0) + clhs96*tangentxislave(0,0) + clhs97*tangentetaslave(0,1) + clhs97*tangentxislave(0,1);
        const double clhs99 =     clhs15*clhs47*clhs64*clhs98*mu[0]*rPenaltyFactor;
        const double clhs100 =     clhs0*clhs99;
        const double clhs101 =     clhs43*clhs47*clhs64*clhs94*mu[0]*rPenaltyFactor;
        const double clhs102 =     clhs15*clhs43*clhs64*clhs72*clhs73*clhs98*mu[0]*rPenaltyFactor;
        const double clhs103 =     DeltaDOperator[5](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs104 =     DeltaDOperator[5](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs105 =     DeltaMOperator[5](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs106 =     DeltaMOperator[5](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs107 =     clhs103*clhs30 + clhs104*clhs31 - clhs105*clhs32 - clhs106*clhs33;
        const double clhs108 =     clhs103*clhs35 + clhs104*clhs36 - clhs105*clhs37 - clhs106*clhs38;
        const double clhs109 =     clhs107*tangentetaslave(1,0) + clhs107*tangentxislave(1,0) + clhs108*tangentetaslave(1,1) + clhs108*tangentxislave(1,1);
        const double clhs110 =     clhs109*clhs15*clhs43*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs111 =     DeltaMOperator[6](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs112 =     DeltaDOperator[6](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs113 =     DeltaDOperator[6](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs114 =     DeltaMOperator[6](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs115 =     normalslave(0,1)*(clhs10*clhs112 + clhs11*clhs113 - clhs111*clhs12 - clhs114*clhs13);
        const double clhs116 =     clhs111*clhs7;
        const double clhs117 =     clhs114*clhs8;
        const double clhs118 =     clhs112*clhs3;
        const double clhs119 =     clhs113*clhs5;
        const double clhs120 =     clhs115 - normalslave(0,0)*(clhs116 + clhs117 - clhs118 - clhs119 + clhs9);
        const double clhs121 =     clhs120*normalslave(0,0)*rPenaltyFactor;
        const double clhs122 =     -clhs111*clhs32 + clhs112*clhs30 + clhs113*clhs31 - clhs114*clhs33;
        const double clhs123 =     -clhs111*clhs37 + clhs112*clhs35 + clhs113*clhs36 - clhs114*clhs38;
        const double clhs124 =     clhs122*tangentetaslave(0,0) + clhs122*tangentxislave(0,0) + clhs123*tangentetaslave(0,1) + clhs123*tangentxislave(0,1);
        const double clhs125 =     clhs124*clhs15*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs126 =     clhs0*clhs125;
        const double clhs127 =     clhs120*clhs43*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs128 =     clhs124*clhs15*clhs43*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs129 =     DeltaDOperator[6](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs130 =     DeltaDOperator[6](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs131 =     DeltaMOperator[6](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs132 =     DeltaMOperator[6](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs133 =     clhs129*clhs30 + clhs130*clhs31 - clhs131*clhs32 - clhs132*clhs33;
        const double clhs134 =     clhs129*clhs35 + clhs130*clhs36 - clhs131*clhs37 - clhs132*clhs38;
        const double clhs135 =     clhs133*tangentetaslave(1,0) + clhs133*tangentxislave(1,0) + clhs134*tangentetaslave(1,1) + clhs134*tangentxislave(1,1);
        const double clhs136 =     clhs135*clhs15*clhs43*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs137 =     DeltaMOperator[7](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs138 =     DeltaDOperator[7](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs139 =     DeltaDOperator[7](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs140 =     DeltaMOperator[7](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs141 =     normalslave(0,0)*(-clhs137*clhs7 + clhs138*clhs3 + clhs139*clhs5 - clhs140*clhs8);
        const double clhs142 =     clhs12*clhs137;
        const double clhs143 =     clhs13*clhs140;
        const double clhs144 =     clhs10*clhs138;
        const double clhs145 =     clhs11*clhs139;
        const double clhs146 =     clhs141 - normalslave(0,1)*(clhs142 + clhs143 - clhs144 - clhs145 + clhs9);
        const double clhs147 =     clhs146*normalslave(0,0)*rPenaltyFactor;
        const double clhs148 =     -clhs137*clhs32 + clhs138*clhs30 + clhs139*clhs31 - clhs140*clhs33;
        const double clhs149 =     -clhs137*clhs37 + clhs138*clhs35 + clhs139*clhs36 - clhs140*clhs38;
        const double clhs150 =     clhs148*tangentetaslave(0,0) + clhs148*tangentxislave(0,0) + clhs149*tangentetaslave(0,1) + clhs149*tangentxislave(0,1);
        const double clhs151 =     clhs15*clhs150*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs152 =     clhs0*clhs151;
        const double clhs153 =     clhs146*clhs43*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs154 =     clhs15*clhs150*clhs43*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs155 =     DeltaDOperator[7](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs156 =     DeltaDOperator[7](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs157 =     DeltaMOperator[7](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs158 =     DeltaMOperator[7](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs159 =     clhs155*clhs30 + clhs156*clhs31 - clhs157*clhs32 - clhs158*clhs33;
        const double clhs160 =     clhs155*clhs35 + clhs156*clhs36 - clhs157*clhs37 - clhs158*clhs38;
        const double clhs161 =     clhs159*tangentetaslave(1,0) + clhs159*tangentxislave(1,0) + clhs160*tangentetaslave(1,1) + clhs160*tangentxislave(1,1);
        const double clhs162 =     clhs15*clhs161*clhs43*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs163 =     DeltaMOperator[0](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs164 =     DeltaDOperator[0](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs165 =     DeltaDOperator[0](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs166 =     DeltaMOperator[0](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs167 =     normalslave(0,0)*(-clhs163*clhs7 + clhs164*clhs3 + clhs165*clhs5 - clhs166*clhs8 + clhs4) + normalslave(0,1)*(clhs10*clhs164 + clhs11*clhs165 - clhs12*clhs163 - clhs13*clhs166);
        const double clhs168 =     clhs167*normalslave(0,0)*rPenaltyFactor;
        const double clhs169 =     -clhs163*clhs32 + clhs164*clhs30 + clhs165*clhs31 - clhs166*clhs33;
        const double clhs170 =     -clhs163*clhs37 + clhs164*clhs35 + clhs165*clhs36 - clhs166*clhs38;
        const double clhs171 =     clhs169*tangentetaslave(0,0) + clhs169*tangentxislave(0,0) + clhs170*tangentetaslave(0,1) + clhs170*tangentxislave(0,1);
        const double clhs172 =     clhs15*clhs171*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs173 =     clhs0*clhs172;
        const double clhs174 =     clhs167*clhs43*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs175 =     clhs15*clhs171*clhs43*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs176 =     DeltaDOperator[0](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs177 =     DeltaDOperator[0](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs178 =     DeltaMOperator[0](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs179 =     DeltaMOperator[0](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs180 =     clhs176*clhs30 + clhs177*clhs31 - clhs178*clhs32 - clhs179*clhs33;
        const double clhs181 =     clhs176*clhs35 + clhs177*clhs36 - clhs178*clhs37 - clhs179*clhs38;
        const double clhs182 =     clhs180*tangentetaslave(1,0) + clhs180*tangentxislave(1,0) + clhs181*tangentetaslave(1,1) + clhs181*tangentxislave(1,1);
        const double clhs183 =     clhs15*clhs182*clhs43*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs184 =     DeltaMOperator[1](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs185 =     DeltaDOperator[1](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs186 =     DeltaDOperator[1](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs187 =     DeltaMOperator[1](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs188 =     normalslave(0,0)*(-clhs184*clhs7 + clhs185*clhs3 + clhs186*clhs5 - clhs187*clhs8) + normalslave(0,1)*(clhs10*clhs185 + clhs11*clhs186 - clhs12*clhs184 - clhs13*clhs187 + clhs4);
        const double clhs189 =     clhs188*normalslave(0,0)*rPenaltyFactor;
        const double clhs190 =     -clhs184*clhs32 + clhs185*clhs30 + clhs186*clhs31 - clhs187*clhs33;
        const double clhs191 =     -clhs184*clhs37 + clhs185*clhs35 + clhs186*clhs36 - clhs187*clhs38;
        const double clhs192 =     clhs190*tangentetaslave(0,0) + clhs190*tangentxislave(0,0) + clhs191*tangentetaslave(0,1) + clhs191*tangentxislave(0,1);
        const double clhs193 =     clhs15*clhs192*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs194 =     clhs0*clhs193;
        const double clhs195 =     clhs188*clhs43*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs196 =     clhs15*clhs192*clhs43*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs197 =     DeltaDOperator[1](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs198 =     DeltaDOperator[1](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs199 =     DeltaMOperator[1](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs200 =     DeltaMOperator[1](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs201 =     clhs197*clhs30 + clhs198*clhs31 - clhs199*clhs32 - clhs200*clhs33;
        const double clhs202 =     clhs197*clhs35 + clhs198*clhs36 - clhs199*clhs37 - clhs200*clhs38;
        const double clhs203 =     clhs201*tangentetaslave(1,0) + clhs201*tangentxislave(1,0) + clhs202*tangentetaslave(1,1) + clhs202*tangentxislave(1,1);
        const double clhs204 =     clhs15*clhs203*clhs43*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs205 =     DeltaMOperator[2](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs206 =     DeltaDOperator[2](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs207 =     DeltaDOperator[2](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs208 =     DeltaMOperator[2](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs209 =     normalslave(0,0)*(-clhs205*clhs7 + clhs206*clhs3 + clhs207*clhs5 - clhs208*clhs8 + clhs6) + normalslave(0,1)*(clhs10*clhs206 + clhs11*clhs207 - clhs12*clhs205 - clhs13*clhs208);
        const double clhs210 =     clhs209*normalslave(0,0)*rPenaltyFactor;
        const double clhs211 =     -clhs205*clhs32 + clhs206*clhs30 + clhs207*clhs31 - clhs208*clhs33;
        const double clhs212 =     -clhs205*clhs37 + clhs206*clhs35 + clhs207*clhs36 - clhs208*clhs38;
        const double clhs213 =     clhs211*tangentetaslave(0,0) + clhs211*tangentxislave(0,0) + clhs212*tangentetaslave(0,1) + clhs212*tangentxislave(0,1);
        const double clhs214 =     clhs15*clhs213*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs215 =     clhs0*clhs214;
        const double clhs216 =     clhs209*clhs43*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs217 =     clhs15*clhs213*clhs43*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs218 =     DeltaDOperator[2](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs219 =     DeltaDOperator[2](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs220 =     DeltaMOperator[2](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs221 =     DeltaMOperator[2](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs222 =     clhs218*clhs30 + clhs219*clhs31 - clhs220*clhs32 - clhs221*clhs33;
        const double clhs223 =     clhs218*clhs35 + clhs219*clhs36 - clhs220*clhs37 - clhs221*clhs38;
        const double clhs224 =     clhs222*tangentetaslave(1,0) + clhs222*tangentxislave(1,0) + clhs223*tangentetaslave(1,1) + clhs223*tangentxislave(1,1);
        const double clhs225 =     clhs15*clhs224*clhs43*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs226 =     DeltaMOperator[3](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs227 =     DeltaDOperator[3](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs228 =     DeltaDOperator[3](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs229 =     DeltaMOperator[3](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs230 =     normalslave(0,0)*(-clhs226*clhs7 + clhs227*clhs3 + clhs228*clhs5 - clhs229*clhs8) + normalslave(0,1)*(clhs10*clhs227 + clhs11*clhs228 - clhs12*clhs226 - clhs13*clhs229 + clhs6);
        const double clhs231 =     clhs230*normalslave(0,0)*rPenaltyFactor;
        const double clhs232 =     -clhs226*clhs32 + clhs227*clhs30 + clhs228*clhs31 - clhs229*clhs33;
        const double clhs233 =     -clhs226*clhs37 + clhs227*clhs35 + clhs228*clhs36 - clhs229*clhs38;
        const double clhs234 =     clhs232*tangentetaslave(0,0) + clhs232*tangentxislave(0,0) + clhs233*tangentetaslave(0,1) + clhs233*tangentxislave(0,1);
        const double clhs235 =     clhs15*clhs234*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs236 =     clhs0*clhs235;
        const double clhs237 =     clhs230*clhs43*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs238 =     clhs15*clhs234*clhs43*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs239 =     DeltaDOperator[3](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs240 =     DeltaDOperator[3](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs241 =     DeltaMOperator[3](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs242 =     DeltaMOperator[3](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs243 =     clhs239*clhs30 + clhs240*clhs31 - clhs241*clhs32 - clhs242*clhs33;
        const double clhs244 =     clhs239*clhs35 + clhs240*clhs36 - clhs241*clhs37 - clhs242*clhs38;
        const double clhs245 =     clhs243*tangentetaslave(1,0) + clhs243*tangentxislave(1,0) + clhs244*tangentetaslave(1,1) + clhs244*tangentxislave(1,1);
        const double clhs246 =     clhs15*clhs245*clhs43*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs247 =     std::pow(normalslave(0,0), 2);
        const double clhs248 =     -clhs14 + clhs2;
        const double clhs249 =     std::pow(tangentetaslave(0,0), 2) + std::pow(tangentxislave(0,0), 2);
        const double clhs250 =     clhs249*clhs47*clhs64*mu[0];
        const double clhs251 =     clhs47*clhs64*mu[0]*normalslave(0,0);
        const double clhs252 =     clhs251*clhs43;
        const double clhs253 =     tangentetaslave(0,0)*tangentetaslave(0,1);
        const double clhs254 =     tangentxislave(0,0)*tangentxislave(0,1);
        const double clhs255 =     clhs253 + clhs254;
        const double clhs256 =     clhs249*clhs43 + clhs255*clhs45;
        const double clhs257 =     clhs248*clhs256*clhs64*clhs73*mu[0];
        const double clhs258 =     rScaleFactor*(clhs247 - clhs248*clhs250 - clhs252 + clhs257*clhs43);
        const double clhs259 =     normalslave(0,0)*normalslave(0,1);
        const double clhs260 =     clhs255*clhs47*clhs64*mu[0];
        const double clhs261 =     -clhs248*clhs260 + clhs259;
        const double clhs262 =     clhs47*clhs64*mu[0]*normalslave(0,1);
        const double clhs263 =     clhs262*clhs43;
        const double clhs264 =     std::pow(tangentetaslave(0,1), 2) + std::pow(tangentxislave(0,1), 2);
        const double clhs265 =     clhs255*clhs43 + clhs264*clhs45;
        const double clhs266 =     clhs248*clhs265*clhs64*clhs73*mu[0];
        const double clhs267 =     rScaleFactor*(clhs261 - clhs263 + clhs266*clhs43);
        const double clhs268 =     tangentetaslave(1,0)*tangentetaslave(1,1) + tangentxislave(1,0)*tangentxislave(1,1);
        const double clhs269 =     clhs268*clhs62 + clhs60*(std::pow(tangentetaslave(1,0), 2) + std::pow(tangentxislave(1,0), 2));
        const double clhs270 =     clhs15*clhs269*clhs43*clhs47*clhs83*mu[0]*rScaleFactor;
        const double clhs271 =     clhs268*clhs60 + clhs62*(std::pow(tangentetaslave(1,1), 2) + std::pow(tangentxislave(1,1), 2));
        const double clhs272 =     clhs15*clhs271*clhs43*clhs47*clhs83*mu[0]*rScaleFactor;
        const double clhs273 =     clhs15*normalslave(0,1);
        const double clhs274 =     clhs25*normalslave(0,1)*rPenaltyFactor;
        const double clhs275 =     clhs15*clhs45*clhs47*clhs64*mu[0];
        const double clhs276 =     clhs25*clhs45*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs277 =     clhs15*clhs45*clhs64*clhs68*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs278 =     clhs15*clhs45*clhs47*clhs81*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs279 =     clhs94*normalslave(0,1)*rPenaltyFactor;
        const double clhs280 =     clhs45*clhs47*clhs64*clhs94*mu[0]*rPenaltyFactor;
        const double clhs281 =     clhs15*clhs45*clhs64*clhs72*clhs73*clhs98*mu[0]*rPenaltyFactor;
        const double clhs282 =     clhs109*clhs15*clhs45*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs283 =     clhs120*normalslave(0,1)*rPenaltyFactor;
        const double clhs284 =     clhs120*clhs45*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs285 =     clhs124*clhs15*clhs45*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs286 =     clhs135*clhs15*clhs45*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs287 =     clhs146*normalslave(0,1)*rPenaltyFactor;
        const double clhs288 =     clhs146*clhs45*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs289 =     clhs15*clhs150*clhs45*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs290 =     clhs15*clhs161*clhs45*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs291 =     clhs167*normalslave(0,1)*rPenaltyFactor;
        const double clhs292 =     clhs167*clhs45*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs293 =     clhs15*clhs171*clhs45*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs294 =     clhs15*clhs182*clhs45*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs295 =     clhs188*normalslave(0,1)*rPenaltyFactor;
        const double clhs296 =     clhs188*clhs45*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs297 =     clhs15*clhs192*clhs45*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs298 =     clhs15*clhs203*clhs45*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs299 =     clhs209*normalslave(0,1)*rPenaltyFactor;
        const double clhs300 =     clhs209*clhs45*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs301 =     clhs15*clhs213*clhs45*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs302 =     clhs15*clhs224*clhs45*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs303 =     clhs230*normalslave(0,1)*rPenaltyFactor;
        const double clhs304 =     clhs230*clhs45*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs305 =     clhs15*clhs234*clhs45*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs306 =     clhs15*clhs245*clhs45*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs307 =     clhs251*clhs45;
        const double clhs308 =     rScaleFactor*(clhs257*clhs45 + clhs261 - clhs307);
        const double clhs309 =     std::pow(normalslave(0,1), 2);
        const double clhs310 =     clhs264*clhs47*clhs64*mu[0];
        const double clhs311 =     clhs262*clhs45;
        const double clhs312 =     rScaleFactor*(-clhs248*clhs310 + clhs266*clhs45 + clhs309 - clhs311);
        const double clhs313 =     clhs15*clhs269*clhs45*clhs47*clhs83*mu[0]*rScaleFactor;
        const double clhs314 =     clhs15*clhs271*clhs45*clhs47*clhs83*mu[0]*rScaleFactor;
        const double clhs315 =     clhs69*clhs9;
        const double clhs316 =     clhs9*clhs99;
        const double clhs317 =     clhs125*clhs9;
        const double clhs318 =     clhs151*clhs9;
        const double clhs319 =     clhs172*clhs9;
        const double clhs320 =     clhs193*clhs9;
        const double clhs321 =     clhs214*clhs9;
        const double clhs322 =     clhs235*clhs9;
        const double clhs323 =     clhs248*normalslave(0,0);
        const double clhs324 =     -clhs0;
        const double clhs325 =     clhs20 + normalslave(0,0)*(-clhs21 - clhs22 + clhs23 + clhs24 + clhs324);
        const double clhs326 =     clhs325*normalslave(0,0)*rPenaltyFactor;
        const double clhs327 =     clhs248*clhs43*clhs47*clhs64*mu[0];
        const double clhs328 =     clhs248*clhs47*clhs64*clhs68*mu[0]*rPenaltyFactor;
        const double clhs329 =     clhs328*clhs4;
        const double clhs330 =     clhs325*clhs43*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs331 =     clhs248*clhs43*clhs64*clhs68*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs332 =     clhs248*clhs43*clhs47*clhs81*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs333 =     clhs89 + normalslave(0,1)*(clhs324 - clhs90 - clhs91 + clhs92 + clhs93);
        const double clhs334 =     clhs333*normalslave(0,0)*rPenaltyFactor;
        const double clhs335 =     clhs248*clhs47*clhs64*clhs98*mu[0]*rPenaltyFactor;
        const double clhs336 =     clhs335*clhs4;
        const double clhs337 =     clhs333*clhs43*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs338 =     clhs248*clhs43*clhs64*clhs72*clhs73*clhs98*mu[0]*rPenaltyFactor;
        const double clhs339 =     clhs109*clhs248*clhs43*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs340 =     -clhs9;
        const double clhs341 =     clhs115 + normalslave(0,0)*(-clhs116 - clhs117 + clhs118 + clhs119 + clhs340);
        const double clhs342 =     clhs341*normalslave(0,0)*rPenaltyFactor;
        const double clhs343 =     clhs124*clhs248*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs344 =     clhs343*clhs4;
        const double clhs345 =     clhs341*clhs43*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs346 =     clhs124*clhs248*clhs43*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs347 =     clhs135*clhs248*clhs43*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs348 =     clhs141 + normalslave(0,1)*(-clhs142 - clhs143 + clhs144 + clhs145 + clhs340);
        const double clhs349 =     clhs348*normalslave(0,0)*rPenaltyFactor;
        const double clhs350 =     clhs150*clhs248*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs351 =     clhs350*clhs4;
        const double clhs352 =     clhs348*clhs43*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs353 =     clhs150*clhs248*clhs43*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs354 =     clhs161*clhs248*clhs43*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs355 =     clhs171*clhs248*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs356 =     clhs355*clhs4;
        const double clhs357 =     clhs171*clhs248*clhs43*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs358 =     clhs182*clhs248*clhs43*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs359 =     clhs192*clhs248*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs360 =     clhs359*clhs4;
        const double clhs361 =     clhs192*clhs248*clhs43*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs362 =     clhs203*clhs248*clhs43*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs363 =     clhs213*clhs248*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs364 =     clhs363*clhs4;
        const double clhs365 =     clhs213*clhs248*clhs43*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs366 =     clhs224*clhs248*clhs43*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs367 =     clhs234*clhs248*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs368 =     clhs367*clhs4;
        const double clhs369 =     clhs234*clhs248*clhs43*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs370 =     clhs245*clhs248*clhs43*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs371 =     clhs15*clhs256*clhs64*clhs73*mu[0];
        const double clhs372 =     -clhs15*clhs250 + clhs252 + clhs371*clhs43;
        const double clhs373 =     rScaleFactor*(-clhs247 + clhs372);
        const double clhs374 =     -clhs15*clhs260;
        const double clhs375 =     -clhs259 + clhs374;
        const double clhs376 =     clhs15*clhs265*clhs64*clhs73*mu[0];
        const double clhs377 =     clhs263 + clhs376*clhs43;
        const double clhs378 =     rScaleFactor*(clhs375 + clhs377);
        const double clhs379 =     clhs248*normalslave(0,1);
        const double clhs380 =     clhs325*normalslave(0,1)*rPenaltyFactor;
        const double clhs381 =     clhs248*clhs45*clhs47*clhs64*mu[0];
        const double clhs382 =     clhs325*clhs45*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs383 =     clhs248*clhs45*clhs64*clhs68*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs384 =     clhs248*clhs45*clhs47*clhs81*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs385 =     clhs333*normalslave(0,1)*rPenaltyFactor;
        const double clhs386 =     clhs333*clhs45*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs387 =     clhs248*clhs45*clhs64*clhs72*clhs73*clhs98*mu[0]*rPenaltyFactor;
        const double clhs388 =     clhs109*clhs248*clhs45*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs389 =     clhs341*normalslave(0,1)*rPenaltyFactor;
        const double clhs390 =     clhs341*clhs45*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs391 =     clhs124*clhs248*clhs45*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs392 =     clhs135*clhs248*clhs45*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs393 =     clhs348*normalslave(0,1)*rPenaltyFactor;
        const double clhs394 =     clhs348*clhs45*clhs47*clhs64*mu[0]*rPenaltyFactor;
        const double clhs395 =     clhs150*clhs248*clhs45*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs396 =     clhs161*clhs248*clhs45*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs397 =     clhs171*clhs248*clhs45*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs398 =     clhs182*clhs248*clhs45*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs399 =     clhs192*clhs248*clhs45*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs400 =     clhs203*clhs248*clhs45*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs401 =     clhs213*clhs248*clhs45*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs402 =     clhs224*clhs248*clhs45*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs403 =     clhs234*clhs248*clhs45*clhs64*clhs72*clhs73*mu[0]*rPenaltyFactor;
        const double clhs404 =     clhs245*clhs248*clhs45*clhs47*clhs82*clhs83*mu[0]*rPenaltyFactor;
        const double clhs405 =     clhs307 + clhs371*clhs45;
        const double clhs406 =     rScaleFactor*(clhs375 + clhs405);
        const double clhs407 =     -clhs15*clhs310 + clhs311 + clhs376*clhs45;
        const double clhs408 =     rScaleFactor*(-clhs309 + clhs407);
        const double clhs409 =     clhs328*clhs6;
        const double clhs410 =     clhs335*clhs6;
        const double clhs411 =     clhs343*clhs6;
        const double clhs412 =     clhs350*clhs6;
        const double clhs413 =     clhs355*clhs6;
        const double clhs414 =     clhs359*clhs6;
        const double clhs415 =     clhs363*clhs6;
        const double clhs416 =     clhs367*clhs6;
        const double clhs417 =     0.5*clhs249*clhs47*clhs64*mu[0];
        const double clhs418 =     clhs248*clhs68;
        const double clhs419 =     -clhs418;
        const double clhs420 =     1.0/clhs46;
        const double clhs421 =     clhs420*clhs43*clhs72;
        const double clhs422 =     1.0/clhs63;
        const double clhs423 =     clhs248*clhs422*clhs43*clhs82;
        const double clhs424 =     clhs325*clhs43 + clhs418*clhs421 + clhs419 + clhs423*clhs81;
        const double clhs425 =     0.5*clhs255*clhs47*clhs64*mu[0];
        const double clhs426 =     clhs420*clhs45*clhs72;
        const double clhs427 =     clhs248*clhs422*clhs45*clhs82;
        const double clhs428 =     clhs325*clhs45 + clhs418*clhs426 + clhs419 + clhs427*clhs81;
        const double clhs429 =     clhs248*clhs98;
        const double clhs430 =     -clhs429;
        const double clhs431 =     clhs109*clhs423 + clhs333*clhs43 + clhs421*clhs429 + clhs430;
        const double clhs432 =     clhs109*clhs427 + clhs333*clhs45 + clhs426*clhs429 + clhs430;
        const double clhs433 =     clhs124*clhs248;
        const double clhs434 =     -clhs433;
        const double clhs435 =     clhs135*clhs423 + clhs341*clhs43 + clhs421*clhs433 + clhs434;
        const double clhs436 =     clhs135*clhs427 + clhs341*clhs45 + clhs426*clhs433 + clhs434;
        const double clhs437 =     clhs150*clhs248;
        const double clhs438 =     -clhs437;
        const double clhs439 =     clhs161*clhs423 + clhs348*clhs43 + clhs421*clhs437 + clhs438;
        const double clhs440 =     clhs161*clhs427 + clhs348*clhs45 + clhs426*clhs437 + clhs438;
        const double clhs441 =     clhs171*clhs248;
        const double clhs442 =     -clhs441;
        const double clhs443 =     clhs167*clhs43 + clhs182*clhs423 + clhs421*clhs441 + clhs442;
        const double clhs444 =     clhs167*clhs45 + clhs182*clhs427 + clhs426*clhs441 + clhs442;
        const double clhs445 =     clhs192*clhs248;
        const double clhs446 =     -clhs445;
        const double clhs447 =     clhs188*clhs43 + clhs203*clhs423 + clhs421*clhs445 + clhs446;
        const double clhs448 =     clhs188*clhs45 + clhs203*clhs427 + clhs426*clhs445 + clhs446;
        const double clhs449 =     clhs213*clhs248;
        const double clhs450 =     -clhs449;
        const double clhs451 =     clhs209*clhs43 + clhs224*clhs423 + clhs421*clhs449 + clhs450;
        const double clhs452 =     clhs209*clhs45 + clhs224*clhs427 + clhs426*clhs449 + clhs450;
        const double clhs453 =     clhs234*clhs248;
        const double clhs454 =     -clhs453;
        const double clhs455 =     clhs230*clhs43 + clhs245*clhs423 + clhs421*clhs453 + clhs454;
        const double clhs456 =     clhs230*clhs45 + clhs245*clhs427 + clhs426*clhs453 + clhs454;
        const double clhs457 =     1.0/rPenaltyFactor;
        const double clhs458 =     std::pow(rScaleFactor, 2);
        const double clhs459 =     0.5*clhs457*clhs458;
        const double clhs460 =     clhs249 + clhs372;
        const double clhs461 =     clhs253 + clhs254 + clhs374;
        const double clhs462 =     clhs405 + clhs461;
        const double clhs463 =     clhs377 + clhs461;
        const double clhs464 =     clhs264 + clhs407;
        const double clhs465 =     0.5*clhs248*clhs256*clhs457*clhs458*clhs47*clhs83*mu[0];
        const double clhs466 =     0.5*clhs264*clhs47*clhs64*mu[0];
        const double clhs467 =     0.5*clhs248*clhs265*clhs457*clhs458*clhs47*clhs83*mu[0];
    
        lhs(0,0)+=-clhs0*clhs26 + clhs0*clhs71 - clhs0*clhs74 - clhs0*clhs84 - clhs1*clhs16 + clhs1*clhs65 + clhs70;
        lhs(0,1)+=clhs0*clhs101 - clhs0*clhs102 - clhs0*clhs110 - clhs0*clhs95 + clhs100 - clhs16*clhs85 + clhs65*clhs85;
        lhs(0,2)+=-clhs0*clhs121 + clhs0*clhs127 - clhs0*clhs128 - clhs0*clhs136 - clhs111*clhs16 + clhs111*clhs65 + clhs126;
        lhs(0,3)+=-clhs0*clhs147 + clhs0*clhs153 - clhs0*clhs154 - clhs0*clhs162 - clhs137*clhs16 + clhs137*clhs65 + clhs152;
        lhs(0,4)+=-clhs0*clhs168 + clhs0*clhs174 - clhs0*clhs175 - clhs0*clhs183 - clhs16*clhs163 + clhs163*clhs65 + clhs173;
        lhs(0,5)+=-clhs0*clhs189 + clhs0*clhs195 - clhs0*clhs196 - clhs0*clhs204 - clhs16*clhs184 + clhs184*clhs65 + clhs194;
        lhs(0,6)+=-clhs0*clhs210 + clhs0*clhs216 - clhs0*clhs217 - clhs0*clhs225 - clhs16*clhs205 + clhs205*clhs65 + clhs215;
        lhs(0,7)+=-clhs0*clhs231 + clhs0*clhs237 - clhs0*clhs238 - clhs0*clhs246 - clhs16*clhs226 + clhs226*clhs65 + clhs236;
        lhs(0,8)+=clhs0*clhs258;
        lhs(0,9)+=clhs0*clhs267;
        lhs(0,10)+=-clhs0*clhs270;
        lhs(0,11)+=-clhs0*clhs272;
        lhs(1,0)+=-clhs0*clhs274 + clhs0*clhs276 - clhs0*clhs277 - clhs0*clhs278 - clhs1*clhs273 + clhs1*clhs275 + clhs70;
        lhs(1,1)+=-clhs0*clhs279 + clhs0*clhs280 - clhs0*clhs281 - clhs0*clhs282 + clhs100 - clhs273*clhs85 + clhs275*clhs85;
        lhs(1,2)+=-clhs0*clhs283 + clhs0*clhs284 - clhs0*clhs285 - clhs0*clhs286 - clhs111*clhs273 + clhs111*clhs275 + clhs126;
        lhs(1,3)+=-clhs0*clhs287 + clhs0*clhs288 - clhs0*clhs289 - clhs0*clhs290 - clhs137*clhs273 + clhs137*clhs275 + clhs152;
        lhs(1,4)+=-clhs0*clhs291 + clhs0*clhs292 - clhs0*clhs293 - clhs0*clhs294 - clhs163*clhs273 + clhs163*clhs275 + clhs173;
        lhs(1,5)+=-clhs0*clhs295 + clhs0*clhs296 - clhs0*clhs297 - clhs0*clhs298 - clhs184*clhs273 + clhs184*clhs275 + clhs194;
        lhs(1,6)+=-clhs0*clhs299 + clhs0*clhs300 - clhs0*clhs301 - clhs0*clhs302 - clhs205*clhs273 + clhs205*clhs275 + clhs215;
        lhs(1,7)+=-clhs0*clhs303 + clhs0*clhs304 - clhs0*clhs305 - clhs0*clhs306 - clhs226*clhs273 + clhs226*clhs275 + clhs236;
        lhs(1,8)+=clhs0*clhs308;
        lhs(1,9)+=clhs0*clhs312;
        lhs(1,10)+=-clhs0*clhs313;
        lhs(1,11)+=-clhs0*clhs314;
        lhs(2,0)+=-clhs16*clhs19 + clhs19*clhs65 - clhs26*clhs9 + clhs315 + clhs71*clhs9 - clhs74*clhs9 - clhs84*clhs9;
        lhs(2,1)+=clhs101*clhs9 - clhs102*clhs9 - clhs110*clhs9 - clhs16*clhs88 + clhs316 + clhs65*clhs88 - clhs9*clhs95;
        lhs(2,2)+=-clhs114*clhs16 + clhs114*clhs65 - clhs121*clhs9 + clhs127*clhs9 - clhs128*clhs9 - clhs136*clhs9 + clhs317;
        lhs(2,3)+=-clhs140*clhs16 + clhs140*clhs65 - clhs147*clhs9 + clhs153*clhs9 - clhs154*clhs9 - clhs162*clhs9 + clhs318;
        lhs(2,4)+=-clhs16*clhs166 + clhs166*clhs65 - clhs168*clhs9 + clhs174*clhs9 - clhs175*clhs9 - clhs183*clhs9 + clhs319;
        lhs(2,5)+=-clhs16*clhs187 + clhs187*clhs65 - clhs189*clhs9 + clhs195*clhs9 - clhs196*clhs9 - clhs204*clhs9 + clhs320;
        lhs(2,6)+=-clhs16*clhs208 + clhs208*clhs65 - clhs210*clhs9 + clhs216*clhs9 - clhs217*clhs9 - clhs225*clhs9 + clhs321;
        lhs(2,7)+=-clhs16*clhs229 + clhs229*clhs65 - clhs231*clhs9 + clhs237*clhs9 - clhs238*clhs9 - clhs246*clhs9 + clhs322;
        lhs(2,8)+=clhs258*clhs9;
        lhs(2,9)+=clhs267*clhs9;
        lhs(2,10)+=-clhs270*clhs9;
        lhs(2,11)+=-clhs272*clhs9;
        lhs(3,0)+=-clhs19*clhs273 + clhs19*clhs275 - clhs274*clhs9 + clhs276*clhs9 - clhs277*clhs9 - clhs278*clhs9 + clhs315;
        lhs(3,1)+=-clhs273*clhs88 + clhs275*clhs88 - clhs279*clhs9 + clhs280*clhs9 - clhs281*clhs9 - clhs282*clhs9 + clhs316;
        lhs(3,2)+=-clhs114*clhs273 + clhs114*clhs275 - clhs283*clhs9 + clhs284*clhs9 - clhs285*clhs9 - clhs286*clhs9 + clhs317;
        lhs(3,3)+=-clhs140*clhs273 + clhs140*clhs275 - clhs287*clhs9 + clhs288*clhs9 - clhs289*clhs9 - clhs290*clhs9 + clhs318;
        lhs(3,4)+=-clhs166*clhs273 + clhs166*clhs275 - clhs291*clhs9 + clhs292*clhs9 - clhs293*clhs9 - clhs294*clhs9 + clhs319;
        lhs(3,5)+=-clhs187*clhs273 + clhs187*clhs275 - clhs295*clhs9 + clhs296*clhs9 - clhs297*clhs9 - clhs298*clhs9 + clhs320;
        lhs(3,6)+=-clhs208*clhs273 + clhs208*clhs275 - clhs299*clhs9 + clhs300*clhs9 - clhs301*clhs9 - clhs302*clhs9 + clhs321;
        lhs(3,7)+=-clhs229*clhs273 + clhs229*clhs275 - clhs303*clhs9 + clhs304*clhs9 - clhs305*clhs9 - clhs306*clhs9 + clhs322;
        lhs(3,8)+=clhs308*clhs9;
        lhs(3,9)+=clhs312*clhs9;
        lhs(3,10)+=-clhs313*clhs9;
        lhs(3,11)+=-clhs314*clhs9;
        lhs(4,0)+=-clhs17*clhs323 + clhs17*clhs327 + clhs326*clhs4 + clhs329 - clhs330*clhs4 - clhs331*clhs4 - clhs332*clhs4;
        lhs(4,1)+=-clhs323*clhs86 + clhs327*clhs86 + clhs334*clhs4 + clhs336 - clhs337*clhs4 - clhs338*clhs4 - clhs339*clhs4;
        lhs(4,2)+=-clhs112*clhs323 + clhs112*clhs327 + clhs342*clhs4 + clhs344 - clhs345*clhs4 - clhs346*clhs4 - clhs347*clhs4;
        lhs(4,3)+=-clhs138*clhs323 + clhs138*clhs327 + clhs349*clhs4 + clhs351 - clhs352*clhs4 - clhs353*clhs4 - clhs354*clhs4;
        lhs(4,4)+=-clhs164*clhs323 + clhs164*clhs327 + clhs168*clhs4 - clhs174*clhs4 + clhs356 - clhs357*clhs4 - clhs358*clhs4;
        lhs(4,5)+=-clhs185*clhs323 + clhs185*clhs327 + clhs189*clhs4 - clhs195*clhs4 + clhs360 - clhs361*clhs4 - clhs362*clhs4;
        lhs(4,6)+=-clhs206*clhs323 + clhs206*clhs327 + clhs210*clhs4 - clhs216*clhs4 + clhs364 - clhs365*clhs4 - clhs366*clhs4;
        lhs(4,7)+=-clhs227*clhs323 + clhs227*clhs327 + clhs231*clhs4 - clhs237*clhs4 + clhs368 - clhs369*clhs4 - clhs370*clhs4;
        lhs(4,8)+=clhs373*clhs4;
        lhs(4,9)+=clhs378*clhs4;
        lhs(4,10)+=clhs270*clhs4;
        lhs(4,11)+=clhs272*clhs4;
        lhs(5,0)+=-clhs17*clhs379 + clhs17*clhs381 + clhs329 + clhs380*clhs4 - clhs382*clhs4 - clhs383*clhs4 - clhs384*clhs4;
        lhs(5,1)+=clhs336 - clhs379*clhs86 + clhs381*clhs86 + clhs385*clhs4 - clhs386*clhs4 - clhs387*clhs4 - clhs388*clhs4;
        lhs(5,2)+=-clhs112*clhs379 + clhs112*clhs381 + clhs344 + clhs389*clhs4 - clhs390*clhs4 - clhs391*clhs4 - clhs392*clhs4;
        lhs(5,3)+=-clhs138*clhs379 + clhs138*clhs381 + clhs351 + clhs393*clhs4 - clhs394*clhs4 - clhs395*clhs4 - clhs396*clhs4;
        lhs(5,4)+=-clhs164*clhs379 + clhs164*clhs381 + clhs291*clhs4 - clhs292*clhs4 + clhs356 - clhs397*clhs4 - clhs398*clhs4;
        lhs(5,5)+=-clhs185*clhs379 + clhs185*clhs381 + clhs295*clhs4 - clhs296*clhs4 + clhs360 - clhs399*clhs4 - clhs4*clhs400;
        lhs(5,6)+=-clhs206*clhs379 + clhs206*clhs381 + clhs299*clhs4 - clhs300*clhs4 + clhs364 - clhs4*clhs401 - clhs4*clhs402;
        lhs(5,7)+=-clhs227*clhs379 + clhs227*clhs381 + clhs303*clhs4 - clhs304*clhs4 + clhs368 - clhs4*clhs403 - clhs4*clhs404;
        lhs(5,8)+=clhs4*clhs406;
        lhs(5,9)+=clhs4*clhs408;
        lhs(5,10)+=clhs313*clhs4;
        lhs(5,11)+=clhs314*clhs4;
        lhs(6,0)+=-clhs18*clhs323 + clhs18*clhs327 + clhs326*clhs6 - clhs330*clhs6 - clhs331*clhs6 - clhs332*clhs6 + clhs409;
        lhs(6,1)+=-clhs323*clhs87 + clhs327*clhs87 + clhs334*clhs6 - clhs337*clhs6 - clhs338*clhs6 - clhs339*clhs6 + clhs410;
        lhs(6,2)+=-clhs113*clhs323 + clhs113*clhs327 + clhs342*clhs6 - clhs345*clhs6 - clhs346*clhs6 - clhs347*clhs6 + clhs411;
        lhs(6,3)+=-clhs139*clhs323 + clhs139*clhs327 + clhs349*clhs6 - clhs352*clhs6 - clhs353*clhs6 - clhs354*clhs6 + clhs412;
        lhs(6,4)+=-clhs165*clhs323 + clhs165*clhs327 + clhs168*clhs6 - clhs174*clhs6 - clhs357*clhs6 - clhs358*clhs6 + clhs413;
        lhs(6,5)+=-clhs186*clhs323 + clhs186*clhs327 + clhs189*clhs6 - clhs195*clhs6 - clhs361*clhs6 - clhs362*clhs6 + clhs414;
        lhs(6,6)+=-clhs207*clhs323 + clhs207*clhs327 + clhs210*clhs6 - clhs216*clhs6 - clhs365*clhs6 - clhs366*clhs6 + clhs415;
        lhs(6,7)+=-clhs228*clhs323 + clhs228*clhs327 + clhs231*clhs6 - clhs237*clhs6 - clhs369*clhs6 - clhs370*clhs6 + clhs416;
        lhs(6,8)+=clhs373*clhs6;
        lhs(6,9)+=clhs378*clhs6;
        lhs(6,10)+=clhs270*clhs6;
        lhs(6,11)+=clhs272*clhs6;
        lhs(7,0)+=-clhs18*clhs379 + clhs18*clhs381 + clhs380*clhs6 - clhs382*clhs6 - clhs383*clhs6 - clhs384*clhs6 + clhs409;
        lhs(7,1)+=-clhs379*clhs87 + clhs381*clhs87 + clhs385*clhs6 - clhs386*clhs6 - clhs387*clhs6 - clhs388*clhs6 + clhs410;
        lhs(7,2)+=-clhs113*clhs379 + clhs113*clhs381 + clhs389*clhs6 - clhs390*clhs6 - clhs391*clhs6 - clhs392*clhs6 + clhs411;
        lhs(7,3)+=-clhs139*clhs379 + clhs139*clhs381 + clhs393*clhs6 - clhs394*clhs6 - clhs395*clhs6 - clhs396*clhs6 + clhs412;
        lhs(7,4)+=-clhs165*clhs379 + clhs165*clhs381 + clhs291*clhs6 - clhs292*clhs6 - clhs397*clhs6 - clhs398*clhs6 + clhs413;
        lhs(7,5)+=-clhs186*clhs379 + clhs186*clhs381 + clhs295*clhs6 - clhs296*clhs6 - clhs399*clhs6 - clhs400*clhs6 + clhs414;
        lhs(7,6)+=-clhs207*clhs379 + clhs207*clhs381 + clhs299*clhs6 - clhs300*clhs6 - clhs401*clhs6 - clhs402*clhs6 + clhs415;
        lhs(7,7)+=-clhs228*clhs379 + clhs228*clhs381 + clhs303*clhs6 - clhs304*clhs6 - clhs403*clhs6 - clhs404*clhs6 + clhs416;
        lhs(7,8)+=clhs406*clhs6;
        lhs(7,9)+=clhs408*clhs6;
        lhs(7,10)+=clhs313*clhs6;
        lhs(7,11)+=clhs314*clhs6;
        lhs(8,0)+=rScaleFactor*(clhs25*normalslave(0,0) - clhs417*clhs424 - clhs425*clhs428);
        lhs(8,1)+=rScaleFactor*(-clhs417*clhs431 - clhs425*clhs432 + clhs94*normalslave(0,0));
        lhs(8,2)+=rScaleFactor*(clhs120*normalslave(0,0) - clhs417*clhs435 - clhs425*clhs436);
        lhs(8,3)+=rScaleFactor*(clhs146*normalslave(0,0) - clhs417*clhs439 - clhs425*clhs440);
        lhs(8,4)+=rScaleFactor*(clhs167*normalslave(0,0) - clhs417*clhs443 - clhs425*clhs444);
        lhs(8,5)+=rScaleFactor*(clhs188*normalslave(0,0) - clhs417*clhs447 - clhs425*clhs448);
        lhs(8,6)+=rScaleFactor*(clhs209*normalslave(0,0) - clhs417*clhs451 - clhs425*clhs452);
        lhs(8,7)+=rScaleFactor*(clhs230*normalslave(0,0) - clhs417*clhs455 - clhs425*clhs456);
        lhs(8,8)+=clhs459*(clhs249*clhs460 + clhs255*clhs462);
        lhs(8,9)+=clhs459*(clhs249*clhs463 + clhs255*clhs464);
        lhs(8,10)+=-clhs269*clhs465;
        lhs(8,11)+=-clhs271*clhs465;
        lhs(9,0)+=rScaleFactor*(clhs25*normalslave(0,1) - clhs424*clhs425 - clhs428*clhs466);
        lhs(9,1)+=rScaleFactor*(-clhs425*clhs431 - clhs432*clhs466 + clhs94*normalslave(0,1));
        lhs(9,2)+=rScaleFactor*(clhs120*normalslave(0,1) - clhs425*clhs435 - clhs436*clhs466);
        lhs(9,3)+=rScaleFactor*(clhs146*normalslave(0,1) - clhs425*clhs439 - clhs440*clhs466);
        lhs(9,4)+=rScaleFactor*(clhs167*normalslave(0,1) - clhs425*clhs443 - clhs444*clhs466);
        lhs(9,5)+=rScaleFactor*(clhs188*normalslave(0,1) - clhs425*clhs447 - clhs448*clhs466);
        lhs(9,6)+=rScaleFactor*(clhs209*normalslave(0,1) - clhs425*clhs451 - clhs452*clhs466);
        lhs(9,7)+=rScaleFactor*(clhs230*normalslave(0,1) - clhs425*clhs455 - clhs456*clhs466);
        lhs(9,8)+=clhs459*(clhs255*clhs460 + clhs264*clhs462);
        lhs(9,9)+=clhs459*(clhs255*clhs463 + clhs264*clhs464);
        lhs(9,10)+=-clhs269*clhs467;
        lhs(9,11)+=-clhs271*clhs467;
    }
    else // ACTIVE-STICK
    {
        const double clhs0 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs1 =     DeltaMOperator[4](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs2 =     X1(0,0) + u1old(0,0);
        const double clhs3 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs4 =     X1(1,0) + u1old(1,0);
        const double clhs5 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs6 =     X2(0,0) + u2old(0,0);
        const double clhs7 =     X2(1,0) + u2old(1,0);
        const double clhs8 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs9 =     -clhs0*clhs6 + clhs2*clhs3 + clhs4*clhs5 - clhs7*clhs8;
        const double clhs10 =     X1(0,1) + u1old(0,1);
        const double clhs11 =     X1(1,1) + u1old(1,1);
        const double clhs12 =     X2(0,1) + u2old(0,1);
        const double clhs13 =     X2(1,1) + u2old(1,1);
        const double clhs14 =     -clhs0*clhs12 + clhs10*clhs3 + clhs11*clhs5 - clhs13*clhs8;
        const double clhs15 =     rPenaltyFactor*(clhs14*tangentetaslave(0,1) + clhs9*tangentetaslave(0,0)) + rScaleFactor*(lm(0,0)*tangentetaslave(0,0) + lm(0,1)*tangentetaslave(0,1));
        const double clhs16 =     clhs15*tangentetaslave(0,0);
        const double clhs17 =     rPenaltyFactor*(clhs14*tangentxislave(0,1) + clhs9*tangentxislave(0,0)) + rScaleFactor*(lm(0,0)*tangentxislave(0,0) + lm(0,1)*tangentxislave(0,1));
        const double clhs18 =     clhs17*tangentxislave(0,0);
        const double clhs19 =     rScaleFactor*(lm(0,0)*normalslave(0,0) + lm(0,1)*normalslave(0,1));
        const double clhs20 =     X1(0,0) + u1(0,0);
        const double clhs21 =     X1(1,0) + u1(1,0);
        const double clhs22 =     X2(0,0) + u2(0,0);
        const double clhs23 =     X2(1,0) + u2(1,0);
        const double clhs24 =     X1(0,1) + u1(0,1);
        const double clhs25 =     X1(1,1) + u1(1,1);
        const double clhs26 =     X2(0,1) + u2(0,1);
        const double clhs27 =     X2(1,1) + u2(1,1);
        const double clhs28 =     rPenaltyFactor*(normalslave(0,0)*(-clhs0*clhs22 + clhs20*clhs3 + clhs21*clhs5 - clhs23*clhs8) + normalslave(0,1)*(-clhs0*clhs26 + clhs24*clhs3 + clhs25*clhs5 - clhs27*clhs8));
        const double clhs29 =     -clhs19 + clhs28;
        const double clhs30 =     clhs29*normalslave(0,0);
        const double clhs31 =     DeltaDOperator[4](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs32 =     DeltaDOperator[4](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs33 =     DeltaMOperator[4](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs34 =     -clhs1*clhs6 + clhs2*clhs31 + clhs32*clhs4 - clhs33*clhs7;
        const double clhs35 =     -clhs1*clhs12 + clhs10*clhs31 + clhs11*clhs32 - clhs13*clhs33;
        const double clhs36 =     clhs34*tangentetaslave(0,0) + clhs35*tangentetaslave(0,1);
        const double clhs37 =     clhs36*rPenaltyFactor*tangentetaslave(0,0);
        const double clhs38 =     clhs34*tangentxislave(0,0) + clhs35*tangentxislave(0,1);
        const double clhs39 =     clhs38*rPenaltyFactor*tangentxislave(0,0);
        const double clhs40 =     normalslave(0,1)*(-clhs1*clhs26 + clhs24*clhs31 + clhs25*clhs32 - clhs27*clhs33);
        const double clhs41 =     clhs1*clhs22;
        const double clhs42 =     clhs23*clhs33;
        const double clhs43 =     clhs20*clhs31;
        const double clhs44 =     clhs21*clhs32;
        const double clhs45 =     clhs40 - normalslave(0,0)*(clhs0 + clhs41 + clhs42 - clhs43 - clhs44);
        const double clhs46 =     clhs45*normalslave(0,0)*rPenaltyFactor;
        const double clhs47 =     DeltaMOperator[5](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs48 =     DeltaDOperator[5](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs49 =     DeltaDOperator[5](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs50 =     DeltaMOperator[5](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs51 =     clhs2*clhs48 + clhs4*clhs49 - clhs47*clhs6 - clhs50*clhs7;
        const double clhs52 =     clhs10*clhs48 + clhs11*clhs49 - clhs12*clhs47 - clhs13*clhs50;
        const double clhs53 =     clhs51*tangentetaslave(0,0) + clhs52*tangentetaslave(0,1);
        const double clhs54 =     clhs53*rPenaltyFactor*tangentetaslave(0,0);
        const double clhs55 =     clhs51*tangentxislave(0,0) + clhs52*tangentxislave(0,1);
        const double clhs56 =     clhs55*rPenaltyFactor*tangentxislave(0,0);
        const double clhs57 =     normalslave(0,0)*(clhs20*clhs48 + clhs21*clhs49 - clhs22*clhs47 - clhs23*clhs50);
        const double clhs58 =     clhs26*clhs47;
        const double clhs59 =     clhs27*clhs50;
        const double clhs60 =     clhs24*clhs48;
        const double clhs61 =     clhs25*clhs49;
        const double clhs62 =     clhs57 - normalslave(0,1)*(clhs0 + clhs58 + clhs59 - clhs60 - clhs61);
        const double clhs63 =     clhs62*normalslave(0,0)*rPenaltyFactor;
        const double clhs64 =     DeltaMOperator[6](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs65 =     DeltaDOperator[6](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs66 =     DeltaDOperator[6](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs67 =     DeltaMOperator[6](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs68 =     clhs2*clhs65 + clhs4*clhs66 - clhs6*clhs64 - clhs67*clhs7;
        const double clhs69 =     clhs10*clhs65 + clhs11*clhs66 - clhs12*clhs64 - clhs13*clhs67;
        const double clhs70 =     clhs68*tangentetaslave(0,0) + clhs69*tangentetaslave(0,1);
        const double clhs71 =     clhs70*rPenaltyFactor*tangentetaslave(0,0);
        const double clhs72 =     clhs68*tangentxislave(0,0) + clhs69*tangentxislave(0,1);
        const double clhs73 =     clhs72*rPenaltyFactor*tangentxislave(0,0);
        const double clhs74 =     normalslave(0,1)*(clhs24*clhs65 + clhs25*clhs66 - clhs26*clhs64 - clhs27*clhs67);
        const double clhs75 =     clhs22*clhs64;
        const double clhs76 =     clhs23*clhs67;
        const double clhs77 =     clhs20*clhs65;
        const double clhs78 =     clhs21*clhs66;
        const double clhs79 =     clhs74 - normalslave(0,0)*(clhs75 + clhs76 - clhs77 - clhs78 + clhs8);
        const double clhs80 =     clhs79*normalslave(0,0)*rPenaltyFactor;
        const double clhs81 =     DeltaMOperator[7](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs82 =     DeltaDOperator[7](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs83 =     DeltaDOperator[7](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs84 =     DeltaMOperator[7](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs85 =     clhs2*clhs82 + clhs4*clhs83 - clhs6*clhs81 - clhs7*clhs84;
        const double clhs86 =     clhs10*clhs82 + clhs11*clhs83 - clhs12*clhs81 - clhs13*clhs84;
        const double clhs87 =     clhs85*tangentetaslave(0,0) + clhs86*tangentetaslave(0,1);
        const double clhs88 =     clhs87*rPenaltyFactor*tangentetaslave(0,0);
        const double clhs89 =     clhs85*tangentxislave(0,0) + clhs86*tangentxislave(0,1);
        const double clhs90 =     clhs89*rPenaltyFactor*tangentxislave(0,0);
        const double clhs91 =     normalslave(0,0)*(clhs20*clhs82 + clhs21*clhs83 - clhs22*clhs81 - clhs23*clhs84);
        const double clhs92 =     clhs26*clhs81;
        const double clhs93 =     clhs27*clhs84;
        const double clhs94 =     clhs24*clhs82;
        const double clhs95 =     clhs25*clhs83;
        const double clhs96 =     clhs91 - normalslave(0,1)*(clhs8 + clhs92 + clhs93 - clhs94 - clhs95);
        const double clhs97 =     clhs96*normalslave(0,0)*rPenaltyFactor;
        const double clhs98 =     DeltaMOperator[0](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs99 =     DeltaDOperator[0](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs100 =     DeltaDOperator[0](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs101 =     DeltaMOperator[0](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs102 =     clhs100*clhs4 - clhs101*clhs7 + clhs2*clhs99 - clhs6*clhs98;
        const double clhs103 =     clhs10*clhs99 + clhs100*clhs11 - clhs101*clhs13 - clhs12*clhs98;
        const double clhs104 =     clhs102*tangentetaslave(0,0) + clhs103*tangentetaslave(0,1);
        const double clhs105 =     clhs104*rPenaltyFactor*tangentetaslave(0,0);
        const double clhs106 =     clhs102*tangentxislave(0,0) + clhs103*tangentxislave(0,1);
        const double clhs107 =     clhs106*rPenaltyFactor*tangentxislave(0,0);
        const double clhs108 =     normalslave(0,0)*(clhs100*clhs21 - clhs101*clhs23 + clhs20*clhs99 - clhs22*clhs98 + clhs3) + normalslave(0,1)*(clhs100*clhs25 - clhs101*clhs27 + clhs24*clhs99 - clhs26*clhs98);
        const double clhs109 =     clhs108*normalslave(0,0)*rPenaltyFactor;
        const double clhs110 =     DeltaMOperator[1](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs111 =     DeltaDOperator[1](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs112 =     DeltaDOperator[1](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs113 =     DeltaMOperator[1](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs114 =     -clhs110*clhs6 + clhs111*clhs2 + clhs112*clhs4 - clhs113*clhs7;
        const double clhs115 =     clhs10*clhs111 + clhs11*clhs112 - clhs110*clhs12 - clhs113*clhs13;
        const double clhs116 =     clhs114*tangentetaslave(0,0) + clhs115*tangentetaslave(0,1);
        const double clhs117 =     clhs116*rPenaltyFactor*tangentetaslave(0,0);
        const double clhs118 =     clhs114*tangentxislave(0,0) + clhs115*tangentxislave(0,1);
        const double clhs119 =     clhs118*rPenaltyFactor*tangentxislave(0,0);
        const double clhs120 =     normalslave(0,0)*(-clhs110*clhs22 + clhs111*clhs20 + clhs112*clhs21 - clhs113*clhs23) + normalslave(0,1)*(-clhs110*clhs26 + clhs111*clhs24 + clhs112*clhs25 - clhs113*clhs27 + clhs3);
        const double clhs121 =     clhs120*normalslave(0,0)*rPenaltyFactor;
        const double clhs122 =     DeltaMOperator[2](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs123 =     DeltaDOperator[2](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs124 =     DeltaDOperator[2](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs125 =     DeltaMOperator[2](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs126 =     -clhs122*clhs6 + clhs123*clhs2 + clhs124*clhs4 - clhs125*clhs7;
        const double clhs127 =     clhs10*clhs123 + clhs11*clhs124 - clhs12*clhs122 - clhs125*clhs13;
        const double clhs128 =     clhs126*tangentetaslave(0,0) + clhs127*tangentetaslave(0,1);
        const double clhs129 =     clhs128*rPenaltyFactor*tangentetaslave(0,0);
        const double clhs130 =     clhs126*tangentxislave(0,0) + clhs127*tangentxislave(0,1);
        const double clhs131 =     clhs130*rPenaltyFactor*tangentxislave(0,0);
        const double clhs132 =     normalslave(0,0)*(-clhs122*clhs22 + clhs123*clhs20 + clhs124*clhs21 - clhs125*clhs23 + clhs5) + normalslave(0,1)*(-clhs122*clhs26 + clhs123*clhs24 + clhs124*clhs25 - clhs125*clhs27);
        const double clhs133 =     clhs132*normalslave(0,0)*rPenaltyFactor;
        const double clhs134 =     DeltaMOperator[3](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs135 =     DeltaDOperator[3](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs136 =     DeltaDOperator[3](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs137 =     DeltaMOperator[3](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs138 =     -clhs134*clhs6 + clhs135*clhs2 + clhs136*clhs4 - clhs137*clhs7;
        const double clhs139 =     clhs10*clhs135 + clhs11*clhs136 - clhs12*clhs134 - clhs13*clhs137;
        const double clhs140 =     clhs138*tangentetaslave(0,0) + clhs139*tangentetaslave(0,1);
        const double clhs141 =     clhs140*rPenaltyFactor*tangentetaslave(0,0);
        const double clhs142 =     clhs138*tangentxislave(0,0) + clhs139*tangentxislave(0,1);
        const double clhs143 =     clhs142*rPenaltyFactor*tangentxislave(0,0);
        const double clhs144 =     normalslave(0,0)*(-clhs134*clhs22 + clhs135*clhs20 + clhs136*clhs21 - clhs137*clhs23) + normalslave(0,1)*(-clhs134*clhs26 + clhs135*clhs24 + clhs136*clhs25 - clhs137*clhs27 + clhs5);
        const double clhs145 =     clhs144*normalslave(0,0)*rPenaltyFactor;
        const double clhs146 =     rScaleFactor*(std::pow(normalslave(0,0), 2) + std::pow(tangentetaslave(0,0), 2) + std::pow(tangentxislave(0,0), 2));
        const double clhs147 =     rScaleFactor*(normalslave(0,0)*normalslave(0,1) + tangentetaslave(0,0)*tangentetaslave(0,1) + tangentxislave(0,0)*tangentxislave(0,1));
        const double clhs148 =     clhs0*clhs147;
        const double clhs149 =     clhs15*tangentetaslave(0,1);
        const double clhs150 =     clhs17*tangentxislave(0,1);
        const double clhs151 =     clhs29*normalslave(0,1);
        const double clhs152 =     clhs36*rPenaltyFactor*tangentetaslave(0,1);
        const double clhs153 =     clhs38*rPenaltyFactor*tangentxislave(0,1);
        const double clhs154 =     clhs45*normalslave(0,1)*rPenaltyFactor;
        const double clhs155 =     clhs53*rPenaltyFactor*tangentetaslave(0,1);
        const double clhs156 =     clhs55*rPenaltyFactor*tangentxislave(0,1);
        const double clhs157 =     clhs62*normalslave(0,1)*rPenaltyFactor;
        const double clhs158 =     clhs70*rPenaltyFactor*tangentetaslave(0,1);
        const double clhs159 =     clhs72*rPenaltyFactor*tangentxislave(0,1);
        const double clhs160 =     clhs79*normalslave(0,1)*rPenaltyFactor;
        const double clhs161 =     clhs87*rPenaltyFactor*tangentetaslave(0,1);
        const double clhs162 =     clhs89*rPenaltyFactor*tangentxislave(0,1);
        const double clhs163 =     clhs96*normalslave(0,1)*rPenaltyFactor;
        const double clhs164 =     clhs104*rPenaltyFactor*tangentetaslave(0,1);
        const double clhs165 =     clhs106*rPenaltyFactor*tangentxislave(0,1);
        const double clhs166 =     clhs108*normalslave(0,1)*rPenaltyFactor;
        const double clhs167 =     clhs116*rPenaltyFactor*tangentetaslave(0,1);
        const double clhs168 =     clhs118*rPenaltyFactor*tangentxislave(0,1);
        const double clhs169 =     clhs120*normalslave(0,1)*rPenaltyFactor;
        const double clhs170 =     clhs128*rPenaltyFactor*tangentetaslave(0,1);
        const double clhs171 =     clhs130*rPenaltyFactor*tangentxislave(0,1);
        const double clhs172 =     clhs132*normalslave(0,1)*rPenaltyFactor;
        const double clhs173 =     clhs140*rPenaltyFactor*tangentetaslave(0,1);
        const double clhs174 =     clhs142*rPenaltyFactor*tangentxislave(0,1);
        const double clhs175 =     clhs144*normalslave(0,1)*rPenaltyFactor;
        const double clhs176 =     rScaleFactor*(std::pow(normalslave(0,1), 2) + std::pow(tangentetaslave(0,1), 2) + std::pow(tangentxislave(0,1), 2));
        const double clhs177 =     clhs147*clhs8;
        const double clhs178 =     clhs19 - clhs28;
        const double clhs179 =     clhs178*normalslave(0,0);
        const double clhs180 =     -clhs0;
        const double clhs181 =     clhs40 + normalslave(0,0)*(clhs180 - clhs41 - clhs42 + clhs43 + clhs44);
        const double clhs182 =     clhs181*normalslave(0,0)*rPenaltyFactor;
        const double clhs183 =     clhs57 + normalslave(0,1)*(clhs180 - clhs58 - clhs59 + clhs60 + clhs61);
        const double clhs184 =     clhs183*normalslave(0,0)*rPenaltyFactor;
        const double clhs185 =     -clhs8;
        const double clhs186 =     clhs74 + normalslave(0,0)*(clhs185 - clhs75 - clhs76 + clhs77 + clhs78);
        const double clhs187 =     clhs186*normalslave(0,0)*rPenaltyFactor;
        const double clhs188 =     clhs91 + normalslave(0,1)*(clhs185 - clhs92 - clhs93 + clhs94 + clhs95);
        const double clhs189 =     clhs188*normalslave(0,0)*rPenaltyFactor;
        const double clhs190 =     -clhs147*clhs3;
        const double clhs191 =     clhs178*normalslave(0,1);
        const double clhs192 =     clhs181*normalslave(0,1)*rPenaltyFactor;
        const double clhs193 =     clhs183*normalslave(0,1)*rPenaltyFactor;
        const double clhs194 =     clhs186*normalslave(0,1)*rPenaltyFactor;
        const double clhs195 =     clhs188*normalslave(0,1)*rPenaltyFactor;
        const double clhs196 =     -clhs147*clhs5;
    
        lhs(0,0)+=clhs0*clhs37 + clhs0*clhs39 - clhs0*clhs46 + clhs1*clhs16 + clhs1*clhs18 - clhs1*clhs30;
        lhs(0,1)+=clhs0*clhs54 + clhs0*clhs56 - clhs0*clhs63 + clhs16*clhs47 + clhs18*clhs47 - clhs30*clhs47;
        lhs(0,2)+=clhs0*clhs71 + clhs0*clhs73 - clhs0*clhs80 + clhs16*clhs64 + clhs18*clhs64 - clhs30*clhs64;
        lhs(0,3)+=clhs0*clhs88 + clhs0*clhs90 - clhs0*clhs97 + clhs16*clhs81 + clhs18*clhs81 - clhs30*clhs81;
        lhs(0,4)+=clhs0*clhs105 + clhs0*clhs107 - clhs0*clhs109 + clhs16*clhs98 + clhs18*clhs98 - clhs30*clhs98;
        lhs(0,5)+=clhs0*clhs117 + clhs0*clhs119 - clhs0*clhs121 + clhs110*clhs16 + clhs110*clhs18 - clhs110*clhs30;
        lhs(0,6)+=clhs0*clhs129 + clhs0*clhs131 - clhs0*clhs133 + clhs122*clhs16 + clhs122*clhs18 - clhs122*clhs30;
        lhs(0,7)+=clhs0*clhs141 + clhs0*clhs143 - clhs0*clhs145 + clhs134*clhs16 + clhs134*clhs18 - clhs134*clhs30;
        lhs(0,8)+=clhs0*clhs146;
        lhs(0,9)+=clhs148;
        lhs(1,0)+=clhs0*clhs152 + clhs0*clhs153 - clhs0*clhs154 + clhs1*clhs149 + clhs1*clhs150 - clhs1*clhs151;
        lhs(1,1)+=clhs0*clhs155 + clhs0*clhs156 - clhs0*clhs157 + clhs149*clhs47 + clhs150*clhs47 - clhs151*clhs47;
        lhs(1,2)+=clhs0*clhs158 + clhs0*clhs159 - clhs0*clhs160 + clhs149*clhs64 + clhs150*clhs64 - clhs151*clhs64;
        lhs(1,3)+=clhs0*clhs161 + clhs0*clhs162 - clhs0*clhs163 + clhs149*clhs81 + clhs150*clhs81 - clhs151*clhs81;
        lhs(1,4)+=clhs0*clhs164 + clhs0*clhs165 - clhs0*clhs166 + clhs149*clhs98 + clhs150*clhs98 - clhs151*clhs98;
        lhs(1,5)+=clhs0*clhs167 + clhs0*clhs168 - clhs0*clhs169 + clhs110*clhs149 + clhs110*clhs150 - clhs110*clhs151;
        lhs(1,6)+=clhs0*clhs170 + clhs0*clhs171 - clhs0*clhs172 + clhs122*clhs149 + clhs122*clhs150 - clhs122*clhs151;
        lhs(1,7)+=clhs0*clhs173 + clhs0*clhs174 - clhs0*clhs175 + clhs134*clhs149 + clhs134*clhs150 - clhs134*clhs151;
        lhs(1,8)+=clhs148;
        lhs(1,9)+=clhs0*clhs176;
        lhs(2,0)+=clhs16*clhs33 + clhs18*clhs33 - clhs30*clhs33 + clhs37*clhs8 + clhs39*clhs8 - clhs46*clhs8;
        lhs(2,1)+=clhs16*clhs50 + clhs18*clhs50 - clhs30*clhs50 + clhs54*clhs8 + clhs56*clhs8 - clhs63*clhs8;
        lhs(2,2)+=clhs16*clhs67 + clhs18*clhs67 - clhs30*clhs67 + clhs71*clhs8 + clhs73*clhs8 - clhs8*clhs80;
        lhs(2,3)+=clhs16*clhs84 + clhs18*clhs84 - clhs30*clhs84 + clhs8*clhs88 + clhs8*clhs90 - clhs8*clhs97;
        lhs(2,4)+=clhs101*clhs16 + clhs101*clhs18 - clhs101*clhs30 + clhs105*clhs8 + clhs107*clhs8 - clhs109*clhs8;
        lhs(2,5)+=clhs113*clhs16 + clhs113*clhs18 - clhs113*clhs30 + clhs117*clhs8 + clhs119*clhs8 - clhs121*clhs8;
        lhs(2,6)+=clhs125*clhs16 + clhs125*clhs18 - clhs125*clhs30 + clhs129*clhs8 + clhs131*clhs8 - clhs133*clhs8;
        lhs(2,7)+=clhs137*clhs16 + clhs137*clhs18 - clhs137*clhs30 + clhs141*clhs8 + clhs143*clhs8 - clhs145*clhs8;
        lhs(2,8)+=clhs146*clhs8;
        lhs(2,9)+=clhs177;
        lhs(3,0)+=clhs149*clhs33 + clhs150*clhs33 - clhs151*clhs33 + clhs152*clhs8 + clhs153*clhs8 - clhs154*clhs8;
        lhs(3,1)+=clhs149*clhs50 + clhs150*clhs50 - clhs151*clhs50 + clhs155*clhs8 + clhs156*clhs8 - clhs157*clhs8;
        lhs(3,2)+=clhs149*clhs67 + clhs150*clhs67 - clhs151*clhs67 + clhs158*clhs8 + clhs159*clhs8 - clhs160*clhs8;
        lhs(3,3)+=clhs149*clhs84 + clhs150*clhs84 - clhs151*clhs84 + clhs161*clhs8 + clhs162*clhs8 - clhs163*clhs8;
        lhs(3,4)+=clhs101*clhs149 + clhs101*clhs150 - clhs101*clhs151 + clhs164*clhs8 + clhs165*clhs8 - clhs166*clhs8;
        lhs(3,5)+=clhs113*clhs149 + clhs113*clhs150 - clhs113*clhs151 + clhs167*clhs8 + clhs168*clhs8 - clhs169*clhs8;
        lhs(3,6)+=clhs125*clhs149 + clhs125*clhs150 - clhs125*clhs151 + clhs170*clhs8 + clhs171*clhs8 - clhs172*clhs8;
        lhs(3,7)+=clhs137*clhs149 + clhs137*clhs150 - clhs137*clhs151 + clhs173*clhs8 + clhs174*clhs8 - clhs175*clhs8;
        lhs(3,8)+=clhs177;
        lhs(3,9)+=clhs176*clhs8;
        lhs(4,0)+=-clhs16*clhs31 - clhs179*clhs31 - clhs18*clhs31 + clhs182*clhs3 - clhs3*clhs37 - clhs3*clhs39;
        lhs(4,1)+=-clhs16*clhs48 - clhs179*clhs48 - clhs18*clhs48 + clhs184*clhs3 - clhs3*clhs54 - clhs3*clhs56;
        lhs(4,2)+=-clhs16*clhs65 - clhs179*clhs65 - clhs18*clhs65 + clhs187*clhs3 - clhs3*clhs71 - clhs3*clhs73;
        lhs(4,3)+=-clhs16*clhs82 - clhs179*clhs82 - clhs18*clhs82 + clhs189*clhs3 - clhs3*clhs88 - clhs3*clhs90;
        lhs(4,4)+=-clhs105*clhs3 - clhs107*clhs3 + clhs109*clhs3 - clhs16*clhs99 - clhs179*clhs99 - clhs18*clhs99;
        lhs(4,5)+=-clhs111*clhs16 - clhs111*clhs179 - clhs111*clhs18 - clhs117*clhs3 - clhs119*clhs3 + clhs121*clhs3;
        lhs(4,6)+=-clhs123*clhs16 - clhs123*clhs179 - clhs123*clhs18 - clhs129*clhs3 - clhs131*clhs3 + clhs133*clhs3;
        lhs(4,7)+=-clhs135*clhs16 - clhs135*clhs179 - clhs135*clhs18 - clhs141*clhs3 - clhs143*clhs3 + clhs145*clhs3;
        lhs(4,8)+=-clhs146*clhs3;
        lhs(4,9)+=clhs190;
        lhs(5,0)+=-clhs149*clhs31 - clhs150*clhs31 - clhs152*clhs3 - clhs153*clhs3 - clhs191*clhs31 + clhs192*clhs3;
        lhs(5,1)+=-clhs149*clhs48 - clhs150*clhs48 - clhs155*clhs3 - clhs156*clhs3 - clhs191*clhs48 + clhs193*clhs3;
        lhs(5,2)+=-clhs149*clhs65 - clhs150*clhs65 - clhs158*clhs3 - clhs159*clhs3 - clhs191*clhs65 + clhs194*clhs3;
        lhs(5,3)+=-clhs149*clhs82 - clhs150*clhs82 - clhs161*clhs3 - clhs162*clhs3 - clhs191*clhs82 + clhs195*clhs3;
        lhs(5,4)+=-clhs149*clhs99 - clhs150*clhs99 - clhs164*clhs3 - clhs165*clhs3 + clhs166*clhs3 - clhs191*clhs99;
        lhs(5,5)+=-clhs111*clhs149 - clhs111*clhs150 - clhs111*clhs191 - clhs167*clhs3 - clhs168*clhs3 + clhs169*clhs3;
        lhs(5,6)+=-clhs123*clhs149 - clhs123*clhs150 - clhs123*clhs191 - clhs170*clhs3 - clhs171*clhs3 + clhs172*clhs3;
        lhs(5,7)+=-clhs135*clhs149 - clhs135*clhs150 - clhs135*clhs191 - clhs173*clhs3 - clhs174*clhs3 + clhs175*clhs3;
        lhs(5,8)+=clhs190;
        lhs(5,9)+=-clhs176*clhs3;
        lhs(6,0)+=-clhs16*clhs32 - clhs179*clhs32 - clhs18*clhs32 + clhs182*clhs5 - clhs37*clhs5 - clhs39*clhs5;
        lhs(6,1)+=-clhs16*clhs49 - clhs179*clhs49 - clhs18*clhs49 + clhs184*clhs5 - clhs5*clhs54 - clhs5*clhs56;
        lhs(6,2)+=-clhs16*clhs66 - clhs179*clhs66 - clhs18*clhs66 + clhs187*clhs5 - clhs5*clhs71 - clhs5*clhs73;
        lhs(6,3)+=-clhs16*clhs83 - clhs179*clhs83 - clhs18*clhs83 + clhs189*clhs5 - clhs5*clhs88 - clhs5*clhs90;
        lhs(6,4)+=-clhs100*clhs16 - clhs100*clhs179 - clhs100*clhs18 - clhs105*clhs5 - clhs107*clhs5 + clhs109*clhs5;
        lhs(6,5)+=-clhs112*clhs16 - clhs112*clhs179 - clhs112*clhs18 - clhs117*clhs5 - clhs119*clhs5 + clhs121*clhs5;
        lhs(6,6)+=-clhs124*clhs16 - clhs124*clhs179 - clhs124*clhs18 - clhs129*clhs5 - clhs131*clhs5 + clhs133*clhs5;
        lhs(6,7)+=-clhs136*clhs16 - clhs136*clhs179 - clhs136*clhs18 - clhs141*clhs5 - clhs143*clhs5 + clhs145*clhs5;
        lhs(6,8)+=-clhs146*clhs5;
        lhs(6,9)+=clhs196;
        lhs(7,0)+=-clhs149*clhs32 - clhs150*clhs32 - clhs152*clhs5 - clhs153*clhs5 - clhs191*clhs32 + clhs192*clhs5;
        lhs(7,1)+=-clhs149*clhs49 - clhs150*clhs49 - clhs155*clhs5 - clhs156*clhs5 - clhs191*clhs49 + clhs193*clhs5;
        lhs(7,2)+=-clhs149*clhs66 - clhs150*clhs66 - clhs158*clhs5 - clhs159*clhs5 - clhs191*clhs66 + clhs194*clhs5;
        lhs(7,3)+=-clhs149*clhs83 - clhs150*clhs83 - clhs161*clhs5 - clhs162*clhs5 - clhs191*clhs83 + clhs195*clhs5;
        lhs(7,4)+=-clhs100*clhs149 - clhs100*clhs150 - clhs100*clhs191 - clhs164*clhs5 - clhs165*clhs5 + clhs166*clhs5;
        lhs(7,5)+=-clhs112*clhs149 - clhs112*clhs150 - clhs112*clhs191 - clhs167*clhs5 - clhs168*clhs5 + clhs169*clhs5;
        lhs(7,6)+=-clhs124*clhs149 - clhs124*clhs150 - clhs124*clhs191 - clhs170*clhs5 - clhs171*clhs5 + clhs172*clhs5;
        lhs(7,7)+=-clhs136*clhs149 - clhs136*clhs150 - clhs136*clhs191 - clhs173*clhs5 - clhs174*clhs5 + clhs175*clhs5;
        lhs(7,8)+=clhs196;
        lhs(7,9)+=-clhs176*clhs5;
        lhs(8,0)+=rScaleFactor*(clhs181*normalslave(0,0) - clhs36*tangentetaslave(0,0) - clhs38*tangentxislave(0,0));
        lhs(8,1)+=rScaleFactor*(clhs183*normalslave(0,0) - clhs53*tangentetaslave(0,0) - clhs55*tangentxislave(0,0));
        lhs(8,2)+=rScaleFactor*(clhs186*normalslave(0,0) - clhs70*tangentetaslave(0,0) - clhs72*tangentxislave(0,0));
        lhs(8,3)+=rScaleFactor*(clhs188*normalslave(0,0) - clhs87*tangentetaslave(0,0) - clhs89*tangentxislave(0,0));
        lhs(8,4)+=rScaleFactor*(-clhs104*tangentetaslave(0,0) - clhs106*tangentxislave(0,0) + clhs108*normalslave(0,0));
        lhs(8,5)+=rScaleFactor*(-clhs116*tangentetaslave(0,0) - clhs118*tangentxislave(0,0) + clhs120*normalslave(0,0));
        lhs(8,6)+=rScaleFactor*(-clhs128*tangentetaslave(0,0) - clhs130*tangentxislave(0,0) + clhs132*normalslave(0,0));
        lhs(8,7)+=rScaleFactor*(-clhs140*tangentetaslave(0,0) - clhs142*tangentxislave(0,0) + clhs144*normalslave(0,0));
        lhs(9,0)+=rScaleFactor*(clhs181*normalslave(0,1) - clhs36*tangentetaslave(0,1) - clhs38*tangentxislave(0,1));
        lhs(9,1)+=rScaleFactor*(clhs183*normalslave(0,1) - clhs53*tangentetaslave(0,1) - clhs55*tangentxislave(0,1));
        lhs(9,2)+=rScaleFactor*(clhs186*normalslave(0,1) - clhs70*tangentetaslave(0,1) - clhs72*tangentxislave(0,1));
        lhs(9,3)+=rScaleFactor*(clhs188*normalslave(0,1) - clhs87*tangentetaslave(0,1) - clhs89*tangentxislave(0,1));
        lhs(9,4)+=rScaleFactor*(-clhs104*tangentetaslave(0,1) - clhs106*tangentxislave(0,1) + clhs108*normalslave(0,1));
        lhs(9,5)+=rScaleFactor*(-clhs116*tangentetaslave(0,1) - clhs118*tangentxislave(0,1) + clhs120*normalslave(0,1));
        lhs(9,6)+=rScaleFactor*(-clhs128*tangentetaslave(0,1) - clhs130*tangentxislave(0,1) + clhs132*normalslave(0,1));
        lhs(9,7)+=rScaleFactor*(-clhs140*tangentetaslave(0,1) - clhs142*tangentxislave(0,1) + clhs144*normalslave(0,1));
    }
    
    // NODE 1
    if (this->GetGeometry()[1].Is(ACTIVE) == false ) // INACTIVE
    {
        const double clhs0 =     0.5*std::pow(rScaleFactor, 2.0)/rPenaltyFactor;
        const double clhs1 =     clhs0*(normalslave(1,0)*normalslave(1,1) + tangentetaslave(1,0)*tangentetaslave(1,1) + tangentxislave(1,0)*tangentxislave(1,1));
    
        lhs(10,10)+=clhs0*(std::pow(normalslave(1,0), 2) + std::pow(tangentetaslave(1,0), 2) + std::pow(tangentxislave(1,0), 2));
        lhs(10,11)+=clhs1;
        lhs(11,10)+=clhs1;
        lhs(11,11)+=clhs0*(std::pow(normalslave(1,1), 2) + std::pow(tangentetaslave(1,1), 2) + std::pow(tangentxislave(1,1), 2));
    }
    else if (this->GetGeometry()[1].Is(SLIP) == true ) // ACTIVE-SLIP
    {
        const double clhs0 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs1 =     DeltaMOperator[4](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs2 =     rScaleFactor*(lm(1,0)*normalslave(1,0) + lm(1,1)*normalslave(1,1));
        const double clhs3 =     X1(0,0) + u1(0,0);
        const double clhs4 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs5 =     X1(1,0) + u1(1,0);
        const double clhs6 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs7 =     X2(0,0) + u2(0,0);
        const double clhs8 =     X2(1,0) + u2(1,0);
        const double clhs9 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs10 =     X1(0,1) + u1(0,1);
        const double clhs11 =     X1(1,1) + u1(1,1);
        const double clhs12 =     X2(0,1) + u2(0,1);
        const double clhs13 =     X2(1,1) + u2(1,1);
        const double clhs14 =     rPenaltyFactor*(normalslave(1,0)*(-clhs0*clhs7 + clhs3*clhs4 + clhs5*clhs6 - clhs8*clhs9) + normalslave(1,1)*(-clhs0*clhs12 + clhs10*clhs4 + clhs11*clhs6 - clhs13*clhs9));
        const double clhs15 =     clhs14 - clhs2;
        const double clhs16 =     clhs15*normalslave(1,0);
        const double clhs17 =     DeltaDOperator[4](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs18 =     DeltaDOperator[4](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs19 =     DeltaMOperator[4](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs20 =     normalslave(1,1)*(-clhs1*clhs12 + clhs10*clhs17 + clhs11*clhs18 - clhs13*clhs19);
        const double clhs21 =     clhs1*clhs7;
        const double clhs22 =     clhs19*clhs8;
        const double clhs23 =     clhs17*clhs3;
        const double clhs24 =     clhs18*clhs5;
        const double clhs25 =     clhs20 - normalslave(1,0)*(clhs0 + clhs21 + clhs22 - clhs23 - clhs24);
        const double clhs26 =     clhs25*normalslave(1,0)*rPenaltyFactor;
        const double clhs27 =     lm(1,0)*tangentetaslave(1,0) + lm(1,1)*tangentetaslave(1,1);
        const double clhs28 =     lm(1,0)*tangentxislave(1,0) + lm(1,1)*tangentxislave(1,1);
        const double clhs29 =     rScaleFactor*(clhs27*tangentetaslave(1,0) + clhs28*tangentxislave(1,0));
        const double clhs30 =     X1(0,0) + u1old(0,0);
        const double clhs31 =     X1(1,0) + u1old(1,0);
        const double clhs32 =     X2(0,0) + u2old(0,0);
        const double clhs33 =     X2(1,0) + u2old(1,0);
        const double clhs34 =     -clhs0*clhs32 + clhs30*clhs4 + clhs31*clhs6 - clhs33*clhs9;
        const double clhs35 =     X1(0,1) + u1old(0,1);
        const double clhs36 =     X1(1,1) + u1old(1,1);
        const double clhs37 =     X2(0,1) + u2old(0,1);
        const double clhs38 =     X2(1,1) + u2old(1,1);
        const double clhs39 =     -clhs0*clhs37 + clhs35*clhs4 + clhs36*clhs6 - clhs38*clhs9;
        const double clhs40 =     rPenaltyFactor*(clhs34*tangentetaslave(1,0) + clhs39*tangentetaslave(1,1));
        const double clhs41 =     rPenaltyFactor*(clhs34*tangentxislave(1,0) + clhs39*tangentxislave(1,1));
        const double clhs42 =     clhs40 + clhs41;
        const double clhs43 =     clhs29 + clhs42;
        const double clhs44 =     lm(0,0)*tangentetaslave(0,0) + lm(0,1)*tangentetaslave(0,1);
        const double clhs45 =     lm(0,0)*tangentxislave(0,0) + lm(0,1)*tangentxislave(0,1);
        const double clhs46 =     rScaleFactor*(clhs44*tangentetaslave(0,0) + clhs45*tangentxislave(0,0));
        const double clhs47 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs48 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs49 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs50 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs51 =     clhs30*clhs47 + clhs31*clhs48 - clhs32*clhs49 - clhs33*clhs50;
        const double clhs52 =     clhs35*clhs47 + clhs36*clhs48 - clhs37*clhs49 - clhs38*clhs50;
        const double clhs53 =     rPenaltyFactor*(clhs51*tangentetaslave(0,0) + clhs52*tangentetaslave(0,1));
        const double clhs54 =     rPenaltyFactor*(clhs51*tangentxislave(0,0) + clhs52*tangentxislave(0,1));
        const double clhs55 =     clhs53 + clhs54;
        const double clhs56 =     clhs46 + clhs55;
        const double clhs57 =     rScaleFactor*(clhs44*tangentetaslave(0,1) + clhs45*tangentxislave(0,1));
        const double clhs58 =     clhs55 + clhs57;
        const double clhs59 =     std::pow(clhs56, 2) + std::pow(clhs58, 2);
        const double clhs60 =     std::pow(clhs59, -1.0L/2.0L);
        const double clhs61 =     rScaleFactor*(clhs27*tangentetaslave(1,1) + clhs28*tangentxislave(1,1));
        const double clhs62 =     clhs42 + clhs61;
        const double clhs63 =     std::pow(clhs43, 2) + std::pow(clhs62, 2);
        const double clhs64 =     std::pow(clhs63, -1.0L/2.0L);
        const double clhs65 =     clhs15*clhs43*clhs60*clhs64*mu[1];
        const double clhs66 =     -clhs1*clhs32 + clhs17*clhs30 + clhs18*clhs31 - clhs19*clhs33;
        const double clhs67 =     -clhs1*clhs37 + clhs17*clhs35 + clhs18*clhs36 - clhs19*clhs38;
        const double clhs68 =     clhs66*tangentetaslave(1,0) + clhs66*tangentxislave(1,0) + clhs67*tangentetaslave(1,1) + clhs67*tangentxislave(1,1);
        const double clhs69 =     clhs15*clhs60*clhs64*clhs68*mu[1]*rPenaltyFactor;
        const double clhs70 =     clhs0*clhs69;
        const double clhs71 =     clhs25*clhs43*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs72 =     DeltaDOperator[4](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs73 =     DeltaDOperator[4](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs74 =     DeltaMOperator[4](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs75 =     DeltaMOperator[4](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs76 =     clhs30*clhs72 + clhs31*clhs73 - clhs32*clhs74 - clhs33*clhs75;
        const double clhs77 =     clhs35*clhs72 + clhs36*clhs73 - clhs37*clhs74 - clhs38*clhs75;
        const double clhs78 =     clhs76*tangentetaslave(0,0) + clhs76*tangentxislave(0,0) + clhs77*tangentetaslave(0,1) + clhs77*tangentxislave(0,1);
        const double clhs79 =     clhs46 + 2*clhs53 + 2*clhs54 + clhs57;
        const double clhs80 =     std::pow(clhs59, -3.0L/2.0L);
        const double clhs81 =     clhs15*clhs43*clhs64*clhs78*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs82 =     clhs29 + 2*clhs40 + 2*clhs41 + clhs61;
        const double clhs83 =     std::pow(clhs63, -3.0L/2.0L);
        const double clhs84 =     clhs15*clhs43*clhs60*clhs68*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs85 =     DeltaMOperator[5](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs86 =     DeltaDOperator[5](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs87 =     DeltaDOperator[5](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs88 =     DeltaMOperator[5](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs89 =     normalslave(1,0)*(clhs3*clhs86 + clhs5*clhs87 - clhs7*clhs85 - clhs8*clhs88);
        const double clhs90 =     clhs12*clhs85;
        const double clhs91 =     clhs13*clhs88;
        const double clhs92 =     clhs10*clhs86;
        const double clhs93 =     clhs11*clhs87;
        const double clhs94 =     clhs89 - normalslave(1,1)*(clhs0 + clhs90 + clhs91 - clhs92 - clhs93);
        const double clhs95 =     clhs94*normalslave(1,0)*rPenaltyFactor;
        const double clhs96 =     clhs30*clhs86 + clhs31*clhs87 - clhs32*clhs85 - clhs33*clhs88;
        const double clhs97 =     clhs35*clhs86 + clhs36*clhs87 - clhs37*clhs85 - clhs38*clhs88;
        const double clhs98 =     clhs96*tangentetaslave(1,0) + clhs96*tangentxislave(1,0) + clhs97*tangentetaslave(1,1) + clhs97*tangentxislave(1,1);
        const double clhs99 =     clhs15*clhs60*clhs64*clhs98*mu[1]*rPenaltyFactor;
        const double clhs100 =     clhs0*clhs99;
        const double clhs101 =     clhs43*clhs60*clhs64*clhs94*mu[1]*rPenaltyFactor;
        const double clhs102 =     DeltaDOperator[5](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs103 =     DeltaDOperator[5](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs104 =     DeltaMOperator[5](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs105 =     DeltaMOperator[5](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs106 =     clhs102*clhs30 + clhs103*clhs31 - clhs104*clhs32 - clhs105*clhs33;
        const double clhs107 =     clhs102*clhs35 + clhs103*clhs36 - clhs104*clhs37 - clhs105*clhs38;
        const double clhs108 =     clhs106*tangentetaslave(0,0) + clhs106*tangentxislave(0,0) + clhs107*tangentetaslave(0,1) + clhs107*tangentxislave(0,1);
        const double clhs109 =     clhs108*clhs15*clhs43*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs110 =     clhs15*clhs43*clhs60*clhs82*clhs83*clhs98*mu[1]*rPenaltyFactor;
        const double clhs111 =     DeltaMOperator[6](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs112 =     DeltaDOperator[6](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs113 =     DeltaDOperator[6](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs114 =     DeltaMOperator[6](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs115 =     normalslave(1,1)*(clhs10*clhs112 + clhs11*clhs113 - clhs111*clhs12 - clhs114*clhs13);
        const double clhs116 =     clhs111*clhs7;
        const double clhs117 =     clhs114*clhs8;
        const double clhs118 =     clhs112*clhs3;
        const double clhs119 =     clhs113*clhs5;
        const double clhs120 =     clhs115 - normalslave(1,0)*(clhs116 + clhs117 - clhs118 - clhs119 + clhs9);
        const double clhs121 =     clhs120*normalslave(1,0)*rPenaltyFactor;
        const double clhs122 =     -clhs111*clhs32 + clhs112*clhs30 + clhs113*clhs31 - clhs114*clhs33;
        const double clhs123 =     -clhs111*clhs37 + clhs112*clhs35 + clhs113*clhs36 - clhs114*clhs38;
        const double clhs124 =     clhs122*tangentetaslave(1,0) + clhs122*tangentxislave(1,0) + clhs123*tangentetaslave(1,1) + clhs123*tangentxislave(1,1);
        const double clhs125 =     clhs124*clhs15*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs126 =     clhs0*clhs125;
        const double clhs127 =     clhs120*clhs43*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs128 =     DeltaDOperator[6](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs129 =     DeltaDOperator[6](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs130 =     DeltaMOperator[6](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs131 =     DeltaMOperator[6](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs132 =     clhs128*clhs30 + clhs129*clhs31 - clhs130*clhs32 - clhs131*clhs33;
        const double clhs133 =     clhs128*clhs35 + clhs129*clhs36 - clhs130*clhs37 - clhs131*clhs38;
        const double clhs134 =     clhs132*tangentetaslave(0,0) + clhs132*tangentxislave(0,0) + clhs133*tangentetaslave(0,1) + clhs133*tangentxislave(0,1);
        const double clhs135 =     clhs134*clhs15*clhs43*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs136 =     clhs124*clhs15*clhs43*clhs60*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs137 =     DeltaMOperator[7](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs138 =     DeltaDOperator[7](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs139 =     DeltaDOperator[7](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs140 =     DeltaMOperator[7](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs141 =     normalslave(1,0)*(-clhs137*clhs7 + clhs138*clhs3 + clhs139*clhs5 - clhs140*clhs8);
        const double clhs142 =     clhs12*clhs137;
        const double clhs143 =     clhs13*clhs140;
        const double clhs144 =     clhs10*clhs138;
        const double clhs145 =     clhs11*clhs139;
        const double clhs146 =     clhs141 - normalslave(1,1)*(clhs142 + clhs143 - clhs144 - clhs145 + clhs9);
        const double clhs147 =     clhs146*normalslave(1,0)*rPenaltyFactor;
        const double clhs148 =     -clhs137*clhs32 + clhs138*clhs30 + clhs139*clhs31 - clhs140*clhs33;
        const double clhs149 =     -clhs137*clhs37 + clhs138*clhs35 + clhs139*clhs36 - clhs140*clhs38;
        const double clhs150 =     clhs148*tangentetaslave(1,0) + clhs148*tangentxislave(1,0) + clhs149*tangentetaslave(1,1) + clhs149*tangentxislave(1,1);
        const double clhs151 =     clhs15*clhs150*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs152 =     clhs0*clhs151;
        const double clhs153 =     clhs146*clhs43*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs154 =     DeltaDOperator[7](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs155 =     DeltaDOperator[7](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs156 =     DeltaMOperator[7](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs157 =     DeltaMOperator[7](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs158 =     clhs154*clhs30 + clhs155*clhs31 - clhs156*clhs32 - clhs157*clhs33;
        const double clhs159 =     clhs154*clhs35 + clhs155*clhs36 - clhs156*clhs37 - clhs157*clhs38;
        const double clhs160 =     clhs158*tangentetaslave(0,0) + clhs158*tangentxislave(0,0) + clhs159*tangentetaslave(0,1) + clhs159*tangentxislave(0,1);
        const double clhs161 =     clhs15*clhs160*clhs43*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs162 =     clhs15*clhs150*clhs43*clhs60*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs163 =     DeltaMOperator[0](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs164 =     DeltaDOperator[0](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs165 =     DeltaDOperator[0](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs166 =     DeltaMOperator[0](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs167 =     normalslave(1,0)*(-clhs163*clhs7 + clhs164*clhs3 + clhs165*clhs5 - clhs166*clhs8 + clhs4) + normalslave(1,1)*(clhs10*clhs164 + clhs11*clhs165 - clhs12*clhs163 - clhs13*clhs166);
        const double clhs168 =     clhs167*normalslave(1,0)*rPenaltyFactor;
        const double clhs169 =     -clhs163*clhs32 + clhs164*clhs30 + clhs165*clhs31 - clhs166*clhs33;
        const double clhs170 =     -clhs163*clhs37 + clhs164*clhs35 + clhs165*clhs36 - clhs166*clhs38;
        const double clhs171 =     clhs169*tangentetaslave(1,0) + clhs169*tangentxislave(1,0) + clhs170*tangentetaslave(1,1) + clhs170*tangentxislave(1,1);
        const double clhs172 =     clhs15*clhs171*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs173 =     clhs0*clhs172;
        const double clhs174 =     clhs167*clhs43*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs175 =     DeltaDOperator[0](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs176 =     DeltaDOperator[0](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs177 =     DeltaMOperator[0](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs178 =     DeltaMOperator[0](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs179 =     clhs175*clhs30 + clhs176*clhs31 - clhs177*clhs32 - clhs178*clhs33;
        const double clhs180 =     clhs175*clhs35 + clhs176*clhs36 - clhs177*clhs37 - clhs178*clhs38;
        const double clhs181 =     clhs179*tangentetaslave(0,0) + clhs179*tangentxislave(0,0) + clhs180*tangentetaslave(0,1) + clhs180*tangentxislave(0,1);
        const double clhs182 =     clhs15*clhs181*clhs43*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs183 =     clhs15*clhs171*clhs43*clhs60*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs184 =     DeltaMOperator[1](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs185 =     DeltaDOperator[1](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs186 =     DeltaDOperator[1](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs187 =     DeltaMOperator[1](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs188 =     normalslave(1,0)*(-clhs184*clhs7 + clhs185*clhs3 + clhs186*clhs5 - clhs187*clhs8) + normalslave(1,1)*(clhs10*clhs185 + clhs11*clhs186 - clhs12*clhs184 - clhs13*clhs187 + clhs4);
        const double clhs189 =     clhs188*normalslave(1,0)*rPenaltyFactor;
        const double clhs190 =     -clhs184*clhs32 + clhs185*clhs30 + clhs186*clhs31 - clhs187*clhs33;
        const double clhs191 =     -clhs184*clhs37 + clhs185*clhs35 + clhs186*clhs36 - clhs187*clhs38;
        const double clhs192 =     clhs190*tangentetaslave(1,0) + clhs190*tangentxislave(1,0) + clhs191*tangentetaslave(1,1) + clhs191*tangentxislave(1,1);
        const double clhs193 =     clhs15*clhs192*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs194 =     clhs0*clhs193;
        const double clhs195 =     clhs188*clhs43*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs196 =     DeltaDOperator[1](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs197 =     DeltaDOperator[1](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs198 =     DeltaMOperator[1](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs199 =     DeltaMOperator[1](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs200 =     clhs196*clhs30 + clhs197*clhs31 - clhs198*clhs32 - clhs199*clhs33;
        const double clhs201 =     clhs196*clhs35 + clhs197*clhs36 - clhs198*clhs37 - clhs199*clhs38;
        const double clhs202 =     clhs200*tangentetaslave(0,0) + clhs200*tangentxislave(0,0) + clhs201*tangentetaslave(0,1) + clhs201*tangentxislave(0,1);
        const double clhs203 =     clhs15*clhs202*clhs43*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs204 =     clhs15*clhs192*clhs43*clhs60*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs205 =     DeltaMOperator[2](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs206 =     DeltaDOperator[2](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs207 =     DeltaDOperator[2](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs208 =     DeltaMOperator[2](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs209 =     normalslave(1,0)*(-clhs205*clhs7 + clhs206*clhs3 + clhs207*clhs5 - clhs208*clhs8 + clhs6) + normalslave(1,1)*(clhs10*clhs206 + clhs11*clhs207 - clhs12*clhs205 - clhs13*clhs208);
        const double clhs210 =     clhs209*normalslave(1,0)*rPenaltyFactor;
        const double clhs211 =     -clhs205*clhs32 + clhs206*clhs30 + clhs207*clhs31 - clhs208*clhs33;
        const double clhs212 =     -clhs205*clhs37 + clhs206*clhs35 + clhs207*clhs36 - clhs208*clhs38;
        const double clhs213 =     clhs211*tangentetaslave(1,0) + clhs211*tangentxislave(1,0) + clhs212*tangentetaslave(1,1) + clhs212*tangentxislave(1,1);
        const double clhs214 =     clhs15*clhs213*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs215 =     clhs0*clhs214;
        const double clhs216 =     clhs209*clhs43*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs217 =     DeltaDOperator[2](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs218 =     DeltaDOperator[2](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs219 =     DeltaMOperator[2](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs220 =     DeltaMOperator[2](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs221 =     clhs217*clhs30 + clhs218*clhs31 - clhs219*clhs32 - clhs220*clhs33;
        const double clhs222 =     clhs217*clhs35 + clhs218*clhs36 - clhs219*clhs37 - clhs220*clhs38;
        const double clhs223 =     clhs221*tangentetaslave(0,0) + clhs221*tangentxislave(0,0) + clhs222*tangentetaslave(0,1) + clhs222*tangentxislave(0,1);
        const double clhs224 =     clhs15*clhs223*clhs43*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs225 =     clhs15*clhs213*clhs43*clhs60*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs226 =     DeltaMOperator[3](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs227 =     DeltaDOperator[3](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs228 =     DeltaDOperator[3](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs229 =     DeltaMOperator[3](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs230 =     normalslave(1,0)*(-clhs226*clhs7 + clhs227*clhs3 + clhs228*clhs5 - clhs229*clhs8) + normalslave(1,1)*(clhs10*clhs227 + clhs11*clhs228 - clhs12*clhs226 - clhs13*clhs229 + clhs6);
        const double clhs231 =     clhs230*normalslave(1,0)*rPenaltyFactor;
        const double clhs232 =     -clhs226*clhs32 + clhs227*clhs30 + clhs228*clhs31 - clhs229*clhs33;
        const double clhs233 =     -clhs226*clhs37 + clhs227*clhs35 + clhs228*clhs36 - clhs229*clhs38;
        const double clhs234 =     clhs232*tangentetaslave(1,0) + clhs232*tangentxislave(1,0) + clhs233*tangentetaslave(1,1) + clhs233*tangentxislave(1,1);
        const double clhs235 =     clhs15*clhs234*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs236 =     clhs0*clhs235;
        const double clhs237 =     clhs230*clhs43*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs238 =     DeltaDOperator[3](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs239 =     DeltaDOperator[3](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs240 =     DeltaMOperator[3](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs241 =     DeltaMOperator[3](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs242 =     clhs238*clhs30 + clhs239*clhs31 - clhs240*clhs32 - clhs241*clhs33;
        const double clhs243 =     clhs238*clhs35 + clhs239*clhs36 - clhs240*clhs37 - clhs241*clhs38;
        const double clhs244 =     clhs242*tangentetaslave(0,0) + clhs242*tangentxislave(0,0) + clhs243*tangentetaslave(0,1) + clhs243*tangentxislave(0,1);
        const double clhs245 =     clhs15*clhs244*clhs43*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs246 =     clhs15*clhs234*clhs43*clhs60*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs247 =     tangentetaslave(0,0)*tangentetaslave(0,1) + tangentxislave(0,0)*tangentxislave(0,1);
        const double clhs248 =     clhs247*clhs58 + clhs56*(std::pow(tangentetaslave(0,0), 2) + std::pow(tangentxislave(0,0), 2));
        const double clhs249 =     clhs15*clhs248*clhs43*clhs64*clhs80*mu[1]*rScaleFactor;
        const double clhs250 =     clhs247*clhs56 + clhs58*(std::pow(tangentetaslave(0,1), 2) + std::pow(tangentxislave(0,1), 2));
        const double clhs251 =     clhs15*clhs250*clhs43*clhs64*clhs80*mu[1]*rScaleFactor;
        const double clhs252 =     std::pow(normalslave(1,0), 2);
        const double clhs253 =     -clhs14 + clhs2;
        const double clhs254 =     std::pow(tangentetaslave(1,0), 2) + std::pow(tangentxislave(1,0), 2);
        const double clhs255 =     clhs254*clhs60*clhs64*mu[1];
        const double clhs256 =     clhs60*clhs64*mu[1]*normalslave(1,0);
        const double clhs257 =     clhs256*clhs43;
        const double clhs258 =     tangentetaslave(1,0)*tangentetaslave(1,1);
        const double clhs259 =     tangentxislave(1,0)*tangentxislave(1,1);
        const double clhs260 =     clhs258 + clhs259;
        const double clhs261 =     clhs254*clhs43 + clhs260*clhs62;
        const double clhs262 =     clhs253*clhs261*clhs60*clhs83*mu[1];
        const double clhs263 =     rScaleFactor*(clhs252 - clhs253*clhs255 - clhs257 + clhs262*clhs43);
        const double clhs264 =     normalslave(1,0)*normalslave(1,1);
        const double clhs265 =     clhs260*clhs60*clhs64*mu[1];
        const double clhs266 =     -clhs253*clhs265 + clhs264;
        const double clhs267 =     clhs60*clhs64*mu[1]*normalslave(1,1);
        const double clhs268 =     clhs267*clhs43;
        const double clhs269 =     std::pow(tangentetaslave(1,1), 2) + std::pow(tangentxislave(1,1), 2);
        const double clhs270 =     clhs260*clhs43 + clhs269*clhs62;
        const double clhs271 =     clhs253*clhs270*clhs60*clhs83*mu[1];
        const double clhs272 =     rScaleFactor*(clhs266 - clhs268 + clhs271*clhs43);
        const double clhs273 =     clhs15*normalslave(1,1);
        const double clhs274 =     clhs25*normalslave(1,1)*rPenaltyFactor;
        const double clhs275 =     clhs15*clhs60*clhs62*clhs64*mu[1];
        const double clhs276 =     clhs25*clhs60*clhs62*clhs64*mu[1]*rPenaltyFactor;
        const double clhs277 =     clhs15*clhs62*clhs64*clhs78*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs278 =     clhs15*clhs60*clhs62*clhs68*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs279 =     clhs94*normalslave(1,1)*rPenaltyFactor;
        const double clhs280 =     clhs60*clhs62*clhs64*clhs94*mu[1]*rPenaltyFactor;
        const double clhs281 =     clhs108*clhs15*clhs62*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs282 =     clhs15*clhs60*clhs62*clhs82*clhs83*clhs98*mu[1]*rPenaltyFactor;
        const double clhs283 =     clhs120*normalslave(1,1)*rPenaltyFactor;
        const double clhs284 =     clhs120*clhs60*clhs62*clhs64*mu[1]*rPenaltyFactor;
        const double clhs285 =     clhs134*clhs15*clhs62*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs286 =     clhs124*clhs15*clhs60*clhs62*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs287 =     clhs146*normalslave(1,1)*rPenaltyFactor;
        const double clhs288 =     clhs146*clhs60*clhs62*clhs64*mu[1]*rPenaltyFactor;
        const double clhs289 =     clhs15*clhs160*clhs62*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs290 =     clhs15*clhs150*clhs60*clhs62*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs291 =     clhs167*normalslave(1,1)*rPenaltyFactor;
        const double clhs292 =     clhs167*clhs60*clhs62*clhs64*mu[1]*rPenaltyFactor;
        const double clhs293 =     clhs15*clhs181*clhs62*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs294 =     clhs15*clhs171*clhs60*clhs62*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs295 =     clhs188*normalslave(1,1)*rPenaltyFactor;
        const double clhs296 =     clhs188*clhs60*clhs62*clhs64*mu[1]*rPenaltyFactor;
        const double clhs297 =     clhs15*clhs202*clhs62*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs298 =     clhs15*clhs192*clhs60*clhs62*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs299 =     clhs209*normalslave(1,1)*rPenaltyFactor;
        const double clhs300 =     clhs209*clhs60*clhs62*clhs64*mu[1]*rPenaltyFactor;
        const double clhs301 =     clhs15*clhs223*clhs62*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs302 =     clhs15*clhs213*clhs60*clhs62*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs303 =     clhs230*normalslave(1,1)*rPenaltyFactor;
        const double clhs304 =     clhs230*clhs60*clhs62*clhs64*mu[1]*rPenaltyFactor;
        const double clhs305 =     clhs15*clhs244*clhs62*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs306 =     clhs15*clhs234*clhs60*clhs62*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs307 =     clhs15*clhs248*clhs62*clhs64*clhs80*mu[1]*rScaleFactor;
        const double clhs308 =     clhs15*clhs250*clhs62*clhs64*clhs80*mu[1]*rScaleFactor;
        const double clhs309 =     clhs256*clhs62;
        const double clhs310 =     rScaleFactor*(clhs262*clhs62 + clhs266 - clhs309);
        const double clhs311 =     std::pow(normalslave(1,1), 2);
        const double clhs312 =     clhs269*clhs60*clhs64*mu[1];
        const double clhs313 =     clhs267*clhs62;
        const double clhs314 =     rScaleFactor*(-clhs253*clhs312 + clhs271*clhs62 + clhs311 - clhs313);
        const double clhs315 =     clhs69*clhs9;
        const double clhs316 =     clhs9*clhs99;
        const double clhs317 =     clhs125*clhs9;
        const double clhs318 =     clhs151*clhs9;
        const double clhs319 =     clhs172*clhs9;
        const double clhs320 =     clhs193*clhs9;
        const double clhs321 =     clhs214*clhs9;
        const double clhs322 =     clhs235*clhs9;
        const double clhs323 =     clhs253*normalslave(1,0);
        const double clhs324 =     -clhs0;
        const double clhs325 =     clhs20 + normalslave(1,0)*(-clhs21 - clhs22 + clhs23 + clhs24 + clhs324);
        const double clhs326 =     clhs325*normalslave(1,0)*rPenaltyFactor;
        const double clhs327 =     clhs253*clhs43*clhs60*clhs64*mu[1];
        const double clhs328 =     clhs253*clhs60*clhs64*clhs68*mu[1]*rPenaltyFactor;
        const double clhs329 =     clhs328*clhs4;
        const double clhs330 =     clhs325*clhs43*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs331 =     clhs253*clhs43*clhs64*clhs78*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs332 =     clhs253*clhs43*clhs60*clhs68*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs333 =     clhs89 + normalslave(1,1)*(clhs324 - clhs90 - clhs91 + clhs92 + clhs93);
        const double clhs334 =     clhs333*normalslave(1,0)*rPenaltyFactor;
        const double clhs335 =     clhs253*clhs60*clhs64*clhs98*mu[1]*rPenaltyFactor;
        const double clhs336 =     clhs335*clhs4;
        const double clhs337 =     clhs333*clhs43*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs338 =     clhs108*clhs253*clhs43*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs339 =     clhs253*clhs43*clhs60*clhs82*clhs83*clhs98*mu[1]*rPenaltyFactor;
        const double clhs340 =     -clhs9;
        const double clhs341 =     clhs115 + normalslave(1,0)*(-clhs116 - clhs117 + clhs118 + clhs119 + clhs340);
        const double clhs342 =     clhs341*normalslave(1,0)*rPenaltyFactor;
        const double clhs343 =     clhs124*clhs253*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs344 =     clhs343*clhs4;
        const double clhs345 =     clhs341*clhs43*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs346 =     clhs134*clhs253*clhs43*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs347 =     clhs124*clhs253*clhs43*clhs60*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs348 =     clhs141 + normalslave(1,1)*(-clhs142 - clhs143 + clhs144 + clhs145 + clhs340);
        const double clhs349 =     clhs348*normalslave(1,0)*rPenaltyFactor;
        const double clhs350 =     clhs150*clhs253*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs351 =     clhs350*clhs4;
        const double clhs352 =     clhs348*clhs43*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs353 =     clhs160*clhs253*clhs43*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs354 =     clhs150*clhs253*clhs43*clhs60*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs355 =     clhs171*clhs253*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs356 =     clhs355*clhs4;
        const double clhs357 =     clhs181*clhs253*clhs43*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs358 =     clhs171*clhs253*clhs43*clhs60*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs359 =     clhs192*clhs253*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs360 =     clhs359*clhs4;
        const double clhs361 =     clhs202*clhs253*clhs43*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs362 =     clhs192*clhs253*clhs43*clhs60*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs363 =     clhs213*clhs253*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs364 =     clhs363*clhs4;
        const double clhs365 =     clhs223*clhs253*clhs43*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs366 =     clhs213*clhs253*clhs43*clhs60*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs367 =     clhs234*clhs253*clhs60*clhs64*mu[1]*rPenaltyFactor;
        const double clhs368 =     clhs367*clhs4;
        const double clhs369 =     clhs244*clhs253*clhs43*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs370 =     clhs234*clhs253*clhs43*clhs60*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs371 =     clhs15*clhs261*clhs60*clhs83*mu[1];
        const double clhs372 =     -clhs15*clhs255 + clhs257 + clhs371*clhs43;
        const double clhs373 =     rScaleFactor*(-clhs252 + clhs372);
        const double clhs374 =     -clhs15*clhs265;
        const double clhs375 =     -clhs264 + clhs374;
        const double clhs376 =     clhs15*clhs270*clhs60*clhs83*mu[1];
        const double clhs377 =     clhs268 + clhs376*clhs43;
        const double clhs378 =     rScaleFactor*(clhs375 + clhs377);
        const double clhs379 =     clhs253*normalslave(1,1);
        const double clhs380 =     clhs325*normalslave(1,1)*rPenaltyFactor;
        const double clhs381 =     clhs253*clhs60*clhs62*clhs64*mu[1];
        const double clhs382 =     clhs325*clhs60*clhs62*clhs64*mu[1]*rPenaltyFactor;
        const double clhs383 =     clhs253*clhs62*clhs64*clhs78*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs384 =     clhs253*clhs60*clhs62*clhs68*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs385 =     clhs333*normalslave(1,1)*rPenaltyFactor;
        const double clhs386 =     clhs333*clhs60*clhs62*clhs64*mu[1]*rPenaltyFactor;
        const double clhs387 =     clhs108*clhs253*clhs62*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs388 =     clhs253*clhs60*clhs62*clhs82*clhs83*clhs98*mu[1]*rPenaltyFactor;
        const double clhs389 =     clhs341*normalslave(1,1)*rPenaltyFactor;
        const double clhs390 =     clhs341*clhs60*clhs62*clhs64*mu[1]*rPenaltyFactor;
        const double clhs391 =     clhs134*clhs253*clhs62*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs392 =     clhs124*clhs253*clhs60*clhs62*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs393 =     clhs348*normalslave(1,1)*rPenaltyFactor;
        const double clhs394 =     clhs348*clhs60*clhs62*clhs64*mu[1]*rPenaltyFactor;
        const double clhs395 =     clhs160*clhs253*clhs62*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs396 =     clhs150*clhs253*clhs60*clhs62*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs397 =     clhs181*clhs253*clhs62*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs398 =     clhs171*clhs253*clhs60*clhs62*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs399 =     clhs202*clhs253*clhs62*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs400 =     clhs192*clhs253*clhs60*clhs62*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs401 =     clhs223*clhs253*clhs62*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs402 =     clhs213*clhs253*clhs60*clhs62*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs403 =     clhs244*clhs253*clhs62*clhs64*clhs79*clhs80*mu[1]*rPenaltyFactor;
        const double clhs404 =     clhs234*clhs253*clhs60*clhs62*clhs82*clhs83*mu[1]*rPenaltyFactor;
        const double clhs405 =     clhs309 + clhs371*clhs62;
        const double clhs406 =     rScaleFactor*(clhs375 + clhs405);
        const double clhs407 =     -clhs15*clhs312 + clhs313 + clhs376*clhs62;
        const double clhs408 =     rScaleFactor*(-clhs311 + clhs407);
        const double clhs409 =     clhs328*clhs6;
        const double clhs410 =     clhs335*clhs6;
        const double clhs411 =     clhs343*clhs6;
        const double clhs412 =     clhs350*clhs6;
        const double clhs413 =     clhs355*clhs6;
        const double clhs414 =     clhs359*clhs6;
        const double clhs415 =     clhs363*clhs6;
        const double clhs416 =     clhs367*clhs6;
        const double clhs417 =     0.5*clhs254*clhs60*clhs64*mu[1];
        const double clhs418 =     clhs253*clhs68;
        const double clhs419 =     -clhs418;
        const double clhs420 =     1.0/clhs59;
        const double clhs421 =     clhs253*clhs420*clhs43*clhs79;
        const double clhs422 =     1.0/clhs63;
        const double clhs423 =     clhs422*clhs43*clhs82;
        const double clhs424 =     clhs325*clhs43 + clhs418*clhs423 + clhs419 + clhs421*clhs78;
        const double clhs425 =     0.5*clhs260*clhs60*clhs64*mu[1];
        const double clhs426 =     clhs253*clhs420*clhs62*clhs79;
        const double clhs427 =     clhs422*clhs62*clhs82;
        const double clhs428 =     clhs325*clhs62 + clhs418*clhs427 + clhs419 + clhs426*clhs78;
        const double clhs429 =     clhs253*clhs98;
        const double clhs430 =     -clhs429;
        const double clhs431 =     clhs108*clhs421 + clhs333*clhs43 + clhs423*clhs429 + clhs430;
        const double clhs432 =     clhs108*clhs426 + clhs333*clhs62 + clhs427*clhs429 + clhs430;
        const double clhs433 =     clhs124*clhs253;
        const double clhs434 =     -clhs433;
        const double clhs435 =     clhs134*clhs421 + clhs341*clhs43 + clhs423*clhs433 + clhs434;
        const double clhs436 =     clhs134*clhs426 + clhs341*clhs62 + clhs427*clhs433 + clhs434;
        const double clhs437 =     clhs150*clhs253;
        const double clhs438 =     -clhs437;
        const double clhs439 =     clhs160*clhs421 + clhs348*clhs43 + clhs423*clhs437 + clhs438;
        const double clhs440 =     clhs160*clhs426 + clhs348*clhs62 + clhs427*clhs437 + clhs438;
        const double clhs441 =     clhs171*clhs253;
        const double clhs442 =     -clhs441;
        const double clhs443 =     clhs167*clhs43 + clhs181*clhs421 + clhs423*clhs441 + clhs442;
        const double clhs444 =     clhs167*clhs62 + clhs181*clhs426 + clhs427*clhs441 + clhs442;
        const double clhs445 =     clhs192*clhs253;
        const double clhs446 =     -clhs445;
        const double clhs447 =     clhs188*clhs43 + clhs202*clhs421 + clhs423*clhs445 + clhs446;
        const double clhs448 =     clhs188*clhs62 + clhs202*clhs426 + clhs427*clhs445 + clhs446;
        const double clhs449 =     clhs213*clhs253;
        const double clhs450 =     -clhs449;
        const double clhs451 =     clhs209*clhs43 + clhs223*clhs421 + clhs423*clhs449 + clhs450;
        const double clhs452 =     clhs209*clhs62 + clhs223*clhs426 + clhs427*clhs449 + clhs450;
        const double clhs453 =     clhs234*clhs253;
        const double clhs454 =     -clhs453;
        const double clhs455 =     clhs230*clhs43 + clhs244*clhs421 + clhs423*clhs453 + clhs454;
        const double clhs456 =     clhs230*clhs62 + clhs244*clhs426 + clhs427*clhs453 + clhs454;
        const double clhs457 =     1.0/rPenaltyFactor;
        const double clhs458 =     std::pow(rScaleFactor, 2);
        const double clhs459 =     0.5*clhs248*clhs253*clhs457*clhs458*clhs64*clhs80*mu[1];
        const double clhs460 =     0.5*clhs250*clhs253*clhs457*clhs458*clhs64*clhs80*mu[1];
        const double clhs461 =     0.5*clhs457*clhs458;
        const double clhs462 =     clhs254 + clhs372;
        const double clhs463 =     clhs258 + clhs259 + clhs374;
        const double clhs464 =     clhs405 + clhs463;
        const double clhs465 =     clhs377 + clhs463;
        const double clhs466 =     clhs269 + clhs407;
        const double clhs467 =     0.5*clhs269*clhs60*clhs64*mu[1];
    
        lhs(0,0)+=-clhs0*clhs26 + clhs0*clhs71 - clhs0*clhs81 - clhs0*clhs84 - clhs1*clhs16 + clhs1*clhs65 + clhs70;
        lhs(0,1)+=clhs0*clhs101 - clhs0*clhs109 - clhs0*clhs110 - clhs0*clhs95 + clhs100 - clhs16*clhs85 + clhs65*clhs85;
        lhs(0,2)+=-clhs0*clhs121 + clhs0*clhs127 - clhs0*clhs135 - clhs0*clhs136 - clhs111*clhs16 + clhs111*clhs65 + clhs126;
        lhs(0,3)+=-clhs0*clhs147 + clhs0*clhs153 - clhs0*clhs161 - clhs0*clhs162 - clhs137*clhs16 + clhs137*clhs65 + clhs152;
        lhs(0,4)+=-clhs0*clhs168 + clhs0*clhs174 - clhs0*clhs182 - clhs0*clhs183 - clhs16*clhs163 + clhs163*clhs65 + clhs173;
        lhs(0,5)+=-clhs0*clhs189 + clhs0*clhs195 - clhs0*clhs203 - clhs0*clhs204 - clhs16*clhs184 + clhs184*clhs65 + clhs194;
        lhs(0,6)+=-clhs0*clhs210 + clhs0*clhs216 - clhs0*clhs224 - clhs0*clhs225 - clhs16*clhs205 + clhs205*clhs65 + clhs215;
        lhs(0,7)+=-clhs0*clhs231 + clhs0*clhs237 - clhs0*clhs245 - clhs0*clhs246 - clhs16*clhs226 + clhs226*clhs65 + clhs236;
        lhs(0,8)+=-clhs0*clhs249;
        lhs(0,9)+=-clhs0*clhs251;
        lhs(0,10)+=clhs0*clhs263;
        lhs(0,11)+=clhs0*clhs272;
        lhs(1,0)+=-clhs0*clhs274 + clhs0*clhs276 - clhs0*clhs277 - clhs0*clhs278 - clhs1*clhs273 + clhs1*clhs275 + clhs70;
        lhs(1,1)+=-clhs0*clhs279 + clhs0*clhs280 - clhs0*clhs281 - clhs0*clhs282 + clhs100 - clhs273*clhs85 + clhs275*clhs85;
        lhs(1,2)+=-clhs0*clhs283 + clhs0*clhs284 - clhs0*clhs285 - clhs0*clhs286 - clhs111*clhs273 + clhs111*clhs275 + clhs126;
        lhs(1,3)+=-clhs0*clhs287 + clhs0*clhs288 - clhs0*clhs289 - clhs0*clhs290 - clhs137*clhs273 + clhs137*clhs275 + clhs152;
        lhs(1,4)+=-clhs0*clhs291 + clhs0*clhs292 - clhs0*clhs293 - clhs0*clhs294 - clhs163*clhs273 + clhs163*clhs275 + clhs173;
        lhs(1,5)+=-clhs0*clhs295 + clhs0*clhs296 - clhs0*clhs297 - clhs0*clhs298 - clhs184*clhs273 + clhs184*clhs275 + clhs194;
        lhs(1,6)+=-clhs0*clhs299 + clhs0*clhs300 - clhs0*clhs301 - clhs0*clhs302 - clhs205*clhs273 + clhs205*clhs275 + clhs215;
        lhs(1,7)+=-clhs0*clhs303 + clhs0*clhs304 - clhs0*clhs305 - clhs0*clhs306 - clhs226*clhs273 + clhs226*clhs275 + clhs236;
        lhs(1,8)+=-clhs0*clhs307;
        lhs(1,9)+=-clhs0*clhs308;
        lhs(1,10)+=clhs0*clhs310;
        lhs(1,11)+=clhs0*clhs314;
        lhs(2,0)+=-clhs16*clhs19 + clhs19*clhs65 - clhs26*clhs9 + clhs315 + clhs71*clhs9 - clhs81*clhs9 - clhs84*clhs9;
        lhs(2,1)+=clhs101*clhs9 - clhs109*clhs9 - clhs110*clhs9 - clhs16*clhs88 + clhs316 + clhs65*clhs88 - clhs9*clhs95;
        lhs(2,2)+=-clhs114*clhs16 + clhs114*clhs65 - clhs121*clhs9 + clhs127*clhs9 - clhs135*clhs9 - clhs136*clhs9 + clhs317;
        lhs(2,3)+=-clhs140*clhs16 + clhs140*clhs65 - clhs147*clhs9 + clhs153*clhs9 - clhs161*clhs9 - clhs162*clhs9 + clhs318;
        lhs(2,4)+=-clhs16*clhs166 + clhs166*clhs65 - clhs168*clhs9 + clhs174*clhs9 - clhs182*clhs9 - clhs183*clhs9 + clhs319;
        lhs(2,5)+=-clhs16*clhs187 + clhs187*clhs65 - clhs189*clhs9 + clhs195*clhs9 - clhs203*clhs9 - clhs204*clhs9 + clhs320;
        lhs(2,6)+=-clhs16*clhs208 + clhs208*clhs65 - clhs210*clhs9 + clhs216*clhs9 - clhs224*clhs9 - clhs225*clhs9 + clhs321;
        lhs(2,7)+=-clhs16*clhs229 + clhs229*clhs65 - clhs231*clhs9 + clhs237*clhs9 - clhs245*clhs9 - clhs246*clhs9 + clhs322;
        lhs(2,8)+=-clhs249*clhs9;
        lhs(2,9)+=-clhs251*clhs9;
        lhs(2,10)+=clhs263*clhs9;
        lhs(2,11)+=clhs272*clhs9;
        lhs(3,0)+=-clhs19*clhs273 + clhs19*clhs275 - clhs274*clhs9 + clhs276*clhs9 - clhs277*clhs9 - clhs278*clhs9 + clhs315;
        lhs(3,1)+=-clhs273*clhs88 + clhs275*clhs88 - clhs279*clhs9 + clhs280*clhs9 - clhs281*clhs9 - clhs282*clhs9 + clhs316;
        lhs(3,2)+=-clhs114*clhs273 + clhs114*clhs275 - clhs283*clhs9 + clhs284*clhs9 - clhs285*clhs9 - clhs286*clhs9 + clhs317;
        lhs(3,3)+=-clhs140*clhs273 + clhs140*clhs275 - clhs287*clhs9 + clhs288*clhs9 - clhs289*clhs9 - clhs290*clhs9 + clhs318;
        lhs(3,4)+=-clhs166*clhs273 + clhs166*clhs275 - clhs291*clhs9 + clhs292*clhs9 - clhs293*clhs9 - clhs294*clhs9 + clhs319;
        lhs(3,5)+=-clhs187*clhs273 + clhs187*clhs275 - clhs295*clhs9 + clhs296*clhs9 - clhs297*clhs9 - clhs298*clhs9 + clhs320;
        lhs(3,6)+=-clhs208*clhs273 + clhs208*clhs275 - clhs299*clhs9 + clhs300*clhs9 - clhs301*clhs9 - clhs302*clhs9 + clhs321;
        lhs(3,7)+=-clhs229*clhs273 + clhs229*clhs275 - clhs303*clhs9 + clhs304*clhs9 - clhs305*clhs9 - clhs306*clhs9 + clhs322;
        lhs(3,8)+=-clhs307*clhs9;
        lhs(3,9)+=-clhs308*clhs9;
        lhs(3,10)+=clhs310*clhs9;
        lhs(3,11)+=clhs314*clhs9;
        lhs(4,0)+=-clhs17*clhs323 + clhs17*clhs327 + clhs326*clhs4 + clhs329 - clhs330*clhs4 - clhs331*clhs4 - clhs332*clhs4;
        lhs(4,1)+=-clhs323*clhs86 + clhs327*clhs86 + clhs334*clhs4 + clhs336 - clhs337*clhs4 - clhs338*clhs4 - clhs339*clhs4;
        lhs(4,2)+=-clhs112*clhs323 + clhs112*clhs327 + clhs342*clhs4 + clhs344 - clhs345*clhs4 - clhs346*clhs4 - clhs347*clhs4;
        lhs(4,3)+=-clhs138*clhs323 + clhs138*clhs327 + clhs349*clhs4 + clhs351 - clhs352*clhs4 - clhs353*clhs4 - clhs354*clhs4;
        lhs(4,4)+=-clhs164*clhs323 + clhs164*clhs327 + clhs168*clhs4 - clhs174*clhs4 + clhs356 - clhs357*clhs4 - clhs358*clhs4;
        lhs(4,5)+=-clhs185*clhs323 + clhs185*clhs327 + clhs189*clhs4 - clhs195*clhs4 + clhs360 - clhs361*clhs4 - clhs362*clhs4;
        lhs(4,6)+=-clhs206*clhs323 + clhs206*clhs327 + clhs210*clhs4 - clhs216*clhs4 + clhs364 - clhs365*clhs4 - clhs366*clhs4;
        lhs(4,7)+=-clhs227*clhs323 + clhs227*clhs327 + clhs231*clhs4 - clhs237*clhs4 + clhs368 - clhs369*clhs4 - clhs370*clhs4;
        lhs(4,8)+=clhs249*clhs4;
        lhs(4,9)+=clhs251*clhs4;
        lhs(4,10)+=clhs373*clhs4;
        lhs(4,11)+=clhs378*clhs4;
        lhs(5,0)+=-clhs17*clhs379 + clhs17*clhs381 + clhs329 + clhs380*clhs4 - clhs382*clhs4 - clhs383*clhs4 - clhs384*clhs4;
        lhs(5,1)+=clhs336 - clhs379*clhs86 + clhs381*clhs86 + clhs385*clhs4 - clhs386*clhs4 - clhs387*clhs4 - clhs388*clhs4;
        lhs(5,2)+=-clhs112*clhs379 + clhs112*clhs381 + clhs344 + clhs389*clhs4 - clhs390*clhs4 - clhs391*clhs4 - clhs392*clhs4;
        lhs(5,3)+=-clhs138*clhs379 + clhs138*clhs381 + clhs351 + clhs393*clhs4 - clhs394*clhs4 - clhs395*clhs4 - clhs396*clhs4;
        lhs(5,4)+=-clhs164*clhs379 + clhs164*clhs381 + clhs291*clhs4 - clhs292*clhs4 + clhs356 - clhs397*clhs4 - clhs398*clhs4;
        lhs(5,5)+=-clhs185*clhs379 + clhs185*clhs381 + clhs295*clhs4 - clhs296*clhs4 + clhs360 - clhs399*clhs4 - clhs4*clhs400;
        lhs(5,6)+=-clhs206*clhs379 + clhs206*clhs381 + clhs299*clhs4 - clhs300*clhs4 + clhs364 - clhs4*clhs401 - clhs4*clhs402;
        lhs(5,7)+=-clhs227*clhs379 + clhs227*clhs381 + clhs303*clhs4 - clhs304*clhs4 + clhs368 - clhs4*clhs403 - clhs4*clhs404;
        lhs(5,8)+=clhs307*clhs4;
        lhs(5,9)+=clhs308*clhs4;
        lhs(5,10)+=clhs4*clhs406;
        lhs(5,11)+=clhs4*clhs408;
        lhs(6,0)+=-clhs18*clhs323 + clhs18*clhs327 + clhs326*clhs6 - clhs330*clhs6 - clhs331*clhs6 - clhs332*clhs6 + clhs409;
        lhs(6,1)+=-clhs323*clhs87 + clhs327*clhs87 + clhs334*clhs6 - clhs337*clhs6 - clhs338*clhs6 - clhs339*clhs6 + clhs410;
        lhs(6,2)+=-clhs113*clhs323 + clhs113*clhs327 + clhs342*clhs6 - clhs345*clhs6 - clhs346*clhs6 - clhs347*clhs6 + clhs411;
        lhs(6,3)+=-clhs139*clhs323 + clhs139*clhs327 + clhs349*clhs6 - clhs352*clhs6 - clhs353*clhs6 - clhs354*clhs6 + clhs412;
        lhs(6,4)+=-clhs165*clhs323 + clhs165*clhs327 + clhs168*clhs6 - clhs174*clhs6 - clhs357*clhs6 - clhs358*clhs6 + clhs413;
        lhs(6,5)+=-clhs186*clhs323 + clhs186*clhs327 + clhs189*clhs6 - clhs195*clhs6 - clhs361*clhs6 - clhs362*clhs6 + clhs414;
        lhs(6,6)+=-clhs207*clhs323 + clhs207*clhs327 + clhs210*clhs6 - clhs216*clhs6 - clhs365*clhs6 - clhs366*clhs6 + clhs415;
        lhs(6,7)+=-clhs228*clhs323 + clhs228*clhs327 + clhs231*clhs6 - clhs237*clhs6 - clhs369*clhs6 - clhs370*clhs6 + clhs416;
        lhs(6,8)+=clhs249*clhs6;
        lhs(6,9)+=clhs251*clhs6;
        lhs(6,10)+=clhs373*clhs6;
        lhs(6,11)+=clhs378*clhs6;
        lhs(7,0)+=-clhs18*clhs379 + clhs18*clhs381 + clhs380*clhs6 - clhs382*clhs6 - clhs383*clhs6 - clhs384*clhs6 + clhs409;
        lhs(7,1)+=-clhs379*clhs87 + clhs381*clhs87 + clhs385*clhs6 - clhs386*clhs6 - clhs387*clhs6 - clhs388*clhs6 + clhs410;
        lhs(7,2)+=-clhs113*clhs379 + clhs113*clhs381 + clhs389*clhs6 - clhs390*clhs6 - clhs391*clhs6 - clhs392*clhs6 + clhs411;
        lhs(7,3)+=-clhs139*clhs379 + clhs139*clhs381 + clhs393*clhs6 - clhs394*clhs6 - clhs395*clhs6 - clhs396*clhs6 + clhs412;
        lhs(7,4)+=-clhs165*clhs379 + clhs165*clhs381 + clhs291*clhs6 - clhs292*clhs6 - clhs397*clhs6 - clhs398*clhs6 + clhs413;
        lhs(7,5)+=-clhs186*clhs379 + clhs186*clhs381 + clhs295*clhs6 - clhs296*clhs6 - clhs399*clhs6 - clhs400*clhs6 + clhs414;
        lhs(7,6)+=-clhs207*clhs379 + clhs207*clhs381 + clhs299*clhs6 - clhs300*clhs6 - clhs401*clhs6 - clhs402*clhs6 + clhs415;
        lhs(7,7)+=-clhs228*clhs379 + clhs228*clhs381 + clhs303*clhs6 - clhs304*clhs6 - clhs403*clhs6 - clhs404*clhs6 + clhs416;
        lhs(7,8)+=clhs307*clhs6;
        lhs(7,9)+=clhs308*clhs6;
        lhs(7,10)+=clhs406*clhs6;
        lhs(7,11)+=clhs408*clhs6;
        lhs(10,0)+=rScaleFactor*(clhs25*normalslave(1,0) - clhs417*clhs424 - clhs425*clhs428);
        lhs(10,1)+=rScaleFactor*(-clhs417*clhs431 - clhs425*clhs432 + clhs94*normalslave(1,0));
        lhs(10,2)+=rScaleFactor*(clhs120*normalslave(1,0) - clhs417*clhs435 - clhs425*clhs436);
        lhs(10,3)+=rScaleFactor*(clhs146*normalslave(1,0) - clhs417*clhs439 - clhs425*clhs440);
        lhs(10,4)+=rScaleFactor*(clhs167*normalslave(1,0) - clhs417*clhs443 - clhs425*clhs444);
        lhs(10,5)+=rScaleFactor*(clhs188*normalslave(1,0) - clhs417*clhs447 - clhs425*clhs448);
        lhs(10,6)+=rScaleFactor*(clhs209*normalslave(1,0) - clhs417*clhs451 - clhs425*clhs452);
        lhs(10,7)+=rScaleFactor*(clhs230*normalslave(1,0) - clhs417*clhs455 - clhs425*clhs456);
        lhs(10,8)+=-clhs261*clhs459;
        lhs(10,9)+=-clhs261*clhs460;
        lhs(10,10)+=clhs461*(clhs254*clhs462 + clhs260*clhs464);
        lhs(10,11)+=clhs461*(clhs254*clhs465 + clhs260*clhs466);
        lhs(11,0)+=rScaleFactor*(clhs25*normalslave(1,1) - clhs424*clhs425 - clhs428*clhs467);
        lhs(11,1)+=rScaleFactor*(-clhs425*clhs431 - clhs432*clhs467 + clhs94*normalslave(1,1));
        lhs(11,2)+=rScaleFactor*(clhs120*normalslave(1,1) - clhs425*clhs435 - clhs436*clhs467);
        lhs(11,3)+=rScaleFactor*(clhs146*normalslave(1,1) - clhs425*clhs439 - clhs440*clhs467);
        lhs(11,4)+=rScaleFactor*(clhs167*normalslave(1,1) - clhs425*clhs443 - clhs444*clhs467);
        lhs(11,5)+=rScaleFactor*(clhs188*normalslave(1,1) - clhs425*clhs447 - clhs448*clhs467);
        lhs(11,6)+=rScaleFactor*(clhs209*normalslave(1,1) - clhs425*clhs451 - clhs452*clhs467);
        lhs(11,7)+=rScaleFactor*(clhs230*normalslave(1,1) - clhs425*clhs455 - clhs456*clhs467);
        lhs(11,8)+=-clhs270*clhs459;
        lhs(11,9)+=-clhs270*clhs460;
        lhs(11,10)+=clhs461*(clhs260*clhs462 + clhs269*clhs464);
        lhs(11,11)+=clhs461*(clhs260*clhs465 + clhs269*clhs466);
    }
    else // ACTIVE-STICK
    {
        const double clhs0 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs1 =     DeltaMOperator[4](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs2 =     X1(0,0) + u1old(0,0);
        const double clhs3 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs4 =     X1(1,0) + u1old(1,0);
        const double clhs5 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs6 =     X2(0,0) + u2old(0,0);
        const double clhs7 =     X2(1,0) + u2old(1,0);
        const double clhs8 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs9 =     -clhs0*clhs6 + clhs2*clhs3 + clhs4*clhs5 - clhs7*clhs8;
        const double clhs10 =     X1(0,1) + u1old(0,1);
        const double clhs11 =     X1(1,1) + u1old(1,1);
        const double clhs12 =     X2(0,1) + u2old(0,1);
        const double clhs13 =     X2(1,1) + u2old(1,1);
        const double clhs14 =     -clhs0*clhs12 + clhs10*clhs3 + clhs11*clhs5 - clhs13*clhs8;
        const double clhs15 =     rPenaltyFactor*(clhs14*tangentetaslave(1,1) + clhs9*tangentetaslave(1,0)) + rScaleFactor*(lm(1,0)*tangentetaslave(1,0) + lm(1,1)*tangentetaslave(1,1));
        const double clhs16 =     clhs15*tangentetaslave(1,0);
        const double clhs17 =     rPenaltyFactor*(clhs14*tangentxislave(1,1) + clhs9*tangentxislave(1,0)) + rScaleFactor*(lm(1,0)*tangentxislave(1,0) + lm(1,1)*tangentxislave(1,1));
        const double clhs18 =     clhs17*tangentxislave(1,0);
        const double clhs19 =     rScaleFactor*(lm(1,0)*normalslave(1,0) + lm(1,1)*normalslave(1,1));
        const double clhs20 =     X1(0,0) + u1(0,0);
        const double clhs21 =     X1(1,0) + u1(1,0);
        const double clhs22 =     X2(0,0) + u2(0,0);
        const double clhs23 =     X2(1,0) + u2(1,0);
        const double clhs24 =     X1(0,1) + u1(0,1);
        const double clhs25 =     X1(1,1) + u1(1,1);
        const double clhs26 =     X2(0,1) + u2(0,1);
        const double clhs27 =     X2(1,1) + u2(1,1);
        const double clhs28 =     rPenaltyFactor*(normalslave(1,0)*(-clhs0*clhs22 + clhs20*clhs3 + clhs21*clhs5 - clhs23*clhs8) + normalslave(1,1)*(-clhs0*clhs26 + clhs24*clhs3 + clhs25*clhs5 - clhs27*clhs8));
        const double clhs29 =     -clhs19 + clhs28;
        const double clhs30 =     clhs29*normalslave(1,0);
        const double clhs31 =     DeltaDOperator[4](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs32 =     DeltaDOperator[4](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs33 =     DeltaMOperator[4](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs34 =     -clhs1*clhs6 + clhs2*clhs31 + clhs32*clhs4 - clhs33*clhs7;
        const double clhs35 =     -clhs1*clhs12 + clhs10*clhs31 + clhs11*clhs32 - clhs13*clhs33;
        const double clhs36 =     clhs34*tangentetaslave(1,0) + clhs35*tangentetaslave(1,1);
        const double clhs37 =     clhs36*rPenaltyFactor*tangentetaslave(1,0);
        const double clhs38 =     clhs34*tangentxislave(1,0) + clhs35*tangentxislave(1,1);
        const double clhs39 =     clhs38*rPenaltyFactor*tangentxislave(1,0);
        const double clhs40 =     normalslave(1,1)*(-clhs1*clhs26 + clhs24*clhs31 + clhs25*clhs32 - clhs27*clhs33);
        const double clhs41 =     clhs1*clhs22;
        const double clhs42 =     clhs23*clhs33;
        const double clhs43 =     clhs20*clhs31;
        const double clhs44 =     clhs21*clhs32;
        const double clhs45 =     clhs40 - normalslave(1,0)*(clhs0 + clhs41 + clhs42 - clhs43 - clhs44);
        const double clhs46 =     clhs45*normalslave(1,0)*rPenaltyFactor;
        const double clhs47 =     DeltaMOperator[5](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs48 =     DeltaDOperator[5](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs49 =     DeltaDOperator[5](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs50 =     DeltaMOperator[5](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs51 =     clhs2*clhs48 + clhs4*clhs49 - clhs47*clhs6 - clhs50*clhs7;
        const double clhs52 =     clhs10*clhs48 + clhs11*clhs49 - clhs12*clhs47 - clhs13*clhs50;
        const double clhs53 =     clhs51*tangentetaslave(1,0) + clhs52*tangentetaslave(1,1);
        const double clhs54 =     clhs53*rPenaltyFactor*tangentetaslave(1,0);
        const double clhs55 =     clhs51*tangentxislave(1,0) + clhs52*tangentxislave(1,1);
        const double clhs56 =     clhs55*rPenaltyFactor*tangentxislave(1,0);
        const double clhs57 =     normalslave(1,0)*(clhs20*clhs48 + clhs21*clhs49 - clhs22*clhs47 - clhs23*clhs50);
        const double clhs58 =     clhs26*clhs47;
        const double clhs59 =     clhs27*clhs50;
        const double clhs60 =     clhs24*clhs48;
        const double clhs61 =     clhs25*clhs49;
        const double clhs62 =     clhs57 - normalslave(1,1)*(clhs0 + clhs58 + clhs59 - clhs60 - clhs61);
        const double clhs63 =     clhs62*normalslave(1,0)*rPenaltyFactor;
        const double clhs64 =     DeltaMOperator[6](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs65 =     DeltaDOperator[6](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs66 =     DeltaDOperator[6](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs67 =     DeltaMOperator[6](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs68 =     clhs2*clhs65 + clhs4*clhs66 - clhs6*clhs64 - clhs67*clhs7;
        const double clhs69 =     clhs10*clhs65 + clhs11*clhs66 - clhs12*clhs64 - clhs13*clhs67;
        const double clhs70 =     clhs68*tangentetaslave(1,0) + clhs69*tangentetaslave(1,1);
        const double clhs71 =     clhs70*rPenaltyFactor*tangentetaslave(1,0);
        const double clhs72 =     clhs68*tangentxislave(1,0) + clhs69*tangentxislave(1,1);
        const double clhs73 =     clhs72*rPenaltyFactor*tangentxislave(1,0);
        const double clhs74 =     normalslave(1,1)*(clhs24*clhs65 + clhs25*clhs66 - clhs26*clhs64 - clhs27*clhs67);
        const double clhs75 =     clhs22*clhs64;
        const double clhs76 =     clhs23*clhs67;
        const double clhs77 =     clhs20*clhs65;
        const double clhs78 =     clhs21*clhs66;
        const double clhs79 =     clhs74 - normalslave(1,0)*(clhs75 + clhs76 - clhs77 - clhs78 + clhs8);
        const double clhs80 =     clhs79*normalslave(1,0)*rPenaltyFactor;
        const double clhs81 =     DeltaMOperator[7](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs82 =     DeltaDOperator[7](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs83 =     DeltaDOperator[7](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs84 =     DeltaMOperator[7](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs85 =     clhs2*clhs82 + clhs4*clhs83 - clhs6*clhs81 - clhs7*clhs84;
        const double clhs86 =     clhs10*clhs82 + clhs11*clhs83 - clhs12*clhs81 - clhs13*clhs84;
        const double clhs87 =     clhs85*tangentetaslave(1,0) + clhs86*tangentetaslave(1,1);
        const double clhs88 =     clhs87*rPenaltyFactor*tangentetaslave(1,0);
        const double clhs89 =     clhs85*tangentxislave(1,0) + clhs86*tangentxislave(1,1);
        const double clhs90 =     clhs89*rPenaltyFactor*tangentxislave(1,0);
        const double clhs91 =     normalslave(1,0)*(clhs20*clhs82 + clhs21*clhs83 - clhs22*clhs81 - clhs23*clhs84);
        const double clhs92 =     clhs26*clhs81;
        const double clhs93 =     clhs27*clhs84;
        const double clhs94 =     clhs24*clhs82;
        const double clhs95 =     clhs25*clhs83;
        const double clhs96 =     clhs91 - normalslave(1,1)*(clhs8 + clhs92 + clhs93 - clhs94 - clhs95);
        const double clhs97 =     clhs96*normalslave(1,0)*rPenaltyFactor;
        const double clhs98 =     DeltaMOperator[0](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs99 =     DeltaDOperator[0](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs100 =     DeltaDOperator[0](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs101 =     DeltaMOperator[0](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs102 =     clhs100*clhs4 - clhs101*clhs7 + clhs2*clhs99 - clhs6*clhs98;
        const double clhs103 =     clhs10*clhs99 + clhs100*clhs11 - clhs101*clhs13 - clhs12*clhs98;
        const double clhs104 =     clhs102*tangentetaslave(1,0) + clhs103*tangentetaslave(1,1);
        const double clhs105 =     clhs104*rPenaltyFactor*tangentetaslave(1,0);
        const double clhs106 =     clhs102*tangentxislave(1,0) + clhs103*tangentxislave(1,1);
        const double clhs107 =     clhs106*rPenaltyFactor*tangentxislave(1,0);
        const double clhs108 =     normalslave(1,0)*(clhs100*clhs21 - clhs101*clhs23 + clhs20*clhs99 - clhs22*clhs98 + clhs3) + normalslave(1,1)*(clhs100*clhs25 - clhs101*clhs27 + clhs24*clhs99 - clhs26*clhs98);
        const double clhs109 =     clhs108*normalslave(1,0)*rPenaltyFactor;
        const double clhs110 =     DeltaMOperator[1](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs111 =     DeltaDOperator[1](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs112 =     DeltaDOperator[1](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs113 =     DeltaMOperator[1](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs114 =     -clhs110*clhs6 + clhs111*clhs2 + clhs112*clhs4 - clhs113*clhs7;
        const double clhs115 =     clhs10*clhs111 + clhs11*clhs112 - clhs110*clhs12 - clhs113*clhs13;
        const double clhs116 =     clhs114*tangentetaslave(1,0) + clhs115*tangentetaslave(1,1);
        const double clhs117 =     clhs116*rPenaltyFactor*tangentetaslave(1,0);
        const double clhs118 =     clhs114*tangentxislave(1,0) + clhs115*tangentxislave(1,1);
        const double clhs119 =     clhs118*rPenaltyFactor*tangentxislave(1,0);
        const double clhs120 =     normalslave(1,0)*(-clhs110*clhs22 + clhs111*clhs20 + clhs112*clhs21 - clhs113*clhs23) + normalslave(1,1)*(-clhs110*clhs26 + clhs111*clhs24 + clhs112*clhs25 - clhs113*clhs27 + clhs3);
        const double clhs121 =     clhs120*normalslave(1,0)*rPenaltyFactor;
        const double clhs122 =     DeltaMOperator[2](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs123 =     DeltaDOperator[2](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs124 =     DeltaDOperator[2](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs125 =     DeltaMOperator[2](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs126 =     -clhs122*clhs6 + clhs123*clhs2 + clhs124*clhs4 - clhs125*clhs7;
        const double clhs127 =     clhs10*clhs123 + clhs11*clhs124 - clhs12*clhs122 - clhs125*clhs13;
        const double clhs128 =     clhs126*tangentetaslave(1,0) + clhs127*tangentetaslave(1,1);
        const double clhs129 =     clhs128*rPenaltyFactor*tangentetaslave(1,0);
        const double clhs130 =     clhs126*tangentxislave(1,0) + clhs127*tangentxislave(1,1);
        const double clhs131 =     clhs130*rPenaltyFactor*tangentxislave(1,0);
        const double clhs132 =     normalslave(1,0)*(-clhs122*clhs22 + clhs123*clhs20 + clhs124*clhs21 - clhs125*clhs23 + clhs5) + normalslave(1,1)*(-clhs122*clhs26 + clhs123*clhs24 + clhs124*clhs25 - clhs125*clhs27);
        const double clhs133 =     clhs132*normalslave(1,0)*rPenaltyFactor;
        const double clhs134 =     DeltaMOperator[3](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs135 =     DeltaDOperator[3](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs136 =     DeltaDOperator[3](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs137 =     DeltaMOperator[3](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs138 =     -clhs134*clhs6 + clhs135*clhs2 + clhs136*clhs4 - clhs137*clhs7;
        const double clhs139 =     clhs10*clhs135 + clhs11*clhs136 - clhs12*clhs134 - clhs13*clhs137;
        const double clhs140 =     clhs138*tangentetaslave(1,0) + clhs139*tangentetaslave(1,1);
        const double clhs141 =     clhs140*rPenaltyFactor*tangentetaslave(1,0);
        const double clhs142 =     clhs138*tangentxislave(1,0) + clhs139*tangentxislave(1,1);
        const double clhs143 =     clhs142*rPenaltyFactor*tangentxislave(1,0);
        const double clhs144 =     normalslave(1,0)*(-clhs134*clhs22 + clhs135*clhs20 + clhs136*clhs21 - clhs137*clhs23) + normalslave(1,1)*(-clhs134*clhs26 + clhs135*clhs24 + clhs136*clhs25 - clhs137*clhs27 + clhs5);
        const double clhs145 =     clhs144*normalslave(1,0)*rPenaltyFactor;
        const double clhs146 =     rScaleFactor*(std::pow(normalslave(1,0), 2) + std::pow(tangentetaslave(1,0), 2) + std::pow(tangentxislave(1,0), 2));
        const double clhs147 =     rScaleFactor*(normalslave(1,0)*normalslave(1,1) + tangentetaslave(1,0)*tangentetaslave(1,1) + tangentxislave(1,0)*tangentxislave(1,1));
        const double clhs148 =     clhs0*clhs147;
        const double clhs149 =     clhs15*tangentetaslave(1,1);
        const double clhs150 =     clhs17*tangentxislave(1,1);
        const double clhs151 =     clhs29*normalslave(1,1);
        const double clhs152 =     clhs36*rPenaltyFactor*tangentetaslave(1,1);
        const double clhs153 =     clhs38*rPenaltyFactor*tangentxislave(1,1);
        const double clhs154 =     clhs45*normalslave(1,1)*rPenaltyFactor;
        const double clhs155 =     clhs53*rPenaltyFactor*tangentetaslave(1,1);
        const double clhs156 =     clhs55*rPenaltyFactor*tangentxislave(1,1);
        const double clhs157 =     clhs62*normalslave(1,1)*rPenaltyFactor;
        const double clhs158 =     clhs70*rPenaltyFactor*tangentetaslave(1,1);
        const double clhs159 =     clhs72*rPenaltyFactor*tangentxislave(1,1);
        const double clhs160 =     clhs79*normalslave(1,1)*rPenaltyFactor;
        const double clhs161 =     clhs87*rPenaltyFactor*tangentetaslave(1,1);
        const double clhs162 =     clhs89*rPenaltyFactor*tangentxislave(1,1);
        const double clhs163 =     clhs96*normalslave(1,1)*rPenaltyFactor;
        const double clhs164 =     clhs104*rPenaltyFactor*tangentetaslave(1,1);
        const double clhs165 =     clhs106*rPenaltyFactor*tangentxislave(1,1);
        const double clhs166 =     clhs108*normalslave(1,1)*rPenaltyFactor;
        const double clhs167 =     clhs116*rPenaltyFactor*tangentetaslave(1,1);
        const double clhs168 =     clhs118*rPenaltyFactor*tangentxislave(1,1);
        const double clhs169 =     clhs120*normalslave(1,1)*rPenaltyFactor;
        const double clhs170 =     clhs128*rPenaltyFactor*tangentetaslave(1,1);
        const double clhs171 =     clhs130*rPenaltyFactor*tangentxislave(1,1);
        const double clhs172 =     clhs132*normalslave(1,1)*rPenaltyFactor;
        const double clhs173 =     clhs140*rPenaltyFactor*tangentetaslave(1,1);
        const double clhs174 =     clhs142*rPenaltyFactor*tangentxislave(1,1);
        const double clhs175 =     clhs144*normalslave(1,1)*rPenaltyFactor;
        const double clhs176 =     rScaleFactor*(std::pow(normalslave(1,1), 2) + std::pow(tangentetaslave(1,1), 2) + std::pow(tangentxislave(1,1), 2));
        const double clhs177 =     clhs147*clhs8;
        const double clhs178 =     clhs19 - clhs28;
        const double clhs179 =     clhs178*normalslave(1,0);
        const double clhs180 =     -clhs0;
        const double clhs181 =     clhs40 + normalslave(1,0)*(clhs180 - clhs41 - clhs42 + clhs43 + clhs44);
        const double clhs182 =     clhs181*normalslave(1,0)*rPenaltyFactor;
        const double clhs183 =     clhs57 + normalslave(1,1)*(clhs180 - clhs58 - clhs59 + clhs60 + clhs61);
        const double clhs184 =     clhs183*normalslave(1,0)*rPenaltyFactor;
        const double clhs185 =     -clhs8;
        const double clhs186 =     clhs74 + normalslave(1,0)*(clhs185 - clhs75 - clhs76 + clhs77 + clhs78);
        const double clhs187 =     clhs186*normalslave(1,0)*rPenaltyFactor;
        const double clhs188 =     clhs91 + normalslave(1,1)*(clhs185 - clhs92 - clhs93 + clhs94 + clhs95);
        const double clhs189 =     clhs188*normalslave(1,0)*rPenaltyFactor;
        const double clhs190 =     -clhs147*clhs3;
        const double clhs191 =     clhs178*normalslave(1,1);
        const double clhs192 =     clhs181*normalslave(1,1)*rPenaltyFactor;
        const double clhs193 =     clhs183*normalslave(1,1)*rPenaltyFactor;
        const double clhs194 =     clhs186*normalslave(1,1)*rPenaltyFactor;
        const double clhs195 =     clhs188*normalslave(1,1)*rPenaltyFactor;
        const double clhs196 =     -clhs147*clhs5;
    
        lhs(0,0)+=clhs0*clhs37 + clhs0*clhs39 - clhs0*clhs46 + clhs1*clhs16 + clhs1*clhs18 - clhs1*clhs30;
        lhs(0,1)+=clhs0*clhs54 + clhs0*clhs56 - clhs0*clhs63 + clhs16*clhs47 + clhs18*clhs47 - clhs30*clhs47;
        lhs(0,2)+=clhs0*clhs71 + clhs0*clhs73 - clhs0*clhs80 + clhs16*clhs64 + clhs18*clhs64 - clhs30*clhs64;
        lhs(0,3)+=clhs0*clhs88 + clhs0*clhs90 - clhs0*clhs97 + clhs16*clhs81 + clhs18*clhs81 - clhs30*clhs81;
        lhs(0,4)+=clhs0*clhs105 + clhs0*clhs107 - clhs0*clhs109 + clhs16*clhs98 + clhs18*clhs98 - clhs30*clhs98;
        lhs(0,5)+=clhs0*clhs117 + clhs0*clhs119 - clhs0*clhs121 + clhs110*clhs16 + clhs110*clhs18 - clhs110*clhs30;
        lhs(0,6)+=clhs0*clhs129 + clhs0*clhs131 - clhs0*clhs133 + clhs122*clhs16 + clhs122*clhs18 - clhs122*clhs30;
        lhs(0,7)+=clhs0*clhs141 + clhs0*clhs143 - clhs0*clhs145 + clhs134*clhs16 + clhs134*clhs18 - clhs134*clhs30;
        lhs(0,10)+=clhs0*clhs146;
        lhs(0,11)+=clhs148;
        lhs(1,0)+=clhs0*clhs152 + clhs0*clhs153 - clhs0*clhs154 + clhs1*clhs149 + clhs1*clhs150 - clhs1*clhs151;
        lhs(1,1)+=clhs0*clhs155 + clhs0*clhs156 - clhs0*clhs157 + clhs149*clhs47 + clhs150*clhs47 - clhs151*clhs47;
        lhs(1,2)+=clhs0*clhs158 + clhs0*clhs159 - clhs0*clhs160 + clhs149*clhs64 + clhs150*clhs64 - clhs151*clhs64;
        lhs(1,3)+=clhs0*clhs161 + clhs0*clhs162 - clhs0*clhs163 + clhs149*clhs81 + clhs150*clhs81 - clhs151*clhs81;
        lhs(1,4)+=clhs0*clhs164 + clhs0*clhs165 - clhs0*clhs166 + clhs149*clhs98 + clhs150*clhs98 - clhs151*clhs98;
        lhs(1,5)+=clhs0*clhs167 + clhs0*clhs168 - clhs0*clhs169 + clhs110*clhs149 + clhs110*clhs150 - clhs110*clhs151;
        lhs(1,6)+=clhs0*clhs170 + clhs0*clhs171 - clhs0*clhs172 + clhs122*clhs149 + clhs122*clhs150 - clhs122*clhs151;
        lhs(1,7)+=clhs0*clhs173 + clhs0*clhs174 - clhs0*clhs175 + clhs134*clhs149 + clhs134*clhs150 - clhs134*clhs151;
        lhs(1,10)+=clhs148;
        lhs(1,11)+=clhs0*clhs176;
        lhs(2,0)+=clhs16*clhs33 + clhs18*clhs33 - clhs30*clhs33 + clhs37*clhs8 + clhs39*clhs8 - clhs46*clhs8;
        lhs(2,1)+=clhs16*clhs50 + clhs18*clhs50 - clhs30*clhs50 + clhs54*clhs8 + clhs56*clhs8 - clhs63*clhs8;
        lhs(2,2)+=clhs16*clhs67 + clhs18*clhs67 - clhs30*clhs67 + clhs71*clhs8 + clhs73*clhs8 - clhs8*clhs80;
        lhs(2,3)+=clhs16*clhs84 + clhs18*clhs84 - clhs30*clhs84 + clhs8*clhs88 + clhs8*clhs90 - clhs8*clhs97;
        lhs(2,4)+=clhs101*clhs16 + clhs101*clhs18 - clhs101*clhs30 + clhs105*clhs8 + clhs107*clhs8 - clhs109*clhs8;
        lhs(2,5)+=clhs113*clhs16 + clhs113*clhs18 - clhs113*clhs30 + clhs117*clhs8 + clhs119*clhs8 - clhs121*clhs8;
        lhs(2,6)+=clhs125*clhs16 + clhs125*clhs18 - clhs125*clhs30 + clhs129*clhs8 + clhs131*clhs8 - clhs133*clhs8;
        lhs(2,7)+=clhs137*clhs16 + clhs137*clhs18 - clhs137*clhs30 + clhs141*clhs8 + clhs143*clhs8 - clhs145*clhs8;
        lhs(2,10)+=clhs146*clhs8;
        lhs(2,11)+=clhs177;
        lhs(3,0)+=clhs149*clhs33 + clhs150*clhs33 - clhs151*clhs33 + clhs152*clhs8 + clhs153*clhs8 - clhs154*clhs8;
        lhs(3,1)+=clhs149*clhs50 + clhs150*clhs50 - clhs151*clhs50 + clhs155*clhs8 + clhs156*clhs8 - clhs157*clhs8;
        lhs(3,2)+=clhs149*clhs67 + clhs150*clhs67 - clhs151*clhs67 + clhs158*clhs8 + clhs159*clhs8 - clhs160*clhs8;
        lhs(3,3)+=clhs149*clhs84 + clhs150*clhs84 - clhs151*clhs84 + clhs161*clhs8 + clhs162*clhs8 - clhs163*clhs8;
        lhs(3,4)+=clhs101*clhs149 + clhs101*clhs150 - clhs101*clhs151 + clhs164*clhs8 + clhs165*clhs8 - clhs166*clhs8;
        lhs(3,5)+=clhs113*clhs149 + clhs113*clhs150 - clhs113*clhs151 + clhs167*clhs8 + clhs168*clhs8 - clhs169*clhs8;
        lhs(3,6)+=clhs125*clhs149 + clhs125*clhs150 - clhs125*clhs151 + clhs170*clhs8 + clhs171*clhs8 - clhs172*clhs8;
        lhs(3,7)+=clhs137*clhs149 + clhs137*clhs150 - clhs137*clhs151 + clhs173*clhs8 + clhs174*clhs8 - clhs175*clhs8;
        lhs(3,10)+=clhs177;
        lhs(3,11)+=clhs176*clhs8;
        lhs(4,0)+=-clhs16*clhs31 - clhs179*clhs31 - clhs18*clhs31 + clhs182*clhs3 - clhs3*clhs37 - clhs3*clhs39;
        lhs(4,1)+=-clhs16*clhs48 - clhs179*clhs48 - clhs18*clhs48 + clhs184*clhs3 - clhs3*clhs54 - clhs3*clhs56;
        lhs(4,2)+=-clhs16*clhs65 - clhs179*clhs65 - clhs18*clhs65 + clhs187*clhs3 - clhs3*clhs71 - clhs3*clhs73;
        lhs(4,3)+=-clhs16*clhs82 - clhs179*clhs82 - clhs18*clhs82 + clhs189*clhs3 - clhs3*clhs88 - clhs3*clhs90;
        lhs(4,4)+=-clhs105*clhs3 - clhs107*clhs3 + clhs109*clhs3 - clhs16*clhs99 - clhs179*clhs99 - clhs18*clhs99;
        lhs(4,5)+=-clhs111*clhs16 - clhs111*clhs179 - clhs111*clhs18 - clhs117*clhs3 - clhs119*clhs3 + clhs121*clhs3;
        lhs(4,6)+=-clhs123*clhs16 - clhs123*clhs179 - clhs123*clhs18 - clhs129*clhs3 - clhs131*clhs3 + clhs133*clhs3;
        lhs(4,7)+=-clhs135*clhs16 - clhs135*clhs179 - clhs135*clhs18 - clhs141*clhs3 - clhs143*clhs3 + clhs145*clhs3;
        lhs(4,10)+=-clhs146*clhs3;
        lhs(4,11)+=clhs190;
        lhs(5,0)+=-clhs149*clhs31 - clhs150*clhs31 - clhs152*clhs3 - clhs153*clhs3 - clhs191*clhs31 + clhs192*clhs3;
        lhs(5,1)+=-clhs149*clhs48 - clhs150*clhs48 - clhs155*clhs3 - clhs156*clhs3 - clhs191*clhs48 + clhs193*clhs3;
        lhs(5,2)+=-clhs149*clhs65 - clhs150*clhs65 - clhs158*clhs3 - clhs159*clhs3 - clhs191*clhs65 + clhs194*clhs3;
        lhs(5,3)+=-clhs149*clhs82 - clhs150*clhs82 - clhs161*clhs3 - clhs162*clhs3 - clhs191*clhs82 + clhs195*clhs3;
        lhs(5,4)+=-clhs149*clhs99 - clhs150*clhs99 - clhs164*clhs3 - clhs165*clhs3 + clhs166*clhs3 - clhs191*clhs99;
        lhs(5,5)+=-clhs111*clhs149 - clhs111*clhs150 - clhs111*clhs191 - clhs167*clhs3 - clhs168*clhs3 + clhs169*clhs3;
        lhs(5,6)+=-clhs123*clhs149 - clhs123*clhs150 - clhs123*clhs191 - clhs170*clhs3 - clhs171*clhs3 + clhs172*clhs3;
        lhs(5,7)+=-clhs135*clhs149 - clhs135*clhs150 - clhs135*clhs191 - clhs173*clhs3 - clhs174*clhs3 + clhs175*clhs3;
        lhs(5,10)+=clhs190;
        lhs(5,11)+=-clhs176*clhs3;
        lhs(6,0)+=-clhs16*clhs32 - clhs179*clhs32 - clhs18*clhs32 + clhs182*clhs5 - clhs37*clhs5 - clhs39*clhs5;
        lhs(6,1)+=-clhs16*clhs49 - clhs179*clhs49 - clhs18*clhs49 + clhs184*clhs5 - clhs5*clhs54 - clhs5*clhs56;
        lhs(6,2)+=-clhs16*clhs66 - clhs179*clhs66 - clhs18*clhs66 + clhs187*clhs5 - clhs5*clhs71 - clhs5*clhs73;
        lhs(6,3)+=-clhs16*clhs83 - clhs179*clhs83 - clhs18*clhs83 + clhs189*clhs5 - clhs5*clhs88 - clhs5*clhs90;
        lhs(6,4)+=-clhs100*clhs16 - clhs100*clhs179 - clhs100*clhs18 - clhs105*clhs5 - clhs107*clhs5 + clhs109*clhs5;
        lhs(6,5)+=-clhs112*clhs16 - clhs112*clhs179 - clhs112*clhs18 - clhs117*clhs5 - clhs119*clhs5 + clhs121*clhs5;
        lhs(6,6)+=-clhs124*clhs16 - clhs124*clhs179 - clhs124*clhs18 - clhs129*clhs5 - clhs131*clhs5 + clhs133*clhs5;
        lhs(6,7)+=-clhs136*clhs16 - clhs136*clhs179 - clhs136*clhs18 - clhs141*clhs5 - clhs143*clhs5 + clhs145*clhs5;
        lhs(6,10)+=-clhs146*clhs5;
        lhs(6,11)+=clhs196;
        lhs(7,0)+=-clhs149*clhs32 - clhs150*clhs32 - clhs152*clhs5 - clhs153*clhs5 - clhs191*clhs32 + clhs192*clhs5;
        lhs(7,1)+=-clhs149*clhs49 - clhs150*clhs49 - clhs155*clhs5 - clhs156*clhs5 - clhs191*clhs49 + clhs193*clhs5;
        lhs(7,2)+=-clhs149*clhs66 - clhs150*clhs66 - clhs158*clhs5 - clhs159*clhs5 - clhs191*clhs66 + clhs194*clhs5;
        lhs(7,3)+=-clhs149*clhs83 - clhs150*clhs83 - clhs161*clhs5 - clhs162*clhs5 - clhs191*clhs83 + clhs195*clhs5;
        lhs(7,4)+=-clhs100*clhs149 - clhs100*clhs150 - clhs100*clhs191 - clhs164*clhs5 - clhs165*clhs5 + clhs166*clhs5;
        lhs(7,5)+=-clhs112*clhs149 - clhs112*clhs150 - clhs112*clhs191 - clhs167*clhs5 - clhs168*clhs5 + clhs169*clhs5;
        lhs(7,6)+=-clhs124*clhs149 - clhs124*clhs150 - clhs124*clhs191 - clhs170*clhs5 - clhs171*clhs5 + clhs172*clhs5;
        lhs(7,7)+=-clhs136*clhs149 - clhs136*clhs150 - clhs136*clhs191 - clhs173*clhs5 - clhs174*clhs5 + clhs175*clhs5;
        lhs(7,10)+=clhs196;
        lhs(7,11)+=-clhs176*clhs5;
        lhs(10,0)+=rScaleFactor*(clhs181*normalslave(1,0) - clhs36*tangentetaslave(1,0) - clhs38*tangentxislave(1,0));
        lhs(10,1)+=rScaleFactor*(clhs183*normalslave(1,0) - clhs53*tangentetaslave(1,0) - clhs55*tangentxislave(1,0));
        lhs(10,2)+=rScaleFactor*(clhs186*normalslave(1,0) - clhs70*tangentetaslave(1,0) - clhs72*tangentxislave(1,0));
        lhs(10,3)+=rScaleFactor*(clhs188*normalslave(1,0) - clhs87*tangentetaslave(1,0) - clhs89*tangentxislave(1,0));
        lhs(10,4)+=rScaleFactor*(-clhs104*tangentetaslave(1,0) - clhs106*tangentxislave(1,0) + clhs108*normalslave(1,0));
        lhs(10,5)+=rScaleFactor*(-clhs116*tangentetaslave(1,0) - clhs118*tangentxislave(1,0) + clhs120*normalslave(1,0));
        lhs(10,6)+=rScaleFactor*(-clhs128*tangentetaslave(1,0) - clhs130*tangentxislave(1,0) + clhs132*normalslave(1,0));
        lhs(10,7)+=rScaleFactor*(-clhs140*tangentetaslave(1,0) - clhs142*tangentxislave(1,0) + clhs144*normalslave(1,0));
        lhs(11,0)+=rScaleFactor*(clhs181*normalslave(1,1) - clhs36*tangentetaslave(1,1) - clhs38*tangentxislave(1,1));
        lhs(11,1)+=rScaleFactor*(clhs183*normalslave(1,1) - clhs53*tangentetaslave(1,1) - clhs55*tangentxislave(1,1));
        lhs(11,2)+=rScaleFactor*(clhs186*normalslave(1,1) - clhs70*tangentetaslave(1,1) - clhs72*tangentxislave(1,1));
        lhs(11,3)+=rScaleFactor*(clhs188*normalslave(1,1) - clhs87*tangentetaslave(1,1) - clhs89*tangentxislave(1,1));
        lhs(11,4)+=rScaleFactor*(-clhs104*tangentetaslave(1,1) - clhs106*tangentxislave(1,1) + clhs108*normalslave(1,1));
        lhs(11,5)+=rScaleFactor*(-clhs116*tangentetaslave(1,1) - clhs118*tangentxislave(1,1) + clhs120*normalslave(1,1));
        lhs(11,6)+=rScaleFactor*(-clhs128*tangentetaslave(1,1) - clhs130*tangentxislave(1,1) + clhs132*normalslave(1,1));
        lhs(11,7)+=rScaleFactor*(-clhs140*tangentetaslave(1,1) - clhs142*tangentxislave(1,1) + clhs144*normalslave(1,1));
    }


    return lhs;
}


/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

template<>
array_1d<double, 12> AugmentedLagrangianMethodFrictionalMortarContactCondition<2,2>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const double& rPenaltyFactor,
        const double& rScaleFactor,
        const unsigned int& rActiveInactive
        )
{
    array_1d<double,12> rhs(0.0,12);

    // Master segment info
    GeometryType& CurrentMasterElement = mThisMasterElements[rMasterElementIndex]->GetGeometry();

    // Initialize values
    const bounded_matrix<double, 2, 2> u1 = ContactUtilities::GetVariableMatrix<2,2>(this->GetGeometry(), DISPLACEMENT, 0);
    const bounded_matrix<double, 2, 2> u1old = ContactUtilities::GetVariableMatrix<2,2>(this->GetGeometry(), DISPLACEMENT, 1);
    const bounded_matrix<double, 2, 2> u2 = ContactUtilities::GetVariableMatrix<2,2>(CurrentMasterElement, DISPLACEMENT, 0);
    const bounded_matrix<double, 2, 2> u2old = ContactUtilities::GetVariableMatrix<2,2>(CurrentMasterElement, DISPLACEMENT, 1);
    const bounded_matrix<double, 2, 2> X1 = ContactUtilities::GetCoordinates<2,2>(this->GetGeometry(), false);
    const bounded_matrix<double, 2, 2> X2 = ContactUtilities::GetCoordinates<2,2>(CurrentMasterElement, false);
    
    const bounded_matrix<double, 2, 2> lm = ContactUtilities::GetVariableMatrix<2,2>(this->GetGeometry(), VECTOR_LAGRANGE_MULTIPLIER, 0); 
    
    const bounded_matrix<double, 2, 2> normalslave = ContactUtilities::GetVariableMatrix<2,2>(this->GetGeometry(),  NORMAL);
    const bounded_matrix<double, 2, 2> tangentxislave = ContactUtilities::GetVariableMatrix<2,2>(this->GetGeometry(),  TANGENT_XI);
    const bounded_matrix<double, 2, 2> tangentetaslave = ContactUtilities::GetVariableMatrix<2,2>(this->GetGeometry(),  TANGENT_ETA);
    
    // Mortar operators
    const bounded_matrix<double, 2, 2> MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 2, 2> DOperator = rMortarConditionMatrices.DOperator;
    // We get the friction coefficient

    const array_1d<double, 2> mu = GetFrictionCoefficient();

    // NODE 0
    if (this->GetGeometry()[0].Is(ACTIVE) == false ) // INACTIVE
    {
        const double crhs0 =     0.5*std::pow(rScaleFactor, 2.0)/rPenaltyFactor;
        const double crhs1 =     lm(0,0)*normalslave(0,0) + lm(0,1)*normalslave(0,1);
        const double crhs2 =     lm(0,0)*tangentetaslave(0,0) + lm(0,1)*tangentetaslave(0,1);
        const double crhs3 =     lm(0,0)*tangentxislave(0,0) + lm(0,1)*tangentxislave(0,1);
    
        rhs[8]+=-crhs0*(crhs1*normalslave(0,0) + crhs2*tangentetaslave(0,0) + crhs3*tangentxislave(0,0));
        rhs[9]+=-crhs0*(crhs1*normalslave(0,1) + crhs2*tangentetaslave(0,1) + crhs3*tangentxislave(0,1));
    }
    else if (this->GetGeometry()[0].Is(SLIP) == true ) // ACTIVE-SLIP
    {
        const double crhs0 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs1 =     rScaleFactor*(lm(0,0)*normalslave(0,0) + lm(0,1)*normalslave(0,1));
        const double crhs2 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs3 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs4 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs5 =     normalslave(0,0)*(-crhs0*(X2(0,0) + u2(0,0)) + crhs2*(X1(0,0) + u1(0,0)) + crhs3*(X1(1,0) + u1(1,0)) - crhs4*(X2(1,0) + u2(1,0))) + normalslave(0,1)*(-crhs0*(X2(0,1) + u2(0,1)) + crhs2*(X1(0,1) + u1(0,1)) + crhs3*(X1(1,1) + u1(1,1)) - crhs4*(X2(1,1) + u2(1,1)));
        const double crhs6 =     crhs5*rPenaltyFactor;
        const double crhs7 =     -crhs1 + crhs6;
        const double crhs8 =     lm(0,0)*tangentetaslave(0,0) + lm(0,1)*tangentetaslave(0,1);
        const double crhs9 =     lm(0,0)*tangentxislave(0,0) + lm(0,1)*tangentxislave(0,1);
        const double crhs10 =     rScaleFactor*(crhs8*tangentetaslave(0,0) + crhs9*tangentxislave(0,0));
        const double crhs11 =     X1(0,0) + u1old(0,0);
        const double crhs12 =     X1(1,0) + u1old(1,0);
        const double crhs13 =     X2(0,0) + u2old(0,0);
        const double crhs14 =     X2(1,0) + u2old(1,0);
        const double crhs15 =     -crhs0*crhs13 + crhs11*crhs2 + crhs12*crhs3 - crhs14*crhs4;
        const double crhs16 =     X1(0,1) + u1old(0,1);
        const double crhs17 =     X1(1,1) + u1old(1,1);
        const double crhs18 =     X2(0,1) + u2old(0,1);
        const double crhs19 =     X2(1,1) + u2old(1,1);
        const double crhs20 =     -crhs0*crhs18 + crhs16*crhs2 + crhs17*crhs3 - crhs19*crhs4;
        const double crhs21 =     rPenaltyFactor*(crhs15*tangentetaslave(0,0) + crhs20*tangentetaslave(0,1)) + rPenaltyFactor*(crhs15*tangentxislave(0,0) + crhs20*tangentxislave(0,1));
        const double crhs22 =     crhs10 + crhs21;
        const double crhs23 =     rScaleFactor*(crhs8*tangentetaslave(0,1) + crhs9*tangentxislave(0,1));
        const double crhs24 =     crhs21 + crhs23;
        const double crhs25 =     lm(1,0)*tangentetaslave(1,0) + lm(1,1)*tangentetaslave(1,1);
        const double crhs26 =     lm(1,0)*tangentxislave(1,0) + lm(1,1)*tangentxislave(1,1);
        const double crhs27 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs28 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs29 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs30 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs31 =     crhs11*crhs27 + crhs12*crhs28 - crhs13*crhs29 - crhs14*crhs30;
        const double crhs32 =     crhs16*crhs27 + crhs17*crhs28 - crhs18*crhs29 - crhs19*crhs30;
        const double crhs33 =     rPenaltyFactor*(crhs31*tangentetaslave(1,0) + crhs32*tangentetaslave(1,1)) + rPenaltyFactor*(crhs31*tangentxislave(1,0) + crhs32*tangentxislave(1,1));
        const double crhs34 =     mu[0]/(std::sqrt(std::pow(crhs22, 2) + std::pow(crhs24, 2))*std::sqrt(std::pow(crhs33 + rScaleFactor*(crhs25*tangentetaslave(1,0) + crhs26*tangentxislave(1,0)), 2) + std::pow(crhs33 + rScaleFactor*(crhs25*tangentetaslave(1,1) + crhs26*tangentxislave(1,1)), 2)));
        const double crhs35 =     crhs22*crhs34;
        const double crhs36 =     -crhs35 + normalslave(0,0);
        const double crhs37 =     crhs36*crhs7;
        const double crhs38 =     crhs24*crhs34;
        const double crhs39 =     -crhs38 + normalslave(0,1);
        const double crhs40 =     crhs39*crhs7;
        const double crhs41 =     crhs1 - crhs6;
        const double crhs42 =     crhs36*crhs41;
        const double crhs43 =     crhs39*crhs41;
        const double crhs44 =     1.0/rPenaltyFactor;
        const double crhs45 =     0.5*crhs44*(crhs10 + crhs35*crhs41);
        const double crhs46 =     tangentetaslave(0,0)*tangentetaslave(0,1) + tangentxislave(0,0)*tangentxislave(0,1);
        const double crhs47 =     0.5*crhs44*(crhs23 + crhs38*crhs41);
    
        rhs[0]+=crhs0*crhs37;
        rhs[1]+=crhs0*crhs40;
        rhs[2]+=crhs37*crhs4;
        rhs[3]+=crhs4*crhs40;
        rhs[4]+=crhs2*crhs42;
        rhs[5]+=crhs2*crhs43;
        rhs[6]+=crhs3*crhs42;
        rhs[7]+=crhs3*crhs43;
        rhs[8]+=-rScaleFactor*(crhs45*(std::pow(tangentetaslave(0,0), 2) + std::pow(tangentxislave(0,0), 2)) + crhs46*crhs47 + crhs5*normalslave(0,0));
        rhs[9]+=-rScaleFactor*(crhs45*crhs46 + crhs47*(std::pow(tangentetaslave(0,1), 2) + std::pow(tangentxislave(0,1), 2)) + crhs5*normalslave(0,1));
    }
    else // ACTIVE-STICK
    {
        const double crhs0 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs1 =     rScaleFactor*(lm(0,0)*normalslave(0,0) + lm(0,1)*normalslave(0,1));
        const double crhs2 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs3 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs4 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs5 =     normalslave(0,0)*(-crhs0*(X2(0,0) + u2(0,0)) + crhs2*(X1(0,0) + u1(0,0)) + crhs3*(X1(1,0) + u1(1,0)) - crhs4*(X2(1,0) + u2(1,0))) + normalslave(0,1)*(-crhs0*(X2(0,1) + u2(0,1)) + crhs2*(X1(0,1) + u1(0,1)) + crhs3*(X1(1,1) + u1(1,1)) - crhs4*(X2(1,1) + u2(1,1)));
        const double crhs6 =     crhs5*rPenaltyFactor;
        const double crhs7 =     crhs1 - crhs6;
        const double crhs8 =     -crhs0*(X2(0,0) + u2old(0,0)) + crhs2*(X1(0,0) + u1old(0,0)) + crhs3*(X1(1,0) + u1old(1,0)) - crhs4*(X2(1,0) + u2old(1,0));
        const double crhs9 =     -crhs0*(X2(0,1) + u2old(0,1)) + crhs2*(X1(0,1) + u1old(0,1)) + crhs3*(X1(1,1) + u1old(1,1)) - crhs4*(X2(1,1) + u2old(1,1));
        const double crhs10 =     crhs8*tangentetaslave(0,0) + crhs9*tangentetaslave(0,1);
        const double crhs11 =     crhs10*rPenaltyFactor + rScaleFactor*(lm(0,0)*tangentetaslave(0,0) + lm(0,1)*tangentetaslave(0,1));
        const double crhs12 =     crhs8*tangentxislave(0,0) + crhs9*tangentxislave(0,1);
        const double crhs13 =     crhs12*rPenaltyFactor + rScaleFactor*(lm(0,0)*tangentxislave(0,0) + lm(0,1)*tangentxislave(0,1));
        const double crhs14 =     crhs11*tangentetaslave(0,0) + crhs13*tangentxislave(0,0);
        const double crhs15 =     crhs14 + crhs7*normalslave(0,0);
        const double crhs16 =     crhs11*tangentetaslave(0,1) + crhs13*tangentxislave(0,1);
        const double crhs17 =     crhs16 + crhs7*normalslave(0,1);
        const double crhs18 =     -crhs1 + crhs6;
        const double crhs19 =     crhs14 - crhs18*normalslave(0,0);
        const double crhs20 =     crhs16 - crhs18*normalslave(0,1);
    
        rhs[0]+=-crhs0*crhs15;
        rhs[1]+=-crhs0*crhs17;
        rhs[2]+=-crhs15*crhs4;
        rhs[3]+=-crhs17*crhs4;
        rhs[4]+=crhs19*crhs2;
        rhs[5]+=crhs2*crhs20;
        rhs[6]+=crhs19*crhs3;
        rhs[7]+=crhs20*crhs3;
        rhs[8]+=rScaleFactor*(crhs10*tangentetaslave(0,0) + crhs12*tangentxislave(0,0) - crhs5*normalslave(0,0));
        rhs[9]+=rScaleFactor*(crhs10*tangentetaslave(0,1) + crhs12*tangentxislave(0,1) - crhs5*normalslave(0,1));
    }
    
    // NODE 1
    if (this->GetGeometry()[1].Is(ACTIVE) == false ) // INACTIVE
    {
        const double crhs0 =     0.5*std::pow(rScaleFactor, 2.0)/rPenaltyFactor;
        const double crhs1 =     lm(1,0)*normalslave(1,0) + lm(1,1)*normalslave(1,1);
        const double crhs2 =     lm(1,0)*tangentetaslave(1,0) + lm(1,1)*tangentetaslave(1,1);
        const double crhs3 =     lm(1,0)*tangentxislave(1,0) + lm(1,1)*tangentxislave(1,1);
    
        rhs[10]+=-crhs0*(crhs1*normalslave(1,0) + crhs2*tangentetaslave(1,0) + crhs3*tangentxislave(1,0));
        rhs[11]+=-crhs0*(crhs1*normalslave(1,1) + crhs2*tangentetaslave(1,1) + crhs3*tangentxislave(1,1));
    }
    else if (this->GetGeometry()[1].Is(SLIP) == true ) // ACTIVE-SLIP
    {
        const double crhs0 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs1 =     rScaleFactor*(lm(1,0)*normalslave(1,0) + lm(1,1)*normalslave(1,1));
        const double crhs2 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs3 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs4 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs5 =     normalslave(1,0)*(-crhs0*(X2(0,0) + u2(0,0)) + crhs2*(X1(0,0) + u1(0,0)) + crhs3*(X1(1,0) + u1(1,0)) - crhs4*(X2(1,0) + u2(1,0))) + normalslave(1,1)*(-crhs0*(X2(0,1) + u2(0,1)) + crhs2*(X1(0,1) + u1(0,1)) + crhs3*(X1(1,1) + u1(1,1)) - crhs4*(X2(1,1) + u2(1,1)));
        const double crhs6 =     crhs5*rPenaltyFactor;
        const double crhs7 =     -crhs1 + crhs6;
        const double crhs8 =     lm(1,0)*tangentetaslave(1,0) + lm(1,1)*tangentetaslave(1,1);
        const double crhs9 =     lm(1,0)*tangentxislave(1,0) + lm(1,1)*tangentxislave(1,1);
        const double crhs10 =     rScaleFactor*(crhs8*tangentetaslave(1,0) + crhs9*tangentxislave(1,0));
        const double crhs11 =     X1(0,0) + u1old(0,0);
        const double crhs12 =     X1(1,0) + u1old(1,0);
        const double crhs13 =     X2(0,0) + u2old(0,0);
        const double crhs14 =     X2(1,0) + u2old(1,0);
        const double crhs15 =     -crhs0*crhs13 + crhs11*crhs2 + crhs12*crhs3 - crhs14*crhs4;
        const double crhs16 =     X1(0,1) + u1old(0,1);
        const double crhs17 =     X1(1,1) + u1old(1,1);
        const double crhs18 =     X2(0,1) + u2old(0,1);
        const double crhs19 =     X2(1,1) + u2old(1,1);
        const double crhs20 =     -crhs0*crhs18 + crhs16*crhs2 + crhs17*crhs3 - crhs19*crhs4;
        const double crhs21 =     rPenaltyFactor*(crhs15*tangentetaslave(1,0) + crhs20*tangentetaslave(1,1)) + rPenaltyFactor*(crhs15*tangentxislave(1,0) + crhs20*tangentxislave(1,1));
        const double crhs22 =     crhs10 + crhs21;
        const double crhs23 =     lm(0,0)*tangentetaslave(0,0) + lm(0,1)*tangentetaslave(0,1);
        const double crhs24 =     lm(0,0)*tangentxislave(0,0) + lm(0,1)*tangentxislave(0,1);
        const double crhs25 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs26 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs27 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs28 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs29 =     crhs11*crhs25 + crhs12*crhs26 - crhs13*crhs27 - crhs14*crhs28;
        const double crhs30 =     crhs16*crhs25 + crhs17*crhs26 - crhs18*crhs27 - crhs19*crhs28;
        const double crhs31 =     rPenaltyFactor*(crhs29*tangentetaslave(0,0) + crhs30*tangentetaslave(0,1)) + rPenaltyFactor*(crhs29*tangentxislave(0,0) + crhs30*tangentxislave(0,1));
        const double crhs32 =     rScaleFactor*(crhs8*tangentetaslave(1,1) + crhs9*tangentxislave(1,1));
        const double crhs33 =     crhs21 + crhs32;
        const double crhs34 =     mu[1]/(std::sqrt(std::pow(crhs22, 2) + std::pow(crhs33, 2))*std::sqrt(std::pow(crhs31 + rScaleFactor*(crhs23*tangentetaslave(0,0) + crhs24*tangentxislave(0,0)), 2) + std::pow(crhs31 + rScaleFactor*(crhs23*tangentetaslave(0,1) + crhs24*tangentxislave(0,1)), 2)));
        const double crhs35 =     crhs22*crhs34;
        const double crhs36 =     -crhs35 + normalslave(1,0);
        const double crhs37 =     crhs36*crhs7;
        const double crhs38 =     crhs33*crhs34;
        const double crhs39 =     -crhs38 + normalslave(1,1);
        const double crhs40 =     crhs39*crhs7;
        const double crhs41 =     crhs1 - crhs6;
        const double crhs42 =     crhs36*crhs41;
        const double crhs43 =     crhs39*crhs41;
        const double crhs44 =     1.0/rPenaltyFactor;
        const double crhs45 =     0.5*crhs44*(crhs10 + crhs35*crhs41);
        const double crhs46 =     tangentetaslave(1,0)*tangentetaslave(1,1) + tangentxislave(1,0)*tangentxislave(1,1);
        const double crhs47 =     0.5*crhs44*(crhs32 + crhs38*crhs41);
    
        rhs[0]+=crhs0*crhs37;
        rhs[1]+=crhs0*crhs40;
        rhs[2]+=crhs37*crhs4;
        rhs[3]+=crhs4*crhs40;
        rhs[4]+=crhs2*crhs42;
        rhs[5]+=crhs2*crhs43;
        rhs[6]+=crhs3*crhs42;
        rhs[7]+=crhs3*crhs43;
        rhs[10]+=-rScaleFactor*(crhs45*(std::pow(tangentetaslave(1,0), 2) + std::pow(tangentxislave(1,0), 2)) + crhs46*crhs47 + crhs5*normalslave(1,0));
        rhs[11]+=-rScaleFactor*(crhs45*crhs46 + crhs47*(std::pow(tangentetaslave(1,1), 2) + std::pow(tangentxislave(1,1), 2)) + crhs5*normalslave(1,1));
    }
    else // ACTIVE-STICK
    {
        const double crhs0 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs1 =     rScaleFactor*(lm(1,0)*normalslave(1,0) + lm(1,1)*normalslave(1,1));
        const double crhs2 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs3 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs4 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs5 =     normalslave(1,0)*(-crhs0*(X2(0,0) + u2(0,0)) + crhs2*(X1(0,0) + u1(0,0)) + crhs3*(X1(1,0) + u1(1,0)) - crhs4*(X2(1,0) + u2(1,0))) + normalslave(1,1)*(-crhs0*(X2(0,1) + u2(0,1)) + crhs2*(X1(0,1) + u1(0,1)) + crhs3*(X1(1,1) + u1(1,1)) - crhs4*(X2(1,1) + u2(1,1)));
        const double crhs6 =     crhs5*rPenaltyFactor;
        const double crhs7 =     crhs1 - crhs6;
        const double crhs8 =     -crhs0*(X2(0,0) + u2old(0,0)) + crhs2*(X1(0,0) + u1old(0,0)) + crhs3*(X1(1,0) + u1old(1,0)) - crhs4*(X2(1,0) + u2old(1,0));
        const double crhs9 =     -crhs0*(X2(0,1) + u2old(0,1)) + crhs2*(X1(0,1) + u1old(0,1)) + crhs3*(X1(1,1) + u1old(1,1)) - crhs4*(X2(1,1) + u2old(1,1));
        const double crhs10 =     crhs8*tangentetaslave(1,0) + crhs9*tangentetaslave(1,1);
        const double crhs11 =     crhs10*rPenaltyFactor + rScaleFactor*(lm(1,0)*tangentetaslave(1,0) + lm(1,1)*tangentetaslave(1,1));
        const double crhs12 =     crhs8*tangentxislave(1,0) + crhs9*tangentxislave(1,1);
        const double crhs13 =     crhs12*rPenaltyFactor + rScaleFactor*(lm(1,0)*tangentxislave(1,0) + lm(1,1)*tangentxislave(1,1));
        const double crhs14 =     crhs11*tangentetaslave(1,0) + crhs13*tangentxislave(1,0);
        const double crhs15 =     crhs14 + crhs7*normalslave(1,0);
        const double crhs16 =     crhs11*tangentetaslave(1,1) + crhs13*tangentxislave(1,1);
        const double crhs17 =     crhs16 + crhs7*normalslave(1,1);
        const double crhs18 =     -crhs1 + crhs6;
        const double crhs19 =     crhs14 - crhs18*normalslave(1,0);
        const double crhs20 =     crhs16 - crhs18*normalslave(1,1);
    
        rhs[0]+=-crhs0*crhs15;
        rhs[1]+=-crhs0*crhs17;
        rhs[2]+=-crhs15*crhs4;
        rhs[3]+=-crhs17*crhs4;
        rhs[4]+=crhs19*crhs2;
        rhs[5]+=crhs2*crhs20;
        rhs[6]+=crhs19*crhs3;
        rhs[7]+=crhs20*crhs3;
        rhs[10]+=rScaleFactor*(crhs10*tangentetaslave(1,0) + crhs12*tangentxislave(1,0) - crhs5*normalslave(1,0));
        rhs[11]+=rScaleFactor*(crhs10*tangentetaslave(1,1) + crhs12*tangentxislave(1,1) - crhs5*normalslave(1,1));
    }


    return rhs;
}


/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo 
    )
{
    KRATOS_TRY;   
    
    boost::shared_ptr<ConditionMap>& AllConditionSets = this->GetValue( CONTACT_SETS );
    
    // Calculates the size of the system
    const unsigned int ConditionSize = (TDim * ( TNumNodes + TNumNodes + TNumNodes) ) * AllConditionSets->size(); 
    
    if (rResult.size() != ConditionSize)
    {
        rResult.resize( ConditionSize, false );
    }
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    for (auto ipair = AllConditionSets->begin(); ipair != AllConditionSets->end(); ++ipair )
    {
        GeometryType& current_master = (ipair->first)->GetGeometry( );
        
        // Master Nodes Displacement Equation IDs
        for ( unsigned int i_master = 0; i_master < TNumNodes; i_master++ ) // NOTE: Assuming same number of nodes for master and slave
        {
            NodeType& master_node = current_master[i_master];
            rResult[index++] = master_node.GetDof( DISPLACEMENT_X ).EquationId( );
            rResult[index++] = master_node.GetDof( DISPLACEMENT_Y ).EquationId( );
            if (TDim == 3)
            {
                rResult[index++] = master_node.GetDof( DISPLACEMENT_Z ).EquationId( );
            }
        }

        // Slave Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
        {
            NodeType& slave_node = this->GetGeometry()[ i_slave ];
            rResult[index++] = slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
            rResult[index++] = slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
            if (TDim == 3)
            {
                rResult[index++] = slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
            }
        }

        // Slave Nodes  Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
        {
            NodeType& slave_node = this->GetGeometry()[ i_slave ];
            rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_X ).EquationId( );
            rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ).EquationId( );
            if (TDim == 3)
            {
                rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ).EquationId( );
            }
        }
        
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim, TNumNodes>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo 
)
{
    KRATOS_TRY;
    
    boost::shared_ptr<ConditionMap>& AllConditionSets = this->GetValue( CONTACT_SETS );
    
    // Calculates the size of the system
    const unsigned int ConditionSize = (TDim * ( TNumNodes + TNumNodes + TNumNodes) ) * AllConditionSets->size(); 
    
    if (rConditionalDofList.size() != ConditionSize)
    {
        rConditionalDofList.resize( ConditionSize );
    }
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    for (auto ipair = AllConditionSets->begin(); ipair != AllConditionSets->end(); ++ipair )
    {
        GeometryType& current_master = (ipair->first)->GetGeometry( );

        // Master Nodes Displacement Equation IDs
        for ( unsigned int i_master = 0; i_master < TNumNodes; i_master++ ) // NOTE: Assuming same number of nodes for master and slave
        {
            NodeType& master_node = current_master[i_master];
            rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_X );
            rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Y );
            if (TDim == 3)
            {
                rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Z );
            }
        }

        // Slave Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
        {
            NodeType& slave_node = this->GetGeometry()[ i_slave ];
            rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_X );
            rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Y );
            if (TDim == 3)
            {
                rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Z );
            }
        }

        // Slave Nodes Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
        {
            NodeType& slave_node = this->GetGeometry()[ i_slave ];
            rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_X );
            rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Y );
            if (TDim == 3)
            {
                rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Z );
            }
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template class AugmentedLagrangianMethodFrictionalMortarContactCondition<2, 2>;
// template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3>;
// template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4>;

} // Namespace Kratos
