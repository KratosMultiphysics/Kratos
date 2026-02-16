// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Alejandro Cornejo
//

// Project includes
#include "non_linear_timoshenko_beam_element_2D2N.h"

namespace Kratos
{

Element::Pointer NonLinearTimoshenkoBeamElement2D2N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    NonLinearTimoshenkoBeamElement2D2N::Pointer p_new_elem = Kratos::make_intrusive<NonLinearTimoshenkoBeamElement2D2N>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::GetShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi) const
{
    VectorType local_N(4);
    BaseType::GetShapeFunctionsValues(local_N, Length, Phi, xi);
    GlobalSizeVector(rN, local_N);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::GetFirstDerivativesShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi) const
{
    VectorType local_N(4);
    BaseType::GetFirstDerivativesShapeFunctionsValues(local_N, Length, Phi, xi);
    GlobalSizeVector(rN, local_N);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::GetSecondDerivativesShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi) const
{
    VectorType local_N(4);
    BaseType::GetSecondDerivativesShapeFunctionsValues(local_N, Length, Phi, xi);
    GlobalSizeVector(rN, local_N);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::GetThirdDerivativesShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi) const
{
    VectorType local_N(4);
    BaseType::GetThirdDerivativesShapeFunctionsValues(local_N, Length, Phi, xi);
    GlobalSizeVector(rN, local_N);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::GetFourthDerivativesShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi) const
{
    VectorType local_N(4);
    BaseType::GetFourthDerivativesShapeFunctionsValues(local_N, Length, Phi, xi);
    GlobalSizeVector(rN, local_N);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::GetNThetaShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi) const
{
    VectorType local_N(4);
    BaseType::GetNThetaShapeFunctionsValues(local_N, Length, Phi, xi);
    GlobalSizeVector(rN, local_N);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::GetFirstDerivativesNThetaShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi) const
{
    VectorType local_N(4);
    BaseType::GetFirstDerivativesNThetaShapeFunctionsValues(local_N, Length, Phi, xi);
    GlobalSizeVector(rN, local_N);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::GetNu0ShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi) const
{
    VectorType local_N(2);
    BaseType::GetNu0ShapeFunctionsValues(local_N, Length, Phi, xi);
    GlobalSizeAxialVector(rN, local_N);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::GetFirstDerivativesNu0ShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi) const
{
    VectorType local_N(2);
    BaseType::GetFirstDerivativesNu0ShapeFunctionsValues(local_N, Length, Phi, xi);
    GlobalSizeAxialVector(rN, local_N);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    VectorType dummy_rhs;
    CalculateAll(rLeftHandSideMatrix, dummy_rhs, rCurrentProcessInfo, true, false);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    MatrixType dummy_lhs;
    CalculateAll(dummy_lhs, rRightHandSideVector, rCurrentProcessInfo, false, true);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLHS,
    const bool ComputeRHS
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

} // namespace Kratos
