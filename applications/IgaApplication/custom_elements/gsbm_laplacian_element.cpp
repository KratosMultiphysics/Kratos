//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mario Caballero
//

// System includes

// External includes

// Project includes
#include "custom_elements/gsbm_laplacian_element.h"
#include "iga_application_variables.h"

namespace Kratos
{

// The volume fraction is constant over the element's knot span, so weighting
// the test functions (v_alpha = alpha v) scales the whole local system

void GsbmLaplacianElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    rLeftHandSideMatrix *= GetValue(VOLUME_FRACTION);

    KRATOS_CATCH("")
}

void GsbmLaplacianElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    rRightHandSideVector *= GetValue(VOLUME_FRACTION);

    KRATOS_CATCH("")
}

void GsbmLaplacianElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::CalculateMassMatrix(rMassMatrix, rCurrentProcessInfo);
    rMassMatrix *= GetValue(VOLUME_FRACTION);

    KRATOS_CATCH("")
}

int GsbmLaplacianElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(Has(VOLUME_FRACTION)) << "GsbmLaplacianElement #" << Id()
        << ": VOLUME_FRACTION is not set. Run VolumeFractionGsbmProcess first."
        << std::endl;

    return BaseType::Check(rCurrentProcessInfo);
}

} // namespace Kratos
