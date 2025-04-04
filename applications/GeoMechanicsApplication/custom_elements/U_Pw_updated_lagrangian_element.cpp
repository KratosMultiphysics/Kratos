// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// External includes

// Project includes
#include "custom_elements/U_Pw_updated_lagrangian_element.hpp"
#include "custom_utilities/math_utilities.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "utilities/math_utils.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwUpdatedLagrangianElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                      NodesArrayType const& ThisNodes,
                                                                      PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwUpdatedLagrangianElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties,
        this->GetStressStatePolicy().Clone(), this->GetIntegrationCoefficientsCalculator().Clone()));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwUpdatedLagrangianElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                      GeometryType::Pointer pGeom,
                                                                      PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwUpdatedLagrangianElement(
        NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone(),
        this->GetIntegrationCoefficientsCalculator().Clone()));
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwUpdatedLagrangianElement<TDim, TNumNodes>::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                                                VectorType& rRightHandSideVector,
                                                                const ProcessInfo& rCurrentProcessInfo,
                                                                bool CalculateStiffnessMatrixFlag,
                                                                bool CalculateResidualVectorFlag)
{
    KRATOS_TRY;

    UPwSmallStrainElement<TDim, TNumNodes>::CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
                                                         rCurrentProcessInfo, CalculateStiffnessMatrixFlag,
                                                         CalculateResidualVectorFlag);
    ElementVariables variables;
    this->InitializeElementVariables(variables, rCurrentProcessInfo);

    if (CalculateStiffnessMatrixFlag && variables.ConsiderGeometricStiffness) {
        const auto& integration_points = this->GetGeometry().IntegrationPoints(mThisIntegrationMethod);
        const auto integration_coefficients =
            this->CalculateIntegrationCoefficients(integration_points, variables.detJContainer);
        for (IndexType GPoint = 0; GPoint < integration_points.size(); ++GPoint) {
            this->CalculateAndAddGeometricStiffnessMatrix(
                rLeftHandSideMatrix, this->mStressVector[GPoint], variables.DN_DXContainer[GPoint],
                integration_coefficients[GPoint]);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwUpdatedLagrangianElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == REFERENCE_DEFORMATION_GRADIENT_DETERMINANT) {
        rOutput = GeoMechanicsMathUtilities::CalculateDeterminants(this->CalculateDeformationGradients());
    } else {
        UPwSmallStrainElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwUpdatedLagrangianElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    rOutput.resize(this->GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod));

    if (rVariable == REFERENCE_DEFORMATION_GRADIENT) {
        rOutput = this->CalculateDeformationGradients();
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR) {
        const auto deformation_gradients = this->CalculateDeformationGradients();
        std::transform(deformation_gradients.begin(), deformation_gradients.end(), rOutput.begin(),
                       [this](const Matrix& rDeformationGradient) {
            return MathUtils<>::StrainVectorToTensor(this->CalculateGreenLagrangeStrain(rDeformationGradient));
        });
    } else {
        UPwSmallStrainElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
std::vector<double> UPwUpdatedLagrangianElement<TDim, TNumNodes>::GetOptionalPermeabilityUpdateFactors(
    const std::vector<Vector>&) const
{
    return {};
}

template class UPwUpdatedLagrangianElement<2, 3>;
template class UPwUpdatedLagrangianElement<2, 4>;
template class UPwUpdatedLagrangianElement<3, 4>;
template class UPwUpdatedLagrangianElement<3, 8>;

template class UPwUpdatedLagrangianElement<2, 6>;
template class UPwUpdatedLagrangianElement<2, 8>;
template class UPwUpdatedLagrangianElement<2, 9>;
template class UPwUpdatedLagrangianElement<2, 10>;
template class UPwUpdatedLagrangianElement<2, 15>;
template class UPwUpdatedLagrangianElement<3, 10>;
template class UPwUpdatedLagrangianElement<3, 20>;
template class UPwUpdatedLagrangianElement<3, 27>;

} // Namespace Kratos
