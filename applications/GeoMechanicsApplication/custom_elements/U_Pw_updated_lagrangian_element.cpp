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
#include "utilities/math_utils.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwUpdatedLagrangianElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                      NodesArrayType const& ThisNodes,
                                                                      PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwUpdatedLagrangianElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties, this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwUpdatedLagrangianElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                      GeometryType::Pointer pGeom,
                                                                      PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwUpdatedLagrangianElement(NewId, pGeom, pProperties,
                                                            this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwUpdatedLagrangianElement<TDim, TNumNodes>::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                                                VectorType& rRightHandSideVector,
                                                                const ProcessInfo& rCurrentProcessInfo,
                                                                bool CalculateStiffnessMatrixFlag,
                                                                bool CalculateResidualVectorFlag)
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        this->GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    // Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(this->GetGeometry(), this->GetProperties(),
                                                       rCurrentProcessInfo);

    // Stiffness matrix is always needed for Biot coefficient
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag) ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(this->GetProperties(), rCurrentProcessInfo);

    const bool hasBiotCoefficient = this->GetProperties().Has(BIOT_COEFFICIENT);
    const auto b_matrices = this->CalculateBMatrices(Variables.DN_DXContainer, Variables.NContainer);
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, Variables.detJContainer);
    const auto deformation_gradients = this->CalculateDeformationGradients();
    const auto determinants_of_deformation_gradients =
        GeoMechanicsMathUtilities::CalculateDeterminants(deformation_gradients);
    const auto strain_vectors = this->CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain);

    // Computing in all integrations points
    for (IndexType GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {
        // Compute element kinematics B, F, GradNpT ...
        this->CalculateKinematics(Variables, GPoint);
        Variables.B = b_matrices[GPoint];

        // Cauchy strain: This needs to be investigated which strain measure should be used
        // In some references, e.g. Bathe, suggested to use Almansi strain measure
        Variables.F            = deformation_gradients[GPoint];
        Variables.detF         = determinants_of_deformation_gradients[GPoint];
        Variables.StrainVector = strain_vectors[GPoint];

        // set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // Compute Np, Nu and BodyAcceleration
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, Variables.NContainer, GPoint);
        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, Variables.NContainer, Variables.VolumeAcceleration, GPoint);

        // Compute constitutive tensor and stresses
        ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        this->CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        // calculate Bulk modulus from stiffness matrix
        this->InitializeBiotCoefficients(Variables, hasBiotCoefficient);

        Variables.IntegrationCoefficient = integration_coefficients[GPoint];

        Variables.IntegrationCoefficientInitialConfiguration = this->CalculateIntegrationCoefficient(
            IntegrationPoints[GPoint], Variables.detJInitialConfiguration);

        if (CalculateStiffnessMatrixFlag) {
            // Contributions to stiffness matrix calculated on the reference config
            /* Material stiffness matrix */
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

            /* Geometric stiffness matrix */
            if (Variables.ConsiderGeometricStiffness)
                this->CalculateAndAddGeometricStiffnessMatrix(rLeftHandSideMatrix, Variables, GPoint);
        }

        if (CalculateResidualVectorFlag) {
            // Contributions to the right hand side
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
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

//----------------------------------------------------------------------------------------
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

//----------------------------------------------------------------------------------------

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
