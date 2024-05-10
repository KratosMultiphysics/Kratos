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
#include "custom_elements/U_Pw_updated_lagrangian_FIC_element.hpp"
#include "utilities/math_utils.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwUpdatedLagrangianFICElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                         NodesArrayType const& ThisNodes,
                                                                         PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwUpdatedLagrangianFICElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties, this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwUpdatedLagrangianFICElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                         GeometryType::Pointer pGeom,
                                                                         PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwUpdatedLagrangianFICElement(
        NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwUpdatedLagrangianFICElement<TDim, TNumNodes>::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                                                   VectorType& rRightHandSideVector,
                                                                   const ProcessInfo& rCurrentProcessInfo,
                                                                   bool CalculateStiffnessMatrixFlag,
                                                                   bool CalculateResidualVectorFlag)
{
    KRATOS_TRY;

    const PropertiesType&                           Prop = this->GetProperties();
    const GeometryType&                             Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        Geom.IntegrationPoints(mThisIntegrationMethod);

    // Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom, Prop, rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag) ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    // Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    FICElementVariables FICVariables;
    this->InitializeFICElementVariables(FICVariables, Variables.DN_DXContainer, Geom, Prop, rCurrentProcessInfo);

    // create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(this->GetProperties(), rCurrentProcessInfo);

    const bool hasBiotCoefficient = Prop.Has(BIOT_COEFFICIENT);

    const auto b_matrices = this->CalculateBMatrices(Variables.DN_DXContainer, Variables.NContainer);

    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, Variables.detJContainer);

    // Computing in all integrations points
    for (IndexType GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {
        // Compute element kinematics B, F, GradNpT ...
        this->CalculateKinematics(Variables, GPoint);
        Variables.B = b_matrices[GPoint];

        // Cauchy strain: This needs to be investigated which strain measure should be used
        // In some references, e.g. Bathe, suggested to use Almansi strain measure
        Variables.F            = this->CalculateDeformationGradient(GPoint);
        Variables.detF         = MathUtils<>::Det(Variables.F);
        Variables.StrainVector = this->CalculateStrain(
            Variables.F, Variables.B, Variables.DisplacementVector, Variables.UseHenckyStrain);

        // set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // Compute Np, Nu and BodyAcceleration
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, Variables.NContainer, GPoint);

        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, Variables.NContainer, Variables.VolumeAcceleration, GPoint);

        // Compute ShapeFunctionsSecondOrderGradients
        this->CalculateShapeFunctionsSecondOrderGradients(FICVariables, Variables);

        // Compute constitutive tensor and stresses
        ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        this->CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        // set shear modulus from stiffness matrix
        FICVariables.ShearModulus = CalculateShearModulus(Variables.ConstitutiveMatrix);

        // calculate Bulk modulus from stiffness matrix
        this->InitializeBiotCoefficients(Variables, hasBiotCoefficient);

        Variables.IntegrationCoefficient = integration_coefficients[GPoint];

        Variables.IntegrationCoefficientInitialConfiguration = this->CalculateIntegrationCoefficient(
            IntegrationPoints[GPoint], Variables.detJInitialConfiguration);

        if (CalculateStiffnessMatrixFlag) {
            // Contributions to stiffness matrix calculated on the reference config
            /* Material stiffness matrix */
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);
            this->CalculateAndAddLHSStabilization(rLeftHandSideMatrix, Variables, FICVariables);

            /* Geometric stiffness matrix */
            if (Variables.ConsiderGeometricStiffness)
                this->CalculateAndAddGeometricStiffnessMatrix(rLeftHandSideMatrix, Variables, GPoint);
        }

        if (CalculateResidualVectorFlag) {
            // Contributions to the right hand side
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
            this->CalculateAndAddRHSStabilization(rRightHandSideVector, Variables, FICVariables);
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwUpdatedLagrangianFICElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == REFERENCE_DEFORMATION_GRADIENT_DETERMINANT) {
        rOutput.clear();
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput.emplace_back(MathUtils<>::Det(this->CalculateDeformationGradient(GPoint)));
        }
    } else {
        UPwSmallStrainFICElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
            rVariable, rOutput, rCurrentProcessInfo);
    }
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwUpdatedLagrangianFICElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    rOutput.resize(this->GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod));

    if (rVariable == REFERENCE_DEFORMATION_GRADIENT) {
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = this->CalculateDeformationGradient(GPoint);
        }
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR) {
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = MathUtils<>::StrainVectorToTensor(
                this->CalculateGreenLagrangeStrain(this->CalculateDeformationGradient(GPoint)));
        }
    } else {
        UPwSmallStrainFICElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
            rVariable, rOutput, rCurrentProcessInfo);
    }
}

//----------------------------------------------------------------------------------------

template class UPwUpdatedLagrangianFICElement<2, 3>;
template class UPwUpdatedLagrangianFICElement<2, 4>;
template class UPwUpdatedLagrangianFICElement<3, 4>;
template class UPwUpdatedLagrangianFICElement<3, 8>;

} // Namespace Kratos
