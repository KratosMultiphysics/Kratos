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
#include "custom_utilities/equation_of_motion_utilities.hpp"
#include "custom_utilities/math_utilities.hpp"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "utilities/math_utils.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwUpdatedLagrangianFICElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                         NodesArrayType const& ThisNodes,
                                                                         PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwUpdatedLagrangianFICElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties,
        this->GetStressStatePolicy().Clone(), this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwUpdatedLagrangianFICElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                         GeometryType::Pointer pGeom,
                                                                         PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwUpdatedLagrangianFICElement(
        NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone(),
        this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwUpdatedLagrangianFICElement<TDim, TNumNodes>::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                                                   VectorType& rRightHandSideVector,
                                                                   const ProcessInfo& rCurrentProcessInfo,
                                                                   bool CalculateStiffnessMatrixFlag,
                                                                   bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

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

    RetentionLaw::Parameters RetentionParameters(this->GetProperties());

    const auto b_matrices = this->CalculateBMatrices(Variables.DN_DXContainer, Variables.NContainer);
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, Variables.detJContainer);
    const auto det_Js_initial_configuration =
        GeoEquationOfMotionUtilities::CalculateDetJsInitialConfiguration(Geom, this->GetIntegrationMethod());
    const auto integration_coefficients_on_initial_configuration =
        this->CalculateIntegrationCoefficients(IntegrationPoints, det_Js_initial_configuration);

    const auto deformation_gradients = this->CalculateDeformationGradients();
    auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        this->GetStressStatePolicy().GetVoigtSize());
    std::vector<Matrix> constitutive_matrices;
    this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                         Variables.NContainer, Variables.DN_DXContainer,
                                         strain_vectors, mStressVector, constitutive_matrices);
    const auto biot_coefficients =
        GeoTransportEquationUtilities::CalculateBiotCoefficients(constitutive_matrices, Prop);
    const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
        Variables.NContainer, Variables.PressureVector);
    const auto relative_permeability_values = this->CalculateRelativePermeabilityValues(fluid_pressures);
    const auto degrees_of_saturation     = this->CalculateDegreesOfSaturation(fluid_pressures);
    const auto derivatives_of_saturation = this->CalculateDerivativesOfSaturation(fluid_pressures);
    const auto biot_moduli_inverse = GeoTransportEquationUtilities::CalculateInverseBiotModuli(
        biot_coefficients, degrees_of_saturation, derivatives_of_saturation, Prop);
    const auto bishop_coefficients = this->CalculateBishopCoefficients(fluid_pressures);

    // Computing in all integrations points
    for (IndexType GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {
        // Compute element kinematics B, F, GradNpT ...
        this->CalculateKinematics(Variables, GPoint);
        Variables.B = b_matrices[GPoint];

        // Cauchy strain: This needs to be investigated which strain measure should be used
        // In some references, e.g. Bathe, suggested to use Almansi strain measure
        Variables.F                  = deformation_gradients[GPoint];
        Variables.StrainVector       = strain_vectors[GPoint];
        Variables.ConstitutiveMatrix = constitutive_matrices[GPoint];

        // Compute Np, Nu and BodyAcceleration
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, Variables.NContainer, GPoint);

        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, Variables.NContainer, Variables.VolumeAcceleration, GPoint);

        // Compute ShapeFunctionsSecondOrderGradients
        this->CalculateShapeFunctionsSecondOrderGradients(FICVariables, Variables);

        Variables.RelativePermeability = relative_permeability_values[GPoint];
        Variables.BishopCoefficient    = bishop_coefficients[GPoint];

        // set shear modulus from stiffness matrix
        FICVariables.ShearModulus = CalculateShearModulus(Variables.ConstitutiveMatrix);

        Variables.BiotCoefficient    = biot_coefficients[GPoint];
        Variables.BiotModulusInverse = biot_moduli_inverse[GPoint];
        Variables.DegreeOfSaturation = degrees_of_saturation[GPoint];

        Variables.IntegrationCoefficient = integration_coefficients[GPoint];

        Variables.IntegrationCoefficientInitialConfiguration =
            integration_coefficients_on_initial_configuration[GPoint];

        if (CalculateStiffnessMatrixFlag) {
            // Contributions to stiffness matrix calculated on the reference config
            /* Material stiffness matrix */
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);
            this->CalculateAndAddLHSStabilization(rLeftHandSideMatrix, Variables, FICVariables);

            if (Variables.ConsiderGeometricStiffness)
                this->CalculateAndAddGeometricStiffnessMatrix(
                    rLeftHandSideMatrix, this->mStressVector[GPoint], Variables.GradNpT,
                    Variables.IntegrationCoefficient);
        }

        if (CalculateResidualVectorFlag) {
            // Contributions to the right hand side
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
            this->CalculateAndAddRHSStabilization(rRightHandSideVector, Variables, FICVariables);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwUpdatedLagrangianFICElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == REFERENCE_DEFORMATION_GRADIENT_DETERMINANT) {
        rOutput = GeoMechanicsMathUtilities::CalculateDeterminants(this->CalculateDeformationGradients());
    } else {
        UPwSmallStrainFICElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
            rVariable, rOutput, rCurrentProcessInfo);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwUpdatedLagrangianFICElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    rOutput.resize(this->GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod));

    if (rVariable == REFERENCE_DEFORMATION_GRADIENT) {
        rOutput = this->CalculateDeformationGradients();
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR) {
        const auto deformation_gradients = this->CalculateDeformationGradients();
        std::ranges::transform(deformation_gradients, rOutput.begin(), [this](const Matrix& rDeformationGradient) {
            return MathUtils<>::StrainVectorToTensor(this->CalculateGreenLagrangeStrain(rDeformationGradient));
        });
    } else {
        UPwSmallStrainFICElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
            rVariable, rOutput, rCurrentProcessInfo);
    }
}

template class UPwUpdatedLagrangianFICElement<2, 3>;
template class UPwUpdatedLagrangianFICElement<2, 4>;
template class UPwUpdatedLagrangianFICElement<3, 4>;
template class UPwUpdatedLagrangianFICElement<3, 8>;

} // Namespace Kratos
