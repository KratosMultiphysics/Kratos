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
#include "custom_elements/updated_lagrangian_U_Pw_diff_order_element.hpp"
#include "custom_utilities/equation_of_motion_utilities.h"
#include "custom_utilities/math_utilities.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "utilities/math_utils.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------
Element::Pointer UpdatedLagrangianUPwDiffOrderElement::Create(IndexType             NewId,
                                                              NodesArrayType const& ThisNodes,
                                                              PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UpdatedLagrangianUPwDiffOrderElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties, this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------
Element::Pointer UpdatedLagrangianUPwDiffOrderElement::Create(IndexType             NewId,
                                                              GeometryType::Pointer pGeom,
                                                              PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UpdatedLagrangianUPwDiffOrderElement(
        NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderElement::CalculateAll(MatrixType&        rLeftHandSideMatrix,
                                                        VectorType&        rRightHandSideVector,
                                                        const ProcessInfo& rCurrentProcessInfo,
                                                        bool CalculateStiffnessMatrixFlag,
                                                        bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const GeometryType&   rGeom = GetGeometry();
    const PropertiesType& rProp = this->GetProperties();

    // Definition of variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    // Stiffness matrix is always needed to calculate Biot coefficient
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag)
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);

    RetentionLaw::Parameters RetentionParameters(rProp);

    // Loop over integration points
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(this->GetIntegrationMethod());

    const auto b_matrices = this->CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
    const auto deformation_gradients = this->CalculateDeformationGradients();
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, Variables.detJuContainer);
    const auto det_Js_initial_configuration =
        GeoEquationOfMotionUtilities::CalculateDetJsInitialConfiguration(rGeom, this->GetIntegrationMethod());
    const auto integration_coefficients_on_initial_configuration =
        this->CalculateIntegrationCoefficients(IntegrationPoints, det_Js_initial_configuration);
    auto strain_vectors = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        this->GetStressStatePolicy().GetVoigtSize());
    std::vector<Matrix> constitutive_matrices;
    this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                         Variables.NuContainer, Variables.DNu_DXContainer,
                                         strain_vectors, mStressVector, constitutive_matrices);

    const auto biot_coefficients = GeoTransportEquationUtilities::CalculateBiotCoefficients(
        constitutive_matrices, this->GetProperties());
    const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
        Variables.NpContainer, Variables.PressureVector);
    const auto relative_permeability_values = CalculateRelativePermeabilityValues(fluid_pressures);
    const auto degrees_of_saturation        = this->CalculateDegreesOfSaturation(fluid_pressures);
    const auto derivatives_of_saturation = this->CalculateDerivativesOfSaturation(fluid_pressures);
    const auto biot_moduli_inverse = GeoTransportEquationUtilities::CalculateInverseBiotModuli(
        biot_coefficients, degrees_of_saturation, derivatives_of_saturation, rProp);
    const auto bishop_coefficients = CalculateBishopCoefficients(fluid_pressures);

    for (IndexType GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {
        this->CalculateKinematics(Variables, GPoint);
        Variables.B                      = b_matrices[GPoint];
        Variables.F                      = deformation_gradients[GPoint];
        Variables.StrainVector           = strain_vectors[GPoint];
        Variables.ConstitutiveMatrix     = constitutive_matrices[GPoint];
        Variables.BiotCoefficient        = biot_coefficients[GPoint];
        Variables.BiotModulusInverse     = biot_moduli_inverse[GPoint];
        Variables.DegreeOfSaturation     = degrees_of_saturation[GPoint];
        Variables.BishopCoefficient      = bishop_coefficients[GPoint];
        Variables.IntegrationCoefficient = integration_coefficients[GPoint];
        Variables.RelativePermeability   = relative_permeability_values[GPoint];
        Variables.IntegrationCoefficientInitialConfiguration =
            integration_coefficients_on_initial_configuration[GPoint];

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
void UpdatedLagrangianUPwDiffOrderElement::CalculateAndAddGeometricStiffnessMatrix(
    MatrixType& rLeftHandSideMatrix, const ElementVariables& rVariables, unsigned int GPoint)
{
    KRATOS_TRY

    const GeometryType& r_geom      = GetGeometry();
    const SizeType      num_U_nodes = r_geom.PointsNumber();
    const SizeType      dimension   = r_geom.WorkingSpaceDimension();

    const auto stress_tensor = MathUtils<double>::StressVectorToTensor(mStressVector[GPoint]);

    const Matrix reduced_Kg_matrix =
        prod(rVariables.DNu_DX,
             rVariables.IntegrationCoefficient * Matrix(prod(stress_tensor, trans(rVariables.DNu_DX)))); // to be optimized

    Matrix geometric_stiffness_matrix = ZeroMatrix(num_U_nodes * dimension, num_U_nodes * dimension);
    MathUtils<double>::ExpandAndAddReducedMatrix(geometric_stiffness_matrix, reduced_Kg_matrix, dimension);

    GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, geometric_stiffness_matrix);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                                                        std::vector<double>& rOutput,
                                                                        const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == REFERENCE_DEFORMATION_GRADIENT_DETERMINANT) {
        rOutput = GeoMechanicsMathUtilities::CalculateDeterminants(this->CalculateDeformationGradients());
    } else {
        SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                        std::vector<Vector>& rOutput,
                                                                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rOutput.resize(this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod()));

    if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
        const auto deformation_gradients = this->CalculateDeformationGradients();

        std::transform(deformation_gradients.begin(), deformation_gradients.end(), rOutput.begin(),
                       [this](const Matrix& rDeformationGradient) {
            return this->CalculateGreenLagrangeStrain(rDeformationGradient);
        });
    } else {
        SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                                        std::vector<Matrix>& rOutput,
                                                                        const ProcessInfo& rCurrentProcessInfo)
{
    rOutput.resize(mConstitutiveLawVector.size());

    if (rVariable == REFERENCE_DEFORMATION_GRADIENT) {
        rOutput = this->CalculateDeformationGradients();
    } else {
        SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

//----------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallStrainUPwDiffOrderElement)
}

//----------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallStrainUPwDiffOrderElement)
}

} // Namespace Kratos
