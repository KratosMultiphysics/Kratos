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
#include "custom_elements/updated_lagrangian_U_Pw_diff_order_element.h"
#include "custom_utilities/math_utilities.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "utilities/math_utils.h"

namespace Kratos
{

Element::Pointer UpdatedLagrangianUPwDiffOrderElement::Create(IndexType             NewId,
                                                              const NodesArrayType& rNodes,
                                                              PropertiesType::Pointer pProperties) const
{
    return Create(NewId, this->GetGeometry().Create(rNodes), pProperties);
}

Element::Pointer UpdatedLagrangianUPwDiffOrderElement::Create(IndexType             NewId,
                                                              GeometryType::Pointer pGeom,
                                                              PropertiesType::Pointer pProperties) const
{
    return make_intrusive<UpdatedLagrangianUPwDiffOrderElement>(
        NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone(),
        this->CloneIntegrationCoefficientModifier());
}

void UpdatedLagrangianUPwDiffOrderElement::CalculateAll(MatrixType&        rLeftHandSideMatrix,
                                                        VectorType&        rRightHandSideVector,
                                                        const ProcessInfo& rCurrentProcessInfo,
                                                        bool CalculateStiffnessMatrixFlag,
                                                        bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    SmallStrainUPwDiffOrderElement::CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
                                                 rCurrentProcessInfo, CalculateStiffnessMatrixFlag,
                                                 CalculateResidualVectorFlag);

    ElementVariables variables;
    this->InitializeElementVariables(variables, rCurrentProcessInfo);

    if (CalculateStiffnessMatrixFlag && variables.ConsiderGeometricStiffness) {
        const auto& integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());
        const auto integration_coefficients =
            this->CalculateIntegrationCoefficients(integration_points, variables.detJuContainer);

        for (IndexType GPoint = 0; GPoint < integration_points.size(); ++GPoint) {
            this->CalculateAndAddGeometricStiffnessMatrix(rLeftHandSideMatrix, mStressVector[GPoint],
                                                          variables.DNu_DXContainer[GPoint],
                                                          integration_coefficients[GPoint]);
        }
    }

    KRATOS_CATCH("")
}

void UpdatedLagrangianUPwDiffOrderElement::CalculateAndAddGeometricStiffnessMatrix(
    MatrixType& rLeftHandSideMatrix, const Vector& rStressVector, const Matrix& rDNuDx, const double IntegrationCoefficient) const
{
    KRATOS_TRY

    const GeometryType& r_geom      = GetGeometry();
    const SizeType      num_U_nodes = r_geom.PointsNumber();
    const SizeType      dimension   = r_geom.WorkingSpaceDimension();

    const Matrix reduced_Kg_matrix =
        prod(rDNuDx, Matrix(prod(MathUtils<double>::StressVectorToTensor(rStressVector), trans(rDNuDx)))) *
        IntegrationCoefficient;

    Matrix geometric_stiffness_matrix = ZeroMatrix(num_U_nodes * dimension, num_U_nodes * dimension);
    MathUtils<double>::ExpandAndAddReducedMatrix(geometric_stiffness_matrix, reduced_Kg_matrix, dimension);

    GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, geometric_stiffness_matrix);

    KRATOS_CATCH("")
}

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

std::vector<double> Kratos::UpdatedLagrangianUPwDiffOrderElement::GetOptionalPermeabilityUpdateFactors(
    const std::vector<Vector>&) const
{
    return {};
}

void UpdatedLagrangianUPwDiffOrderElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallStrainUPwDiffOrderElement)
}

void UpdatedLagrangianUPwDiffOrderElement::load(Serializer& rSerializer){
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallStrainUPwDiffOrderElement)}

std::string UpdatedLagrangianUPwDiffOrderElement::Info() const
{
    const std::string constitutive_info =
        !mConstitutiveLawVector.empty() ? mConstitutiveLawVector[0]->Info() : "not defined";
    return "Updated Lagrangian U-Pw different order Element #" + std::to_string(this->Id()) +
           "\nConstitutive law: " + constitutive_info;
}

void UpdatedLagrangianUPwDiffOrderElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

void UpdatedLagrangianUPwDiffOrderElement::PrintData(std::ostream& rOStream) const
{
    this->pGetGeometry()->PrintData(rOStream);
}
} // Namespace Kratos
