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

    SmallStrainUPwDiffOrderElement::CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
                                                 rCurrentProcessInfo, CalculateStiffnessMatrixFlag,
                                                 CalculateResidualVectorFlag);

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    if (CalculateStiffnessMatrixFlag && Variables.ConsiderGeometricStiffness) {
        const GeometryType& rGeom = GetGeometry();

        const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
            rGeom.IntegrationPoints(this->GetIntegrationMethod());
        const auto integration_coefficients =
            this->CalculateIntegrationCoefficients(IntegrationPoints, Variables.detJuContainer);

        for (IndexType GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {
            this->CalculateAndAddGeometricStiffnessMatrix(rLeftHandSideMatrix, GPoint,
                                                          Variables.DNu_DXContainer[GPoint],
                                                          integration_coefficients[GPoint]);
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderElement::CalculateAndAddGeometricStiffnessMatrix(
    MatrixType& rLeftHandSideMatrix, unsigned int GPoint, const Matrix& rDNuDx, const double IntegrationCoefficient)
{
    KRATOS_TRY

    const GeometryType& r_geom      = GetGeometry();
    const SizeType      num_U_nodes = r_geom.PointsNumber();
    const SizeType      dimension   = r_geom.WorkingSpaceDimension();

    const auto stress_tensor = MathUtils<double>::StressVectorToTensor(mStressVector[GPoint]);

    const Matrix reduced_Kg_matrix =
        prod(rDNuDx, IntegrationCoefficient * Matrix(prod(stress_tensor, trans(rDNuDx))));

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

std::vector<double> Kratos::UpdatedLagrangianUPwDiffOrderElement::GetPermeabilityUpdateFactors(const std::vector<Vector>&) const
{
    return {};
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
