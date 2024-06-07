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

    SmallStrainUPwDiffOrderElement::CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                                                 CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

    const GeometryType&   rGeom = GetGeometry();

    // Definition of variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Loop over integration points
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(this->GetIntegrationMethod());

    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, Variables.detJuContainer);

    for (IndexType GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {
        this->CalculateKinematics(Variables, GPoint);
        Variables.IntegrationCoefficient = integration_coefficients[GPoint];

        if (CalculateStiffnessMatrixFlag) {
            if (Variables.ConsiderGeometricStiffness)
                this->CalculateAndAddGeometricStiffnessMatrix(rLeftHandSideMatrix, Variables, GPoint);
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderElement::CalculateAndAddGeometricStiffnessMatrix(MatrixType& rLeftHandSideMatrix,
                                                                                   ElementVariables& rVariables,
                                                                                   unsigned int GPoint)
{
    KRATOS_TRY

    const GeometryType& rGeom     = GetGeometry();
    const SizeType      NumUNodes = rGeom.PointsNumber();
    const SizeType      Dim       = rGeom.WorkingSpaceDimension();

    Matrix StressTensor = MathUtils<double>::StressVectorToTensor(mStressVector[GPoint]);

    Matrix ReducedKgMatrix =
        prod(rVariables.DNu_DX,
             rVariables.IntegrationCoefficient * Matrix(prod(StressTensor, trans(rVariables.DNu_DX)))); // to be optimized

    Matrix UUMatrix(NumUNodes * Dim, NumUNodes * Dim);
    noalias(UUMatrix) = ZeroMatrix(NumUNodes * Dim, NumUNodes * Dim);
    MathUtils<double>::ExpandAndAddReducedMatrix(UUMatrix, ReducedKgMatrix, Dim);

    // Distribute stiffness block matrix into the elemental matrix
    GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, UUMatrix);

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

std::vector<double> Kratos::UpdatedLagrangianUPwDiffOrderElement::GetPermeabilityUpdateFactors(
    const std::vector<Vector>&) const
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
