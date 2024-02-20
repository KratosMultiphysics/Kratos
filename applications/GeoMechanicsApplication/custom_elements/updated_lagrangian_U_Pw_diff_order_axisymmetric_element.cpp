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

// Application includes
#include "custom_elements/updated_lagrangian_U_Pw_diff_order_axisymmetric_element.hpp"

namespace Kratos
{

//----------------------------------------------------------------------------------------
Element::Pointer UpdatedLagrangianUPwDiffOrderAxisymmetricElement::Create(IndexType NewId,
                                                                          NodesArrayType const& ThisNodes,
                                                                          PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UpdatedLagrangianUPwDiffOrderAxisymmetricElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------
Element::Pointer UpdatedLagrangianUPwDiffOrderAxisymmetricElement::Create(IndexType NewId,
                                                                          GeometryType::Pointer pGeom,
                                                                          PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UpdatedLagrangianUPwDiffOrderAxisymmetricElement(NewId, pGeom, pProperties));
}

//----------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderAxisymmetricElement::CalculateBMatrix(Matrix&       rB,
                                                                        const Matrix& GradNpT,
                                                                        const Vector& Np)
{
    KRATOS_TRY

    AxisymmetricStressState stress_state;
    rB = stress_state.CalculateBMatrix(GradNpT, Np, this->GetGeometry());

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
double UpdatedLagrangianUPwDiffOrderAxisymmetricElement::CalculateIntegrationCoefficient(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints, unsigned int PointNumber, double detJ)

{
    AxisymmetricStressState stress_state;
    return stress_state.CalculateIntegrationCoefficient(IntegrationPoints[PointNumber], detJ, this->GetGeometry());
}

Vector UpdatedLagrangianUPwDiffOrderAxisymmetricElement::CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient)
{
    AxisymmetricStressState stress_state;
    return stress_state.CalculateGreenLagrangeStrain(rDeformationGradient);
}

//----------------------------------------------------------------------------------------------------

} // Namespace Kratos
