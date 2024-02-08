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
#include "custom_elements/small_strain_U_Pw_diff_order_axisymmetric_element.hpp"
#include "axisymmetric_stress_state_policy.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------
Element::Pointer SmallStrainUPwDiffOrderAxisymmetricElement::Create(IndexType             NewId,
                                                                    NodesArrayType const& ThisNodes,
                                                                    PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SmallStrainUPwDiffOrderAxisymmetricElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------
Element::Pointer SmallStrainUPwDiffOrderAxisymmetricElement::Create(IndexType             NewId,
                                                                    GeometryType::Pointer pGeom,
                                                                    PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SmallStrainUPwDiffOrderAxisymmetricElement(NewId, pGeom, pProperties));
}

//----------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderAxisymmetricElement::CalculateBMatrix(Matrix& rB, const Matrix& GradNpT, const Vector& Np)
{
    KRATOS_TRY

    AxisymmetricStressState stress_state;
    rB = stress_state.CalculateBMatrix(GradNpT, Np, this->GetGeometry());

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
double SmallStrainUPwDiffOrderAxisymmetricElement::CalculateIntegrationCoefficient(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints, unsigned int PointNumber, double detJ)

{
    AxisymmetricStressState stress_state;
    return stress_state.CalculateIntegrationCoefficient(IntegrationPoints[PointNumber], detJ, this->GetGeometry());
}

void SmallStrainUPwDiffOrderAxisymmetricElement::CalculateCauchyGreenStrain(SmallStrainUPwDiffOrderElement::ElementVariables& rVariables)
{
    AxisymmetricStressState stress_state;
    stress_state.CalculateGreenLagrangeStrain(rVariables.F);
}

//----------------------------------------------------------------------------------------------------

} // Namespace Kratos
