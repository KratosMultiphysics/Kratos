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
#include "custom_elements/updated_lagrangian_U_Pw_axisymmetric_element.hpp"
#include "axisymmetric_stress_state_policy.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwUpdatedLagrangianAxisymmetricElement<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwUpdatedLagrangianAxisymmetricElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwUpdatedLagrangianAxisymmetricElement<TDim, TNumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwUpdatedLagrangianAxisymmetricElement(NewId, pGeom, pProperties));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwUpdatedLagrangianAxisymmetricElement<TDim, TNumNodes>::CalculateBMatrix(Matrix& rB,
                                                                                const Matrix& GradNpT,
                                                                                const Vector& Np)
{
    KRATOS_TRY

    AxisymmetricStressState stress_state;
    rB = stress_state.CalculateBMatrix(GradNpT, Np, this->GetGeometry());

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
double UPwUpdatedLagrangianAxisymmetricElement<TDim, TNumNodes>::CalculateIntegrationCoefficient(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints, unsigned int PointNumber, double detJ)

{
    AxisymmetricStressState stress_state;
    return stress_state.CalculateIntegrationCoefficient(IntegrationPoints[PointNumber], detJ,
                                                        this->GetGeometry());
}

//----------------------------------------------------------------------------------------------------

template class UPwUpdatedLagrangianAxisymmetricElement<2, 3>;
template class UPwUpdatedLagrangianAxisymmetricElement<2, 4>;

template class UPwUpdatedLagrangianAxisymmetricElement<2, 6>;
template class UPwUpdatedLagrangianAxisymmetricElement<2, 8>;
template class UPwUpdatedLagrangianAxisymmetricElement<2, 9>;
template class UPwUpdatedLagrangianAxisymmetricElement<2, 10>;
template class UPwUpdatedLagrangianAxisymmetricElement<2, 15>;

} // Namespace Kratos
