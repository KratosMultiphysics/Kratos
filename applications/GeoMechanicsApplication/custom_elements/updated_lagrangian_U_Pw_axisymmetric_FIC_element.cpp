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
#include "custom_elements/updated_lagrangian_U_Pw_axisymmetric_FIC_element.hpp"
#include "axisymmetric_stress_state.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwUpdatedLagrangianAxisymmetricFICElement<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwUpdatedLagrangianAxisymmetricFICElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwUpdatedLagrangianAxisymmetricFICElement<TDim, TNumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwUpdatedLagrangianAxisymmetricFICElement(NewId, pGeom, pProperties));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwUpdatedLagrangianAxisymmetricFICElement<TDim, TNumNodes>::CalculateBMatrix(Matrix& rB,
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
double UPwUpdatedLagrangianAxisymmetricFICElement<TDim, TNumNodes>::CalculateIntegrationCoefficient(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints, unsigned int PointNumber, double detJ)

{
    AxisymmetricStressState stress_state;
    return stress_state.CalculateIntegrationCoefficient(IntegrationPoints[PointNumber], detJ, this->GetGeometry());
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwUpdatedLagrangianAxisymmetricFICElement<TDim, TNumNodes>::CalculateGreenLagrangeStrain(ElementVariables& rVariables)
{
    AxisymmetricStressState stress_state;
    stress_state.CalculateGreenLagrangeStrain(rVariables.F);
}

//----------------------------------------------------------------------------------------------------

template class UPwUpdatedLagrangianAxisymmetricFICElement<2, 3>;
template class UPwUpdatedLagrangianAxisymmetricFICElement<2, 4>;

} // Namespace Kratos
