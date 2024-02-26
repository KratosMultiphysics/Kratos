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
#include "custom_elements/U_Pw_small_strain_axisymmetric_FIC_element.hpp"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainAxisymmetricFICElement<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwSmallStrainAxisymmetricFICElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainAxisymmetricFICElement<TDim, TNumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwSmallStrainAxisymmetricFICElement(NewId, pGeom, pProperties));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainAxisymmetricFICElement<TDim, TNumNodes>::CalculateBMatrix(Matrix&       rB,
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
double UPwSmallStrainAxisymmetricFICElement<TDim, TNumNodes>::CalculateIntegrationCoefficient(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints, unsigned int PointNumber, double detJ)

{
    AxisymmetricStressState stress_state;
    return stress_state.CalculateIntegrationCoefficient(IntegrationPoints[PointNumber], detJ, this->GetGeometry());
}

template <unsigned int TDim, unsigned int TNumNodes>
Vector UPwSmallStrainAxisymmetricFICElement<TDim, TNumNodes>::CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient)
{
    AxisymmetricStressState stress_state;
    return stress_state.CalculateGreenLagrangeStrain(rDeformationGradient);
}

//----------------------------------------------------------------------------------------------------

template class UPwSmallStrainAxisymmetricFICElement<2, 3>;
template class UPwSmallStrainAxisymmetricFICElement<2, 4>;

} // Namespace Kratos
