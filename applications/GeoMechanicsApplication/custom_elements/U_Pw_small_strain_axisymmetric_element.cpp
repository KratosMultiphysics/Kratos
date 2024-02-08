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
#include "custom_elements/U_Pw_small_strain_axisymmetric_element.hpp"
#include "axisymmetric_stress_state_policy.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainAxisymmetricElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                            NodesArrayType const& ThisNodes,
                                                                            PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwSmallStrainAxisymmetricElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainAxisymmetricElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                            GeometryType::Pointer pGeom,
                                                                            PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwSmallStrainAxisymmetricElement(NewId, pGeom, pProperties));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainAxisymmetricElement<TDim, TNumNodes>::CalculateBMatrix(Matrix&       rB,
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
double UPwSmallStrainAxisymmetricElement<TDim, TNumNodes>::CalculateIntegrationCoefficient(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints, unsigned int PointNumber, double detJ)

{
    AxisymmetricStressState stress_state;
    return stress_state.CalculateIntegrationCoefficient(IntegrationPoints[PointNumber], detJ,
                                                        this->GetGeometry());
}

//----------------------------------------------------------------------------------------------------

template class UPwSmallStrainAxisymmetricElement<2, 3>;
template class UPwSmallStrainAxisymmetricElement<2, 4>;

template class UPwSmallStrainAxisymmetricElement<2, 6>;
template class UPwSmallStrainAxisymmetricElement<2, 8>;
template class UPwSmallStrainAxisymmetricElement<2, 9>;
template class UPwSmallStrainAxisymmetricElement<2, 10>;
template class UPwSmallStrainAxisymmetricElement<2, 15>;

} // Namespace Kratos
