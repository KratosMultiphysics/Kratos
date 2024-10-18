// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Mohamed Nabi
//                   John van Esch
//

#include "custom_conditions/Pw_point_flux_condition.hpp"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
PwPointFluxCondition<TDim, TNumNodes>::PwPointFluxCondition() : PwCondition<TDim, TNumNodes>()
{
}

template <unsigned int TDim, unsigned int TNumNodes>
PwPointFluxCondition<TDim, TNumNodes>::PwPointFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : PwCondition<TDim, TNumNodes>(NewId, pGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
PwPointFluxCondition<TDim, TNumNodes>::PwPointFluxCondition(IndexType               NewId,
                                                            GeometryType::Pointer   pGeometry,
                                                            PropertiesType::Pointer pProperties)
    : PwCondition<TDim, TNumNodes>(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer PwPointFluxCondition<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                 NodesArrayType const& rThisNodes,
                                                                 PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<PwPointFluxCondition>(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
void PwPointFluxCondition<TDim, TNumNodes>::CalculateRHS(Vector& rRightHandSideVector, const ProcessInfo&)
{
    rRightHandSideVector[0] = this->GetGeometry()[0].FastGetSolutionStepValue(NORMAL_FLUID_FLUX);
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string PwPointFluxCondition<TDim, TNumNodes>::Info() const
{
    return "PwPointFluxCondition";
}

template class PwPointFluxCondition<2, 1>;
template class PwPointFluxCondition<3, 1>;

} // Namespace Kratos
