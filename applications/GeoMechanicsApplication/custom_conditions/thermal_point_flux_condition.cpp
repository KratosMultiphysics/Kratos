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

#include "custom_conditions/thermal_point_flux_condition.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
GeoThermalPointFluxCondition<TDim, TNumNodes>::GeoThermalPointFluxCondition()
    : GeoTCondition<TDim, TNumNodes>()
{
}

template <unsigned int TDim, unsigned int TNumNodes>
GeoThermalPointFluxCondition<TDim, TNumNodes>::GeoThermalPointFluxCondition(IndexType NewId,
                                                                            GeometryType::Pointer pGeometry)
    : GeoTCondition<TDim, TNumNodes>(NewId, pGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
GeoThermalPointFluxCondition<TDim, TNumNodes>::GeoThermalPointFluxCondition(IndexType NewId,
                                                                            GeometryType::Pointer pGeometry,
                                                                            PropertiesType::Pointer pProperties)
    : GeoTCondition<TDim, TNumNodes>(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer GeoThermalPointFluxCondition<TDim, TNumNodes>::Create(IndexType NewId,
                                                                         const NodesArrayType& rThisNodes,
                                                                         PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoThermalPointFluxCondition>(
        NewId, this->GetGeometry().Create(rThisNodes), pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoThermalPointFluxCondition<TDim, TNumNodes>::CalculateRHS(Vector& rRightHandSideVector, const ProcessInfo&)
{
    rRightHandSideVector[0] = this->GetGeometry()[0].FastGetSolutionStepValue(NORMAL_HEAT_FLUX);
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string GeoThermalPointFluxCondition<TDim, TNumNodes>::Info() const
{
    return "GeoThermalPointFluxCondition";
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoThermalPointFluxCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoThermalPointFluxCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
}

template class GeoThermalPointFluxCondition<2, 1>;
template class GeoThermalPointFluxCondition<3, 1>;

} // namespace Kratos
