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

#pragma once

#include "custom_conditions/T_condition.h"
#include "includes/serializer.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoThermalPointFluxCondition
    : public GeoTCondition<TDim, TNumNodes>
{
public:
    using GeometryType   = Geometry<Node>;
    using PropertiesType = Properties;
    using NodesArrayType = GeometryType::PointsArrayType;
    using BaseType       = GeoTCondition<TDim, TNumNodes>;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoThermalPointFluxCondition);

    GeoThermalPointFluxCondition();

    GeoThermalPointFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    GeoThermalPointFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    Condition::Pointer Create(IndexType NewId, const NodesArrayType& rThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GeoThermalPointFluxCondition>(
            NewId, this->GetGeometry().Create(rThisNodes), pProperties);
    }

    std::string Info() const override;

protected:
    void CalculateRHS(Vector& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo) override;

private:
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    }
};

} // namespace Kratos
