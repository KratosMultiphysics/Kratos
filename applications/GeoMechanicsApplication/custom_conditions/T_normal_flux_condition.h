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
#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "includes/serializer.h"

namespace Kratos {

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoTNormalFluxCondition
    : public GeoTCondition<TDim, TNumNodes> {
public:
    using GeometryType = Geometry<Node>;
    using PropertiesType = Properties;
    using NodesArrayType = GeometryType::PointsArrayType;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoTNormalFluxCondition);

    GeoTNormalFluxCondition();

    GeoTNormalFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    GeoTNormalFluxCondition(IndexType NewId,
                            GeometryType::Pointer pGeometry,
                            PropertiesType::Pointer pProperties);

    ~GeoTNormalFluxCondition() override;

    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& rThisNodes,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GeoTNormalFluxCondition>(
            NewId, this->GetGeometry().Create(rThisNodes), pProperties);
    }

protected:
    void CalculateRHS(Vector& rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo) override;

private:
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
    }
};

} // namespace Kratos
