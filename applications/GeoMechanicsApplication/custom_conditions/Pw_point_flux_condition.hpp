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

#include "custom_conditions/Pw_condition.hpp"
#include "includes/serializer.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION) PwPointFluxCondition : public PwCondition<TDim,TNumNodes>
{
public:
    using GeometryType   = Geometry<Node>;
    using PropertiesType = Properties;
    using NodesArrayType = GeometryType::PointsArrayType;
    using BaseType       = PwCondition<TDim, TNumNodes>;
    
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(PwPointFluxCondition);

    PwPointFluxCondition();

    PwPointFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    PwPointFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~PwPointFluxCondition() override = default;

    Condition::Pointer Create(IndexType               NewId,
                              NodesArrayType const&   rThisNodes,
                              PropertiesType::Pointer pProperties) const override;
 
protected:
    void CalculateRHS(Vector& rRightHandSideVector, const ProcessInfo&) override;

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