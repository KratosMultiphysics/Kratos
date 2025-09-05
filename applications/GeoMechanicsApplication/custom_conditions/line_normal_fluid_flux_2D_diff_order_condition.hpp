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
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

// Project includes
#include "custom_conditions/line_load_2D_diff_order_condition.hpp"
#include "includes/serializer.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) LineNormalFluidFlux2DDiffOrderCondition : public LineLoad2DDiffOrderCondition
{
public:
    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using GeometryType   = Geometry<Node>;
    using NodesArrayType = GeometryType::PointsArrayType;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LineNormalFluidFlux2DDiffOrderCondition);

    // Default constructor
    LineNormalFluidFlux2DDiffOrderCondition();

    // Constructor 1
    LineNormalFluidFlux2DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    // Constructor 2
    LineNormalFluidFlux2DDiffOrderCondition(IndexType               NewId,
                                            GeometryType::Pointer   pGeometry,
                                            PropertiesType::Pointer pProperties);

    Condition::Pointer Create(IndexType               NewId,
                              NodesArrayType const&   ThisNodes,
                              PropertiesType::Pointer pProperties) const override;
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    std::string Info() const override;

protected:
    // Member Variables

    void CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber) override;

    void CalculateAndAddConditionForce(Vector& rRightHandSideVector, ConditionVariables& rVariables) override;

private:
    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LineLoad2DDiffOrderCondition)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LineLoad2DDiffOrderCondition)
    }

}; // class LineNormalFluidFlux2DDiffOrderCondition.

} // namespace Kratos.