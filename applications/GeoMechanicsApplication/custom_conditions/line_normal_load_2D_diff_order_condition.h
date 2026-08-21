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
#include "custom_conditions/line_load_2D_diff_order_condition.h"
#include "includes/serializer.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) LineNormalLoad2DDiffOrderCondition : public LineLoad2DDiffOrderCondition
{
public:
    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using GeometryType   = Geometry<Node>;
    using NodesArrayType = GeometryType::PointsArrayType;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LineNormalLoad2DDiffOrderCondition);

    LineNormalLoad2DDiffOrderCondition();

    LineNormalLoad2DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    LineNormalLoad2DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    Condition::Pointer Create(IndexType               NewId,
                              NodesArrayType const&   ThisNodes,
                              PropertiesType::Pointer pProperties) const override;
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    std::string Info() const override;

protected:
    void CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber) override;
    double CalculateIntegrationCoefficient(IndexType                          PointNumber,
                                           const GeometryType::JacobiansType& JContainer,
                                           const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const override;

    void CalculateAndAddConditionForce(Vector& rRightHandSideVector, ConditionVariables& rVariables) override;

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

}; // class LineNormalLoad2DDiffOrderCondition.

} // namespace Kratos.