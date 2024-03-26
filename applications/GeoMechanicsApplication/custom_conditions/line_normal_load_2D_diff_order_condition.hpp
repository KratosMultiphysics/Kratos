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
#include "includes/serializer.h"
#include "custom_conditions/line_load_2D_diff_order_condition.hpp"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION)
    LineNormalLoad2DDiffOrderCondition : public LineLoad2DDiffOrderCondition
{

public:

    using IndexType = std::size_t;
    using PropertiesType = Properties;
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType = Vector;
    using MatrixType = Matrix;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( LineNormalLoad2DDiffOrderCondition );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    LineNormalLoad2DDiffOrderCondition();

    // Constructor 1
    LineNormalLoad2DDiffOrderCondition( IndexType NewId,
                                        GeometryType::Pointer pGeometry );

    // Constructor 2
    LineNormalLoad2DDiffOrderCondition( IndexType NewId,
                                        GeometryType::Pointer pGeometry,
                                        PropertiesType::Pointer pProperties );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer
        Create(IndexType NewId,
               NodesArrayType const& ThisNodes,
               PropertiesType::Pointer pProperties ) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateConditionVector(ConditionVariables& rVariables,
                                  unsigned int PointNumber) override;
    double CalculateIntegrationCoefficient(const IndexType PointNumber,
                                           const GeometryType::JacobiansType& JContainer,
                                           const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const override;

    void CalculateAndAddConditionForce(VectorType& rRightHandSideVector,
                                       ConditionVariables& rVariables) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LineLoad2DDiffOrderCondition )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LineLoad2DDiffOrderCondition )
    }

}; // class LineNormalLoad2DDiffOrderCondition.

} // namespace Kratos.