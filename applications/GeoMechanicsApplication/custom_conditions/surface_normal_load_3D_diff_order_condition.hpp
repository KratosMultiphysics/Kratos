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
#include "custom_conditions/surface_load_3D_diff_order_condition.hpp"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) SurfaceNormalLoad3DDiffOrderCondition : public SurfaceLoad3DDiffOrderCondition
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( SurfaceNormalLoad3DDiffOrderCondition );
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    SurfaceNormalLoad3DDiffOrderCondition();
    
    // Constructor 1
    SurfaceNormalLoad3DDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry );
    
    // Constructor 2
    SurfaceNormalLoad3DDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber) override;

    double CalculateIntegrationCoefficient(const IndexType PointNumber,
                                           const GeometryType::JacobiansType& JContainer,
                                           const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const override;

    void CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;
    
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SurfaceLoad3DDiffOrderCondition )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SurfaceLoad3DDiffOrderCondition )
    }
    
}; // class SurfaceNormalLoad3DDiffOrderCondition.

} // namespace Kratos.
