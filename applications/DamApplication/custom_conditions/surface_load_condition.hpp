//
//   Project Name:         $
//   Last modified by:    $Author:             $
//   Date:                $Date:               $
//   Revision:            $Revision:            $
//

#if !defined(KRATOS_SURFACE_LOAD_CONDITION_H_INCLUDED )
#define  KRATOS_SURFACE_LOAD_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_conditions/general_condition.hpp"

namespace Kratos
{

class SurfaceLoadCondition : public GeneralCondition
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( SurfaceLoadCondition );
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    SurfaceLoadCondition();
    
    // Constructor 1
    SurfaceLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry );
    
    // Constructor 2
    SurfaceLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~SurfaceLoadCondition();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber);

    void CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight);
    
    void CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;
    
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, GeneralCondition )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, GeneralCondition )
    }

}; // class SurfaceLoadCondition.

} // namespace Kratos.

#endif // KRATOS_SURFACE_LOAD_CONDITION_H_INCLUDED defined 
