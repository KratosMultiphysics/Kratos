//
//   Project Name:         $
//   Last modified by:    $Author:             $
//   Date:                $Date:                 $
//   Revision:            $Revision:         $
//

#if !defined(KRATOS_LINE_LOAD_CONDITION_H_INCLUDED )
#define  KRATOS_LINE_LOAD_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_conditions/general_condition.hpp"

namespace Kratos
{

class LineLoadCondition : public GeneralCondition
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( LineLoadCondition );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    LineLoadCondition();
    
    // Constructor 1
    LineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor 2
    LineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~LineLoadCondition();

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
    
}; // class LineLoadCondition.

} // namespace Kratos.

#endif // KRATOS_LINE_LOAD_CONDITION_H_INCLUDED defined 
