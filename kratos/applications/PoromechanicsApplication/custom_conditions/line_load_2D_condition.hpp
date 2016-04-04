//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_LINE_LOAD_2D_CONDITION_H_INCLUDED )
#define  KRATOS_LINE_LOAD_2D_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/general_U_Pw_condition.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class LineLoad2DCondition : public GeneralUPwCondition
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( LineLoad2DCondition );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    LineLoad2DCondition();
    
    // Constructor 1
    LineLoad2DCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor 2
    LineLoad2DCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~LineLoad2DCondition();

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, GeneralUPwCondition )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, GeneralUPwCondition )
    }
    
}; // class LineLoad2DCondition.

} // namespace Kratos.

#endif // KRATOS_LINE_LOAD_2D_CONDITION_H_INCLUDED defined 
