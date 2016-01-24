//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_LINE_NORMAL_LOAD_2D_CONDITION_H_INCLUDED )
#define  KRATOS_LINE_NORMAL_LOAD_2D_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_conditions/general_U_Pw_condition.hpp"

namespace Kratos
{

class LineNormalLoad2DCondition : public GeneralUPwCondition
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( LineNormalLoad2DCondition );
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    LineNormalLoad2DCondition();
    
    // Constructor 1
    LineNormalLoad2DCondition( IndexType NewId, GeometryType::Pointer pGeometry );
    
    // Constructor 2
    LineNormalLoad2DCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~LineNormalLoad2DCondition();

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

}; // class LineNormalLoad2DCondition.

} // namespace Kratos.

#endif // KRATOS_LINE_NORMAL_LOAD_2D_CONDITION_H_INCLUDED defined 
