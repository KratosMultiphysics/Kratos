//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_LINE_NORMAL_FLUID_FLUX_2D_DIFF_ORDER_CONDITION_H_INCLUDED )
#define  KRATOS_LINE_NORMAL_FLUID_FLUX_2D_DIFF_ORDER_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_conditions/general_U_Pw_diff_order_condition.hpp"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) LineNormalFluidFlux2DDiffOrderCondition : public GeneralUPwDiffOrderCondition
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( LineNormalFluidFlux2DDiffOrderCondition );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    LineNormalFluidFlux2DDiffOrderCondition();
    
    // Constructor 1
    LineNormalFluidFlux2DDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor 2
    LineNormalFluidFlux2DDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~LineNormalFluidFlux2DDiffOrderCondition();

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, GeneralUPwDiffOrderCondition )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, GeneralUPwDiffOrderCondition )
    }
    
}; // class LineNormalFluidFlux2DDiffOrderCondition.

} // namespace Kratos.

#endif // KRATOS_LINE_NORMAL_FLUID_FLUX_2D_DIFF_ORDER_CONDITION_H_INCLUDED defined 
