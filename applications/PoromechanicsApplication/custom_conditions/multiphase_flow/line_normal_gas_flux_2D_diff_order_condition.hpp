//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_LINE_NORMAL_GAS_FLUX_2D_DIFF_ORDER_CONDITION_H_INCLUDED )
#define  KRATOS_LINE_NORMAL_GAS_FLUX_2D_DIFF_ORDER_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_conditions/multiphase_flow/general_U_Pw_Pg_diff_order_condition.hpp"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) LineNormalGasFlux2DDiffOrderCondition : public GeneralUPwPgDiffOrderCondition
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( LineNormalGasFlux2DDiffOrderCondition );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    LineNormalGasFlux2DDiffOrderCondition();
    
    // Constructor 1
    LineNormalGasFlux2DDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor 2
    LineNormalGasFlux2DDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    ~LineNormalGasFlux2DDiffOrderCondition() override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber) override;

    void CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight) override;
    
    void CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;
    
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, GeneralUPwPgDiffOrderCondition )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, GeneralUPwPgDiffOrderCondition )
    }
    
}; // class LineNormalGasFlux2DDiffOrderCondition.

} // namespace Kratos.

#endif // KRATOS_LINE_NORMAL_GAS_FLUX_2D_DIFF_ORDER_CONDITION_H_INCLUDED defined 
