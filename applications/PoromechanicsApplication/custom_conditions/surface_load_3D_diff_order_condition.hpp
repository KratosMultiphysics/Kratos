//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_SURFACE_LOAD_3D_DIFF_ORDER_CONDITION_H_INCLUDED )
#define  KRATOS_SURFACE_LOAD_3D_DIFF_ORDER_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_conditions/general_U_Pw_diff_order_condition.hpp"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) SurfaceLoad3DDiffOrderCondition : public GeneralUPwDiffOrderCondition
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( SurfaceLoad3DDiffOrderCondition );
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    SurfaceLoad3DDiffOrderCondition();
    
    // Constructor 1
    SurfaceLoad3DDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry );
    
    // Constructor 2
    SurfaceLoad3DDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    ~SurfaceLoad3DDiffOrderCondition() override;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, GeneralUPwDiffOrderCondition )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, GeneralUPwDiffOrderCondition )
    }

}; // class SurfaceLoad3DDiffOrderCondition.

} // namespace Kratos.

#endif // KRATOS_SURFACE_LOAD_3D_DIFF_ORDER_CONDITION_H_INCLUDED defined 
