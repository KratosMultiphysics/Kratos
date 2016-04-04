//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_SURFACE_LOAD_3D_CONDITION_H_INCLUDED )
#define  KRATOS_SURFACE_LOAD_3D_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_conditions/general_U_Pw_condition.hpp"

#include "poromechanics_application_variables.h"


namespace Kratos
{

class SurfaceLoad3DCondition : public GeneralUPwCondition
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( SurfaceLoad3DCondition );
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    SurfaceLoad3DCondition();
    
    // Constructor 1
    SurfaceLoad3DCondition( IndexType NewId, GeometryType::Pointer pGeometry );
    
    // Constructor 2
    SurfaceLoad3DCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~SurfaceLoad3DCondition();

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

}; // class SurfaceLoad3DCondition.

} // namespace Kratos.

#endif // KRATOS_SURFACE_LOAD_3D_CONDITION_H_INCLUDED defined 
