//
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_LINE_NORMAL_LIQUID_FLUX_2D_DIFF_ORDER_CONDITION_H_INCLUDED )
#define  KRATOS_LINE_NORMAL_LIQUID_FLUX_2D_DIFF_ORDER_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_conditions/one-phase_flow/general_U_Pl_diff_order_condition.hpp"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) LineNormalLiquidFlux2DDiffOrderCondition : public GeneralUPlDiffOrderCondition
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( LineNormalLiquidFlux2DDiffOrderCondition );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
<<<<<<< HEAD:applications/PoromechanicsApplication/custom_conditions/line_normal_fluid_flux_2D_diff_order_condition.hpp
    LineNormalFluidFlux2DDiffOrderCondition();

=======
    LineNormalLiquidFlux2DDiffOrderCondition();
    
>>>>>>> master:applications/PoromechanicsApplication/custom_conditions/one-phase_flow/line_normal_liquid_flux_2D_diff_order_condition.hpp
    // Constructor 1
    LineNormalLiquidFlux2DDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor 2
    LineNormalLiquidFlux2DDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    ~LineNormalLiquidFlux2DDiffOrderCondition() override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber) override;

    void CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, GeneralUPlDiffOrderCondition )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, GeneralUPlDiffOrderCondition )
    }
<<<<<<< HEAD:applications/PoromechanicsApplication/custom_conditions/line_normal_fluid_flux_2D_diff_order_condition.hpp

}; // class LineNormalFluidFlux2DDiffOrderCondition.

} // namespace Kratos.

#endif // KRATOS_LINE_NORMAL_FLUID_FLUX_2D_DIFF_ORDER_CONDITION_H_INCLUDED defined
=======
    
}; // class LineNormalLiquidFlux2DDiffOrderCondition.

} // namespace Kratos.

#endif // KRATOS_LINE_NORMAL_LIQUID_FLUX_2D_DIFF_ORDER_CONDITION_H_INCLUDED defined 
>>>>>>> master:applications/PoromechanicsApplication/custom_conditions/one-phase_flow/line_normal_liquid_flux_2D_diff_order_condition.hpp
