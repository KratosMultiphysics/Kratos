//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_LINE_NORMAL_FLUID_FLUX_2D_FIC_CONDITION_H_INCLUDED )
#define  KRATOS_LINE_NORMAL_FLUID_FLUX_2D_FIC_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_conditions/line_normal_fluid_flux_2D_condition.hpp"

namespace Kratos
{

class LineNormalFluidFlux2DFICCondition : public LineNormalFluidFlux2DCondition
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( LineNormalFluidFlux2DFICCondition );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    LineNormalFluidFlux2DFICCondition();
    
    // Constructor 1
    LineNormalFluidFlux2DFICCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor 2
    LineNormalFluidFlux2DFICCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~LineNormalFluidFlux2DFICCondition();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables);
    
    void CalculateAndAddBoundaryMassMatrix(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables);
    

    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ConditionVariables& rVariables);
    
    void CalculateAndAddBoundaryMassFlow(VectorType& rRightHandSideVector, ConditionVariables& rVariables);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;
    
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LineNormalFluidFlux2DCondition )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LineNormalFluidFlux2DCondition )
    }
    
}; // class LineNormalFluidFlux2DFICCondition.

} // namespace Kratos.

#endif // KRATOS_LINE_NORMAL_FLUID_FLUX_2D_FIC_CONDITION_H_INCLUDED defined 
