//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_SURFACE_NORMAL_FLUID_FLUX_3D_FIC_CONDITION_H_INCLUDED )
#define  KRATOS_SURFACE_NORMAL_FLUID_FLUX_3D_FIC_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_conditions/surface_normal_fluid_flux_3D_condition.hpp"

namespace Kratos
{

class SurfaceNormalFluidFlux3DFICCondition : public SurfaceNormalFluidFlux3DCondition
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( SurfaceNormalFluidFlux3DFICCondition );
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    SurfaceNormalFluidFlux3DFICCondition();
    
    // Constructor 1
    SurfaceNormalFluidFlux3DFICCondition( IndexType NewId, GeometryType::Pointer pGeometry );
    
    // Constructor 2
    SurfaceNormalFluidFlux3DFICCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~SurfaceNormalFluidFlux3DFICCondition();

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SurfaceNormalFluidFlux3DCondition )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SurfaceNormalFluidFlux3DCondition )
    }

}; // class SurfaceNormalFluidFlux3DFICCondition.

} // namespace Kratos.

#endif // KRATOS_SURFACE_NORMAL_FLUID_FLUX_3D_FIC_CONDITION_H_INCLUDED defined 
