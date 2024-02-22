//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana and Danilo Cavalcanti
//


#if !defined(KRATOS_U_PL_PG_LIQUID_DISCHARGE_CONDITION_H_INCLUDED )
#define  KRATOS_U_PL_PG_LIQUID_DISCHARGE_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/two-phase_flow/U_Pl_Pg_condition.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPlPgLiquidDischargeCondition : public UPlPgCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPlPgLiquidDischargeCondition );
    
    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPlPgCondition<TDim,TNumNodes>::mThisIntegrationMethod;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPlPgLiquidDischargeCondition() : UPlPgCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    UPlPgLiquidDischargeCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPlPgCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    UPlPgLiquidDischargeCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPlPgCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

    // Destructor
    ~UPlPgLiquidDischargeCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;
  
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:   
    
    // Member Variables
        
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;
        
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Serialization
    
    friend class Serializer;
    
    void save(Serializer& rSerializer) const override
    {
        typedef UPlPgCondition<TDim,TNumNodes> BaseCondition;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseCondition )
    }

    void load(Serializer& rSerializer) override
    {
        typedef UPlPgCondition<TDim,TNumNodes> BaseCondition;
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseCondition )
    }
    
}; // class UPlPgLiquidDischargeCondition.

} // namespace Kratos.

#endif // KRATOS_U_PL_PG_LIQUID_DISCHARGE_CONDITION_H_INCLUDED defined 
