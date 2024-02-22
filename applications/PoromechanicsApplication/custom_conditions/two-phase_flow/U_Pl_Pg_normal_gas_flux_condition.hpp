//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Danilo Cavalcanti
//


#if !defined(KRATOS_U_PL_PG_NORMAL_GAS_FLUX_CONDITION_H_INCLUDED )
#define  KRATOS_U_PL_PG_NORMAL_GAS_FLUX_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/two-phase_flow/U_Pl_Pg_condition.hpp"
#include "custom_conditions/two-phase_flow/U_Pl_Pg_normal_liquid_flux_condition.hpp"
#include "custom_utilities/poro_element_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPlPgNormalGasFluxCondition : public UPlPgNormalLiquidFluxCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPlPgNormalGasFluxCondition );
    
    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPlPgCondition<TDim,TNumNodes>::mThisIntegrationMethod;
    typedef typename UPlPgNormalLiquidFluxCondition<TDim,TNumNodes>::NormalFluidFluxVariables NormalFluidFluxVariables;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPlPgNormalGasFluxCondition() : UPlPgNormalLiquidFluxCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    UPlPgNormalGasFluxCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPlPgNormalLiquidFluxCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    UPlPgNormalGasFluxCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPlPgNormalLiquidFluxCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

    // Destructor
    ~UPlPgNormalGasFluxCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;
 
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:   
        
    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetNormalFluidFluxVector(array_1d<double,TNumNodes>& rNormalFluidFluxVector, const GeometryType& Geom) override;

    void CalculateAndAddRHS(VectorType& rRightHandSideVector, NormalFluidFluxVariables& rVariables) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Serialization
    
    friend class Serializer;
    
    void save(Serializer& rSerializer) const override
    {
        typedef UPlPgNormalLiquidFluxCondition<TDim,TNumNodes> BaseCondition;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseCondition )
    }

    void load(Serializer& rSerializer) override
    {
        typedef UPlPgNormalLiquidFluxCondition<TDim,TNumNodes> BaseCondition;
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseCondition )
    }
    
}; // class UPlPgNormalGasFluxCondition.

} // namespace Kratos.

#endif // KRATOS_U_PL_PG_NORMAL_GAS_FLUX_CONDITION_H_INCLUDED defined 
