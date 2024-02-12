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
//


#if !defined(KRATOS_U_PL_PG_NORMAL_LIQUID_FLUX_INTERFACE_CONDITION_H_INCLUDED )
#define  KRATOS_U_PL_PG_NORMAL_LIQUID_FLUX_INTERFACE_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/multiphase_flow/U_Pl_Pg_condition.hpp"
#include "custom_conditions/multiphase_flow//U_Pl_Pg_face_load_interface_condition.hpp"
#include "custom_utilities/poro_element_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPlPgNormalLiquidFluxInterfaceCondition : public UPlPgFaceLoadInterfaceCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( UPlPgNormalLiquidFluxInterfaceCondition );
    
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
    UPlPgNormalLiquidFluxInterfaceCondition() : UPlPgFaceLoadInterfaceCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    UPlPgNormalLiquidFluxInterfaceCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPlPgFaceLoadInterfaceCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    UPlPgNormalLiquidFluxInterfaceCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPlPgFaceLoadInterfaceCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

    // Destructor
    ~UPlPgNormalLiquidFluxInterfaceCondition() override {}

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
    }
    
}; // class UPlPgNormalLiquidFluxInterfaceCondition.

} // namespace Kratos.

#endif // KRATOS_U_PL_PG_NORMAL_LIQUID_FLUX_INTERFACE_CONDITION_H_INCLUDED defined 
