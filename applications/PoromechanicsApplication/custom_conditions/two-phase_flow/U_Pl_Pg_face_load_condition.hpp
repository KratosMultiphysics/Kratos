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


#if !defined(KRATOS_U_PL_PG_FACE_LOAD_CONDITION_H_INCLUDED )
#define  KRATOS_U_PL_PG_FACE_LOAD_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/two-phase_flow/U_Pl_Pg_condition.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPlPgFaceLoadCondition : public UPlPgCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPlPgFaceLoadCondition );
    
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
    UPlPgFaceLoadCondition() : UPlPgCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    UPlPgFaceLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPlPgCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    UPlPgFaceLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPlPgCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

    // Destructor
    ~UPlPgFaceLoadCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;
  
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:   
    
    // Member Variables
        
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& Weight);
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Serialization
    
    friend class Serializer;
    
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, UPlPgCondition )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, UPlPgCondition )
    }
    
}; // class UPlPgFaceLoadCondition.

} // namespace Kratos.

#endif // KRATOS_U_PL_PG_FACE_LOAD_CONDITION_H_INCLUDED defined 
