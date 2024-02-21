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


#if !defined(KRATOS_U_PL_PG_CONDITION_H_INCLUDED )
#define  KRATOS_U_PL_PG_CONDITION_H_INCLUDED

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_utilities/poro_condition_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPlPgCondition : public Condition
{

public:

    /// We define the base class Condition
    typedef Condition BaseType;

    /// Definition of the size type
    typedef BaseType::SizeType SizeType;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPlPgCondition );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPlPgCondition() : Condition() {}

    // Constructor 1
    UPlPgCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : Condition(NewId, pGeometry) {}

    // Constructor 2
    UPlPgCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : Condition(NewId, pGeometry, pProperties)
    {
        mThisIntegrationMethod = this->GetIntegrationMethod();
    }

    // Destructor
    virtual ~UPlPgCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

    void GetDofList(DofsVectorType& rConditionDofList,const ProcessInfo& rCurrentProcessInfo ) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,const ProcessInfo& rCurrentProcessInfo ) override;

    void EquationIdVector(EquationIdVectorType& rResult,const ProcessInfo& rCurrentProcessInfo ) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    GeometryData::IntegrationMethod mThisIntegrationMethod;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo );

    virtual void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

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

}; // class UPlPgCondition.

} // namespace Kratos.

#endif // KRATOS_U_PL_PG_CONDITION_H_INCLUDED defined
