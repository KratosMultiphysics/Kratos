//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//


#if !defined(KRATOS_ADDED_MASS_CONDITION_H_INCLUDED )
#define  KRATOS_ADDED_MASS_CONDITION_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/process_info.h"

// Application includes
#include "custom_utilities/poro_condition_utilities.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class AddedMassCondition : public Condition
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( AddedMassCondition );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    AddedMassCondition() : Condition() {}

    // Constructor 1
    AddedMassCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : Condition(NewId, pGeometry) {}

    // Constructor 2
    AddedMassCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : Condition(NewId, pGeometry, pProperties)
    {
        mThisIntegrationMethod = this->GetGeometry().GetDefaultIntegrationMethod();
    }

    // Destructor
    virtual ~AddedMassCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

    void GetDofList(DofsVectorType& rConditionDofList,ProcessInfo& rCurrentProcessInfo ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo ) override;

    void EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    GeometryData::IntegrationMethod mThisIntegrationMethod;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo );

    virtual void CalculateLHS( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo );

    virtual void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo );

    void CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& weight);

    void GetAccelerationVector( Vector& rValues, int Step );

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

}; // class AddedMassCondition.

} // namespace Kratos.

#endif // KRATOS_ADDED_MASS_CONDITION_H_INCLUDED defined
