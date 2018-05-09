//
//   Project Name:        			KratosDamApplication $
//   Last Modified by:    $Author:    	  Lorenzo Gracia $
//   Date:                $Date:           	January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_FREE_SURFACE_CONDITION_H_INCLUDED )
#define  KRATOS_FREE_SURFACE_CONDITION_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/process_info.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class FreeSurfaceCondition : public Condition
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( FreeSurfaceCondition );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    FreeSurfaceCondition() : Condition() {}

    // Constructor 1
    FreeSurfaceCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : Condition(NewId, pGeometry) {}

    // Constructor 2
    FreeSurfaceCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : Condition(NewId, pGeometry, pProperties)
    {
        mThisIntegrationMethod = this->GetGeometry().GetDefaultIntegrationMethod();
    }

    // Destructor
    virtual ~FreeSurfaceCondition() {}

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

}; // class FreeSurfaceCondition.

} // namespace Kratos.

#endif // KRATOS_FREE_SURFACE_CONDITION_H_INCLUDED defined
