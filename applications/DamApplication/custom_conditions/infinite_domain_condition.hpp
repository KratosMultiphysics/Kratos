//
//   Project Name:        			KratosDamApplication $
//   Last Modified by:    $Author:    	  Lorenzo Gracia $
//   Date:                $Date:           	January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_INFINITE_DOMAIN_CONDITION_H_INCLUDED )
#define  KRATOS_INFINITE_DOMAIN_CONDITION_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "custom_conditions/free_surface_condition.hpp"
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/process_info.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class InfiniteDomainCondition : public  FreeSurfaceCondition<TDim,TNumNodes>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( InfiniteDomainCondition );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using FreeSurfaceCondition<TDim,TNumNodes>::mThisIntegrationMethod;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    InfiniteDomainCondition() : FreeSurfaceCondition<TDim,TNumNodes>() {}

    // Constructor 1
    InfiniteDomainCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : FreeSurfaceCondition<TDim,TNumNodes>(NewId, pGeometry) {}

    // Constructor 2
    InfiniteDomainCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : FreeSurfaceCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties)
    {
    }

    // Destructor
    virtual ~InfiniteDomainCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


    void CalculateLHS( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

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

}; // class InfiniteDomainCondition.

} // namespace Kratos.

#endif // KRATOS_INFINITE_DOMAIN_CONDITION_H_INCLUDED defined
