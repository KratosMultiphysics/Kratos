//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_U_PW_NORMAL_FLUX_INTERFACE_CONDITION_H_INCLUDED )
#define  KRATOS_U_PW_NORMAL_FLUX_INTERFACE_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/U_Pw_condition.hpp"
#include "custom_conditions/U_Pw_face_load_interface_condition.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/interface_element_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPwNormalFluxInterfaceCondition : public UPwFaceLoadInterfaceCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( UPwNormalFluxInterfaceCondition );
    
    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPwCondition<TDim,TNumNodes>::mThisIntegrationMethod;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPwNormalFluxInterfaceCondition() : UPwFaceLoadInterfaceCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    UPwNormalFluxInterfaceCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPwFaceLoadInterfaceCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    UPwNormalFluxInterfaceCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPwFaceLoadInterfaceCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

    // Destructor
    virtual ~UPwNormalFluxInterfaceCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const;
 
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:   
    
    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                    
    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo );
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Serialization
    
    friend class Serializer;
    
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
    }
    
}; // class UPwNormalFluxInterfaceCondition.

} // namespace Kratos.

#endif // KRATOS_U_PW_NORMAL_FLUX_INTERFACE_CONDITION_H_INCLUDED defined 
