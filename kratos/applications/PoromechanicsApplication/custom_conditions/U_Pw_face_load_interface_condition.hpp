//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_U_PW_FACE_LOAD_INTERFACE_CONDITION_H_INCLUDED )
#define  KRATOS_U_PW_FACE_LOAD_INTERFACE_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/U_Pw_condition.hpp"
#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/interface_element_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPwFaceLoadInterfaceCondition : public UPwCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( UPwFaceLoadInterfaceCondition );
    
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
    UPwFaceLoadInterfaceCondition() : UPwCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    UPwFaceLoadInterfaceCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPwCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    UPwFaceLoadInterfaceCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPwCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties)
    {
        // Lobatto integration method with the integration points located at the "mid plane nodes" of the interface
        mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
    }

    // Destructor
    virtual ~UPwFaceLoadInterfaceCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const;
 
    void Initialize();
 
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:   
    
    // Member Variables
    
    Vector mInitialGap;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateInitialGap(const GeometryType& Geom);
    
    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo );

    void CheckJointWidth(double& rJointWidth, bool& rComputeJointWidth, boost::numeric::ublas::bounded_matrix<double,TDim,TDim>& rRotationMatrix,
                            const double& MinimumJointWidth, const GeometryType& Geom);

    void CalculateJointWidth( double& rJointWidth, const boost::numeric::ublas::bounded_matrix<double,TDim,TDim*TNumNodes>& Nu,
                                const array_1d<double,TDim*TNumNodes>& DisplacementVector, array_1d<double,TDim>& rRelDispVector,
                                const boost::numeric::ublas::bounded_matrix<double,TDim,TDim>& RotationMatrix,
                                array_1d<double,TDim>& rLocalRelDispVector, const double& MinimumJointWidth, const unsigned int& GPoint );

    void CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& Weight, const double& JointWidth);
        
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
    
}; // class UPwFaceLoadInterfaceCondition.

} // namespace Kratos.

#endif // KRATOS_U_PW_FACE_LOAD_INTERFACE_CONDITION_H_INCLUDED defined 
