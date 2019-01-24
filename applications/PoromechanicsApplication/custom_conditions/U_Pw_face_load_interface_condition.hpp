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


#if !defined(KRATOS_U_PW_FACE_LOAD_INTERFACE_CONDITION_H_INCLUDED )
#define  KRATOS_U_PW_FACE_LOAD_INTERFACE_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/U_Pw_condition.hpp"
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
    ~UPwFaceLoadInterfaceCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;
 
    void Initialize() override;
 
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:   
    
    // Member Variables
    
    std::vector<double> mInitialGap;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateInitialGap(const GeometryType& Geom);
    
    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo ) override;

    void CheckJointWidth(double& rJointWidth, bool& rComputeJointWidth, BoundedMatrix<double,TDim,TDim>& rRotationMatrix,
                            const double& MinimumJointWidth, const GeometryType& Geom);

    void CalculateJointWidth( double& rJointWidth, const BoundedMatrix<double,TDim,TDim*TNumNodes>& Nu,
                                const array_1d<double,TDim*TNumNodes>& DisplacementVector, array_1d<double,TDim>& rRelDispVector,
                                const BoundedMatrix<double,TDim,TDim>& RotationMatrix,
                                array_1d<double,TDim>& rLocalRelDispVector, const double& MinimumJointWidth, const unsigned int& GPoint );

    void CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& Weight, const double& JointWidth);
        
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
    
}; // class UPwFaceLoadInterfaceCondition.

} // namespace Kratos.

#endif // KRATOS_U_PW_FACE_LOAD_INTERFACE_CONDITION_H_INCLUDED defined 
