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


#if !defined(KRATOS_U_PL_FACE_LOAD_INTERFACE_CONDITION_H_INCLUDED )
#define  KRATOS_U_PL_FACE_LOAD_INTERFACE_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/one-phase_flow/U_Pl_condition.hpp"
#include "custom_utilities/interface_element_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPlFaceLoadInterfaceCondition : public UPlCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPlFaceLoadInterfaceCondition );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPlCondition<TDim,TNumNodes>::mThisIntegrationMethod;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPlFaceLoadInterfaceCondition() : UPlCondition<TDim,TNumNodes>() {}

    // Constructor 1
    UPlFaceLoadInterfaceCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPlCondition<TDim,TNumNodes>(NewId, pGeometry) {}

    // Constructor 2
    UPlFaceLoadInterfaceCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPlCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties)
    {
        // Lobatto integration method with the integration points located at the "mid plane nodes" of the interface
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
    }

    // Destructor
    ~UPlFaceLoadInterfaceCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    std::vector<double> mInitialGap;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateInitialGap(const GeometryType& Geom);

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

    void CheckJointWidth(double& rJointWidth, bool& rComputeJointWidth, BoundedMatrix<double,TDim,TDim>& rRotationMatrix,
                            const double& InitialJointWidth, const GeometryType& Geom);

    void CalculateJointWidth( double& rJointWidth, const BoundedMatrix<double,TDim,TDim*TNumNodes>& Nu,
                                const array_1d<double,TDim*TNumNodes>& DisplacementVector, array_1d<double,TDim>& rRelDispVector,
                                const BoundedMatrix<double,TDim,TDim>& RotationMatrix,
                                array_1d<double,TDim>& rLocalRelDispVector, const double& InitialJointWidth, const unsigned int& GPoint );

    void CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& Weight, const double& JointWidth);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        typedef UPlCondition<TDim,TNumNodes> BaseCondition;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseCondition )
    }

    void load(Serializer& rSerializer) override
    {
        typedef UPlCondition<TDim,TNumNodes> BaseCondition;
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseCondition )
    }

}; // class UPlFaceLoadInterfaceCondition.

} // namespace Kratos.

#endif // KRATOS_U_PL_FACE_LOAD_INTERFACE_CONDITION_H_INCLUDED defined
