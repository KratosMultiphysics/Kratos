// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/U_Pw_condition.hpp"
#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/interface_element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwFaceLoadInterfaceCondition : public UPwCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPwFaceLoadInterfaceCondition );

    using IndexType = std::size_t;
    using PropertiesType = Properties;
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType = Vector;
    using MatrixType = Matrix;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    UPwFaceLoadInterfaceCondition() : UPwFaceLoadInterfaceCondition(0, nullptr, nullptr) {}

    UPwFaceLoadInterfaceCondition( IndexType               NewId,
                                   GeometryType::Pointer   pGeometry )
        : UPwFaceLoadInterfaceCondition(NewId, pGeometry, nullptr)
    {}

    UPwFaceLoadInterfaceCondition( IndexType               NewId,
                                   GeometryType::Pointer   pGeometry,
                                   PropertiesType::Pointer pProperties )
        : UPwCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties)
    {
        // Lobatto integration method with the integration points located at the "mid plane nodes" of the interface
        this->SetIntegrationMethod(GeometryData::IntegrationMethod::GI_GAUSS_1);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    void CalculateInitialGap(const GeometryType& Geom);

    void CalculateRHS(VectorType& rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo) override;

    void CheckJointWidth(double& rJointWidth, bool& rComputeJointWidth, BoundedMatrix<double,TDim,TDim>& rRotationMatrix,
                            const double& MinimumJointWidth, const GeometryType& Geom);

    void CalculateJointWidth( double& rJointWidth, const BoundedMatrix<double,TDim,TDim*TNumNodes>& Nu,
                                const array_1d<double,TDim*TNumNodes>& DisplacementVector, array_1d<double,TDim>& rRelDispVector,
                                const BoundedMatrix<double,TDim,TDim>& RotationMatrix,
                                array_1d<double,TDim>& rLocalRelDispVector, const double& MinimumJointWidth, const unsigned int& GPoint );

    double CalculateIntegrationCoefficient(const Matrix& Jacobian, const double& Weight, const double& JointWidth);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    Vector mInitialGap;

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
