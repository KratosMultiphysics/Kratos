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
//  Main authors:    Vahid Galavi
//

#pragma once

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/U_Pw_normal_face_load_condition.hpp"
#include "custom_utilities/condition_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION)
    AxisymmetricUPwNormalFaceLoadCondition : public UPwNormalFaceLoadCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( AxisymmetricUPwNormalFaceLoadCondition );

    using IndexType = std::size_t;
    using PropertiesType = Properties;
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType = Vector;
    using MatrixType = Matrix;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    AxisymmetricUPwNormalFaceLoadCondition() : AxisymmetricUPwNormalFaceLoadCondition(0, nullptr, nullptr) {}

    AxisymmetricUPwNormalFaceLoadCondition( IndexType               NewId,
                                            GeometryType::Pointer   pGeometry )
        : AxisymmetricUPwNormalFaceLoadCondition(NewId, pGeometry, nullptr)
    {}

    AxisymmetricUPwNormalFaceLoadCondition( IndexType               NewId,
                                            GeometryType::Pointer   pGeometry,
                                            PropertiesType::Pointer pProperties )
        : UPwNormalFaceLoadCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties)
    {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double CalculateIntegrationCoefficient( const IndexType PointNumber,
                                            const GeometryType::IntegrationPointsArrayType& IntegrationPoints ) const override;

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

}; // class AxisymmetricUPwNormalFaceLoadCondition.

} // namespace Kratos.
