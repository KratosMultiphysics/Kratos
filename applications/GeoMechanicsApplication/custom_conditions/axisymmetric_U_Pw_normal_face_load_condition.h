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
#include "custom_conditions/U_Pw_normal_face_load_condition.h"
#include "custom_utilities/condition_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) AxisymmetricUPwNormalFaceLoadCondition
    : public UPwNormalFaceLoadCondition<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AxisymmetricUPwNormalFaceLoadCondition);

    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using GeometryType   = Geometry<Node>;
    using NodesArrayType = GeometryType::PointsArrayType;

    AxisymmetricUPwNormalFaceLoadCondition();

    AxisymmetricUPwNormalFaceLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    AxisymmetricUPwNormalFaceLoadCondition(IndexType               NewId,
                                           GeometryType::Pointer   pGeometry,
                                           PropertiesType::Pointer pProperties);

    Condition::Pointer Create(IndexType               NewId,
                              NodesArrayType const&   ThisNodes,
                              PropertiesType::Pointer pProperties) const override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    std::string Info() const override;

protected:
    double CalculateIntegrationCoefficient(const IndexType PointNumber,
                                           const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const override;

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

}; // class AxisymmetricUPwNormalFaceLoadCondition.

} // namespace Kratos.
