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

// Project includes
#include "custom_conditions/axisymmetric_line_normal_fluid_flux_2D_diff_order_condition.h"
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{

AxisymmetricLineNormalFluidFlux2DDiffOrderCondition::AxisymmetricLineNormalFluidFlux2DDiffOrderCondition()
    : LineNormalFluidFlux2DDiffOrderCondition()
{
}

AxisymmetricLineNormalFluidFlux2DDiffOrderCondition::AxisymmetricLineNormalFluidFlux2DDiffOrderCondition(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : LineNormalFluidFlux2DDiffOrderCondition(NewId, std::move(pGeometry))
{
}

AxisymmetricLineNormalFluidFlux2DDiffOrderCondition::AxisymmetricLineNormalFluidFlux2DDiffOrderCondition(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : LineNormalFluidFlux2DDiffOrderCondition(NewId, std::move(pGeometry), std::move(pProperties))
{
}

Condition::Pointer AxisymmetricLineNormalFluidFlux2DDiffOrderCondition::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Create(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer AxisymmetricLineNormalFluidFlux2DDiffOrderCondition::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return make_intrusive<AxisymmetricLineNormalFluidFlux2DDiffOrderCondition>(NewId, pGeom, pProperties);
}

double AxisymmetricLineNormalFluidFlux2DDiffOrderCondition::CalculateIntegrationCoefficient(
    IndexType                                       PointNumber,
    const GeometryType::JacobiansType&              JContainer,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const

{
    KRATOS_TRY

    const double dx_dxi = JContainer[PointNumber](0, 0);
    const double dy_dxi = JContainer[PointNumber](1, 0);

    const double ds = sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);

    Vector N;
    N = this->GetGeometry().ShapeFunctionsValues(N, IntegrationPoints[PointNumber].Coordinates());
    const double radiusWeight =
        GeoElementUtilities::CalculateAxisymmetricCircumference(N, this->GetGeometry());

    return ds * IntegrationPoints[PointNumber].Weight() * radiusWeight;

    KRATOS_CATCH("")
}

std::string AxisymmetricLineNormalFluidFlux2DDiffOrderCondition::Info() const
{
    return "AxisymmetricLineNormalFluidFlux2DDiffOrderCondition";
}

void AxisymmetricLineNormalFluidFlux2DDiffOrderCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LineNormalFluidFlux2DDiffOrderCondition)
}

void AxisymmetricLineNormalFluidFlux2DDiffOrderCondition::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LineNormalFluidFlux2DDiffOrderCondition)
}

} // Namespace Kratos.
