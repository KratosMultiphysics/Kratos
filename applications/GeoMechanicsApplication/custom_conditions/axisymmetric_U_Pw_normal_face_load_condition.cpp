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

// Application includes
#include "custom_conditions/axisymmetric_U_Pw_normal_face_load_condition.h"
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
AxisymmetricUPwNormalFaceLoadCondition<TDim, TNumNodes>::AxisymmetricUPwNormalFaceLoadCondition()
    : AxisymmetricUPwNormalFaceLoadCondition(0, nullptr, nullptr)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
AxisymmetricUPwNormalFaceLoadCondition<TDim, TNumNodes>::AxisymmetricUPwNormalFaceLoadCondition(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : AxisymmetricUPwNormalFaceLoadCondition(NewId, pGeometry, nullptr)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
AxisymmetricUPwNormalFaceLoadCondition<TDim, TNumNodes>::AxisymmetricUPwNormalFaceLoadCondition(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : UPwNormalFaceLoadCondition<TDim, TNumNodes>(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer AxisymmetricUPwNormalFaceLoadCondition<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new AxisymmetricUPwNormalFaceLoadCondition(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
double AxisymmetricUPwNormalFaceLoadCondition<TDim, TNumNodes>::CalculateIntegrationCoefficient(
    const IndexType PointNumber, const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const
{
    Vector N;
    N = this->GetGeometry().ShapeFunctionsValues(N, IntegrationPoints[PointNumber].Coordinates());
    const double radiusWeight =
        GeoElementUtilities::CalculateAxisymmetricCircumference(N, this->GetGeometry());

    return IntegrationPoints[PointNumber].Weight() * radiusWeight;
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod AxisymmetricUPwNormalFaceLoadCondition<TDim, TNumNodes>::GetIntegrationMethod() const
{
    switch (this->GetGeometry().GetGeometryOrderType()) {
        using enum GeometryData::KratosGeometryOrderType;
        using enum GeometryData::IntegrationMethod;
    case Kratos_Cubic_Order:
        return GI_GAUSS_3;
    case Kratos_Quartic_Order:
        return GI_GAUSS_5;
    default:
        return GI_GAUSS_2;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string AxisymmetricUPwNormalFaceLoadCondition<TDim, TNumNodes>::Info() const
{
    return "AxisymmetricUPwNormalFaceLoadCondition";
}

template <unsigned int TDim, unsigned int TNumNodes>
void AxisymmetricUPwNormalFaceLoadCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
}

template <unsigned int TDim, unsigned int TNumNodes>
void AxisymmetricUPwNormalFaceLoadCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
}

template class AxisymmetricUPwNormalFaceLoadCondition<2, 2>;
template class AxisymmetricUPwNormalFaceLoadCondition<2, 3>;
template class AxisymmetricUPwNormalFaceLoadCondition<2, 4>;
template class AxisymmetricUPwNormalFaceLoadCondition<2, 5>;

} // Namespace Kratos.
