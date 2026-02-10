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

// Application includes
#include "custom_conditions/U_Pw_normal_flux_condition.h"
#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/variables_utilities.hpp"

#include <numeric>

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
UPwNormalFluxCondition<TDim, TNumNodes>::UPwNormalFluxCondition()
    : UPwFaceLoadCondition<TDim, TNumNodes>()
{
}

template <unsigned int TDim, unsigned int TNumNodes>
UPwNormalFluxCondition<TDim, TNumNodes>::UPwNormalFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : UPwFaceLoadCondition<TDim, TNumNodes>(NewId, pGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
UPwNormalFluxCondition<TDim, TNumNodes>::UPwNormalFluxCondition(IndexType               NewId,
                                                                GeometryType::Pointer   pGeometry,
                                                                PropertiesType::Pointer pProperties)
    : UPwFaceLoadCondition<TDim, TNumNodes>(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer UPwNormalFluxCondition<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                   NodesArrayType const& ThisNodes,
                                                                   PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPwNormalFluxCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFluxCondition<TDim, TNumNodes>::CalculateRHS(Vector& rRightHandSideVector, const ProcessInfo&)
{
    const auto& r_geometry           = this->GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int number_of_integration_points = r_integration_points.size();

    // Containers of variables at all integration points
    const Matrix& r_n_container = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    GeometryType::JacobiansType j_container(number_of_integration_points);
    for (auto& j : j_container)
        j.resize(TDim, r_geometry.LocalSpaceDimension(), false);
    r_geometry.Jacobian(j_container, this->GetIntegrationMethod());

    // Condition variables
    const auto normal_flux_vector = VariablesUtilities::GetNodalValues(r_geometry, NORMAL_FLUID_FLUX);

    for (unsigned int integration_point = 0; integration_point < number_of_integration_points; ++integration_point) {
        const auto shape_function_values = row(r_n_container, integration_point);
        const auto normal_flux           = std::inner_product(
            shape_function_values.begin(), shape_function_values.end(), normal_flux_vector.cbegin(), 0.0);

        // Compute weighting coefficient for integration
        auto integration_coefficient = ConditionUtilities::CalculateIntegrationCoefficient(
            j_container[integration_point], r_integration_points[integration_point].Weight());

        // Contributions to the right hand side
        GeoElementUtilities::AssemblePBlockVector(
            rRightHandSideVector, -normal_flux * shape_function_values * integration_coefficient);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string UPwNormalFluxCondition<TDim, TNumNodes>::Info() const
{
    return "UPwNormalFluxCondition";
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFluxCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFluxCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
}
template class UPwNormalFluxCondition<2, 2>;
template class UPwNormalFluxCondition<2, 3>;
template class UPwNormalFluxCondition<2, 4>;
template class UPwNormalFluxCondition<2, 5>;
template class UPwNormalFluxCondition<3, 3>;
template class UPwNormalFluxCondition<3, 4>;

} // Namespace Kratos.
