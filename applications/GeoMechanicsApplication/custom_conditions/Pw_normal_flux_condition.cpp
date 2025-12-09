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
//                   Vahid Galavi,
//                   Aron Noordam
//

// Application includes
#include "custom_conditions/Pw_normal_flux_condition.h"
#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/variables_utilities.hpp"

#include <numeric>

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer PwNormalFluxCondition<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                  NodesArrayType const& ThisNodes,
                                                                  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new PwNormalFluxCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
void PwNormalFluxCondition<TDim, TNumNodes>::CalculateRHS(Vector&            rRightHandSideVector,
                                                          const ProcessInfo& CurrentProcessInfo)
{
    const auto& r_geometry           = this->GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int number_of_integration_points = r_integration_points.size();

    // Containers of variables at all integration points
    const Matrix& r_n_container = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    GeometryType::JacobiansType j_container(number_of_integration_points);
    for (auto& j : j_container) {
        j.resize(TDim, r_geometry.LocalSpaceDimension(), false);
    }
    r_geometry.Jacobian(j_container, this->GetIntegrationMethod());

    // Condition variables
    Vector normal_flux_vector(TNumNodes);
    VariablesUtilities::GetNodalValues(r_geometry, NORMAL_FLUID_FLUX, normal_flux_vector.begin());

    for (unsigned int integration_point = 0; integration_point < number_of_integration_points; ++integration_point) {
        // Interpolation of nodal normal flux to integration point normal flux.
        const auto shape_function_values = row(r_n_container, integration_point);
        const auto normal_flux           = std::inner_product(
            shape_function_values.begin(), shape_function_values.end(), normal_flux_vector.cbegin(), 0.0);

        // Compute weighting coefficient for integration
        auto integration_coefficient = ConditionUtilities::CalculateIntegrationCoefficient(
            j_container[integration_point], r_integration_points[integration_point].Weight());

        // Contributions to the right hand side
        rRightHandSideVector -= normal_flux * row(r_n_container, integration_point) * integration_coefficient;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void PwNormalFluxCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
}

template <unsigned int TDim, unsigned int TNumNodes>
void PwNormalFluxCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string PwNormalFluxCondition<TDim, TNumNodes>::Info() const
{
    return "PwNormalFluxCondition";
}

template class PwNormalFluxCondition<2, 2>;
template class PwNormalFluxCondition<2, 3>;
template class PwNormalFluxCondition<2, 4>;
template class PwNormalFluxCondition<2, 5>;
template class PwNormalFluxCondition<3, 3>;
template class PwNormalFluxCondition<3, 4>;

} // Namespace Kratos.
