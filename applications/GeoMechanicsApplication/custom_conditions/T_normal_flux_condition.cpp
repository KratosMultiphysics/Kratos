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
//  Main authors:    Mohamed Nabi
//                   John van Esch
//

// Application includes
#include "custom_conditions/T_normal_flux_condition.h"
#include "custom_utilities/condition_utilities.hpp"
#include "utilities/math_utils.h"
#include <custom_utilities/variables_utilities.hpp>

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
GeoTNormalFluxCondition<TDim, TNumNodes>::GeoTNormalFluxCondition()
    : GeoTCondition<TDim, TNumNodes>()
{
}

template <unsigned int TDim, unsigned int TNumNodes>
GeoTNormalFluxCondition<TDim, TNumNodes>::GeoTNormalFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : GeoTCondition<TDim, TNumNodes>(NewId, pGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
GeoTNormalFluxCondition<TDim, TNumNodes>::GeoTNormalFluxCondition(IndexType             NewId,
                                                                  GeometryType::Pointer pGeometry,
                                                                  PropertiesType::Pointer pProperties)
    : GeoTCondition<TDim, TNumNodes>(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTNormalFluxCondition<TDim, TNumNodes>::CalculateRHS(Vector&            rRightHandSideVector,
                                                            const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType&                             r_geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points =
        r_geom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int num_integration_points = r_integration_points.size();

    const Matrix& r_N_container = r_geom.ShapeFunctionsValues(this->GetIntegrationMethod());

    GeometryType::JacobiansType j_container(num_integration_points);
    for (auto& j : j_container) {
        j.resize(TDim, r_geom.LocalSpaceDimension(), false);
    }
    r_geom.Jacobian(j_container, this->GetIntegrationMethod());

    Vector normal_flux_vector(TNumNodes);
    VariablesUtilities::GetNodalValues(r_geom, NORMAL_HEAT_FLUX, normal_flux_vector.begin());

    for (unsigned int integration_point = 0; integration_point < num_integration_points; ++integration_point) {
        // Obtain N
        auto N = row(r_N_container, integration_point);

        // Interpolation of nodal normal flux to integration point normal flux.
        auto normal_flux_on_integration_point = MathUtils<>::Dot(N, normal_flux_vector);

        auto weighting_integration_coefficient = ConditionUtilities::CalculateIntegrationCoefficient(
            j_container[integration_point], r_integration_points[integration_point].Weight());

        // Contributions to the right hand side
        rRightHandSideVector += normal_flux_on_integration_point * N * weighting_integration_coefficient;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string GeoTNormalFluxCondition<TDim, TNumNodes>::Info() const
{
    return "GeoTNormalFluxCondition";
}

template class GeoTNormalFluxCondition<2, 2>;
template class GeoTNormalFluxCondition<2, 3>;
template class GeoTNormalFluxCondition<2, 4>;
template class GeoTNormalFluxCondition<2, 5>;
template class GeoTNormalFluxCondition<3, 3>;
template class GeoTNormalFluxCondition<3, 4>;
template class GeoTNormalFluxCondition<3, 6>;
template class GeoTNormalFluxCondition<3, 8>;
template class GeoTNormalFluxCondition<3, 9>;

} // namespace Kratos
