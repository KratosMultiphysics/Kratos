#include "custom_utilities/boussinesq_coupling_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include <cmath>
#include <algorithm>

namespace Kratos
{

double BoussinesqCouplingUtilities::ComputeRelativeResidual(
    ModelPart& rModelPart,
    const Variable<double>& rVarNew,
    const Variable<double>& rVarOld)
{
    const auto& r_elements = rModelPart.Elements();

    double diff_norm_sq = block_for_each<SumReduction<double>>(
        r_elements,
        [&](const Element& rElement) -> double {
            double local_diff_norm_sq = 0.0;
            const auto& r_geometry = rElement.GetGeometry();
            const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();
            const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
            const unsigned int number_of_integration_points = r_integration_points.size();
            const unsigned int number_of_nodes = r_geometry.size();

            Vector detJ0(number_of_integration_points);
            r_geometry.DeterminantOfJacobian(detJ0, r_integration_method);
            const Matrix& N_container = r_geometry.ShapeFunctionsValues(r_integration_method);

            for (unsigned int g = 0; g < number_of_integration_points; ++g) {
                double v_new_gauss = 0.0;
                double v_old_gauss = 0.0;
                for (unsigned int i = 0; i < number_of_nodes; ++i) {
                    const double N_i = N_container(g, i);
                    v_new_gauss += N_i * r_geometry[i].FastGetSolutionStepValue(rVarNew);
                    v_old_gauss += N_i * r_geometry[i].GetValue(rVarOld);
                }
                const double weight = r_integration_points[g].Weight() * detJ0[g];
                local_diff_norm_sq += std::pow(v_new_gauss - v_old_gauss, 2) * weight;
            }
            return local_diff_norm_sq;
        }
    );

    double old_norm_sq = block_for_each<SumReduction<double>>(
        r_elements,
        [&](const Element& rElement) -> double {
            double local_old_norm_sq = 0.0;
            const auto& r_geometry = rElement.GetGeometry();
            const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();
            const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
            const unsigned int number_of_integration_points = r_integration_points.size();
            const unsigned int number_of_nodes = r_geometry.size();

            Vector detJ0(number_of_integration_points);
            r_geometry.DeterminantOfJacobian(detJ0, r_integration_method);
            const Matrix& N_container = r_geometry.ShapeFunctionsValues(r_integration_method);

            for (unsigned int g = 0; g < number_of_integration_points; ++g) {
                double v_old_gauss = 0.0;
                for (unsigned int i = 0; i < number_of_nodes; ++i) {
                    const double N_i = N_container(g, i);
                    v_old_gauss += N_i * r_geometry[i].GetValue(rVarOld);
                }
                const double weight = r_integration_points[g].Weight() * detJ0[g];
                local_old_norm_sq += std::pow(v_old_gauss, 2) * weight;
            }
            return local_old_norm_sq;
        }
    );

    return std::sqrt(diff_norm_sq) / std::max(std::sqrt(old_norm_sq), 1.0e-10);
}

void BoussinesqCouplingUtilities::ApplyRelaxation(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const Variable<double>& rOldVariable,
    const double Omega)
{
    const double one_minus_omega = 1.0 - Omega;
    auto& r_nodes = rModelPart.Nodes();

    block_for_each(r_nodes, [&](Node& rNode) {
        const double v_new = rNode.FastGetSolutionStepValue(rVariable);
        const double v_old = rNode.GetValue(rOldVariable);
        rNode.FastGetSolutionStepValue(rVariable) = Omega * v_new + one_minus_omega * v_old;
    });
}

} // namespace Kratos
