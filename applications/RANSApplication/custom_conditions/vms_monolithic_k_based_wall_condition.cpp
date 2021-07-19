//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <string>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "vms_monolithic_k_based_wall_condition.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
int VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::Check(
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    int check = BaseType::Check(rCurrentProcessInfo);

    const auto& r_geometry = this->GetGeometry();

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const auto& r_node = r_geometry[i_node];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
    }

    return check;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (RansCalculationUtilities::IsWallFunctionActive(*this)) {
        const array_1d<double, 3>& r_normal = this->GetValue(NORMAL);
        KRATOS_ERROR_IF(norm_2(r_normal) == 0.0)
            << "NORMAL must be calculated before using this " << this->Info() << "\n";

        KRATOS_ERROR_IF(this->GetValue(NEIGHBOUR_ELEMENTS).size() == 0)
            << this->Info() << " cannot find parent element\n";

        mWallHeight = RansCalculationUtilities::CalculateWallHeight(*this, r_normal);

        KRATOS_ERROR_IF(mWallHeight == 0.0) << this->Info() << " has zero wall height.\n";
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "VMSMonolithicKBasedWallCondition" << TDim << "D";
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "VMSMonolithicKBasedWallCondition";
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::ApplyWallLaw(
    MatrixType& rLocalMatrix,
    VectorType& rLocalVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    using namespace RansCalculationUtilities;

    if (IsWallFunctionActive(*this)) {
        const auto& r_geometry = this->GetGeometry();
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        CalculateConditionGeometryData(r_geometry, this->GetIntegrationMethod(),
                                       gauss_weights, shape_functions);
        const IndexType num_gauss_points = gauss_weights.size();

        const size_t block_size = TDim + 1;

        const double eps = std::numeric_limits<double>::epsilon();
        const double c_mu_25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
        const double kappa = rCurrentProcessInfo[VON_KARMAN];

        // get parent element
        auto& r_parent_element = this->GetValue(NEIGHBOUR_ELEMENTS)[0];
        auto p_constitutive_law = r_parent_element.GetValue(CONSTITUTIVE_LAW);

        // get fluid properties from parent element
        const auto& r_elem_properties = r_parent_element.GetProperties();
        const double rho = r_elem_properties[DENSITY];
        ConstitutiveLaw::Parameters cl_parameters(r_geometry, r_elem_properties, rCurrentProcessInfo);

        // get surface properties from condition
        const PropertiesType& r_cond_properties = this->GetProperties();
        const double beta = r_cond_properties.GetValue(WALL_SMOOTHNESS_BETA);
        const double y_plus_limit = r_cond_properties.GetValue(RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT);
        const double inv_kappa = 1.0 / kappa;

        double condition_y_plus{0.0};
        array_1d<double, 3> condition_u_tau = ZeroVector(3);

        double tke, nu;
        array_1d<double, 3> wall_velocity;

        for (size_t g = 0; g < num_gauss_points; ++g) {
            const Vector& gauss_shape_functions = row(shape_functions, g);

            cl_parameters.SetShapeFunctionsValues(gauss_shape_functions);
            p_constitutive_law->CalculateValue(cl_parameters, EFFECTIVE_VISCOSITY, nu);
            nu /= rho;

            FluidCalculationUtilities::EvaluateInPoint(
                r_geometry, gauss_shape_functions,
                std::tie(tke, TURBULENT_KINETIC_ENERGY),
                std::tie(wall_velocity, VELOCITY));

            const double wall_velocity_magnitude = norm_2(wall_velocity);

            double y_plus{0.0}, u_tau{0.0};
            CalculateYPlusAndUtau(y_plus, u_tau, wall_velocity_magnitude,
                                  mWallHeight, nu, kappa, beta);
            y_plus = std::max(y_plus, y_plus_limit);

            condition_y_plus += y_plus;

            if (wall_velocity_magnitude > eps) {
                const double u_tau = SoftMax(
                    c_mu_25 * std::sqrt(std::max(tke, 0.0)),
                    wall_velocity_magnitude / (inv_kappa * std::log(y_plus) + beta));

                noalias(condition_u_tau) += wall_velocity * u_tau / wall_velocity_magnitude;

                const double value = rho * std::pow(u_tau, 2) *
                                     gauss_weights[g] / wall_velocity_magnitude;

                for (IndexType a = 0; a < r_geometry.PointsNumber(); ++a) {
                    for (IndexType dim = 0; dim < TDim; ++dim) {
                        for (IndexType b = 0; b < r_geometry.PointsNumber(); ++b) {
                            rLocalMatrix(a * block_size + dim, b * block_size + dim) +=
                                gauss_shape_functions[a] * gauss_shape_functions[b] * value;
                        }
                        rLocalVector[a * block_size + dim] -=
                            gauss_shape_functions[a] * value * wall_velocity[dim];
                    }
                }
            }
        }

        const double inv_number_of_gauss_points =
            static_cast<double>(1.0 / num_gauss_points);
        this->SetValue(RANS_Y_PLUS, condition_y_plus * inv_number_of_gauss_points);
        this->SetValue(FRICTION_VELOCITY, condition_u_tau * inv_number_of_gauss_points);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

// template instantiations

template class VMSMonolithicKBasedWallCondition<2, 2>;
template class VMSMonolithicKBasedWallCondition<3, 3>;

} // namespace Kratos.
