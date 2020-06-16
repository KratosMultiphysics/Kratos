//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_vms_monolithic_k_based_wall_condition.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
RansVMSMonolithicKBasedWallCondition<TDim, TNumNodes>& RansVMSMonolithicKBasedWallCondition<TDim, TNumNodes>::operator=(
    RansVMSMonolithicKBasedWallCondition<TDim, TNumNodes> const& rOther)
{
    Condition::operator=(rOther);

    return *this;
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer RansVMSMonolithicKBasedWallCondition<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<RansVMSMonolithicKBasedWallCondition>(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer RansVMSMonolithicKBasedWallCondition<TDim, TNumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<RansVMSMonolithicKBasedWallCondition>(
        NewId, pGeom, pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer RansVMSMonolithicKBasedWallCondition<TDim, TNumNodes>::Clone(
    IndexType NewId, NodesArrayType const& rThisNodes) const
{
    Condition::Pointer pNewCondition = Create(
        NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

    pNewCondition->SetData(this->GetData());
    pNewCondition->SetFlags(this->GetFlags());

    return pNewCondition;
}

template <unsigned int TDim, unsigned int TNumNodes>
int RansVMSMonolithicKBasedWallCondition<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int Check = BaseType::Check(rCurrentProcessInfo);

    if (Check != 0)
    {
        return Check;
    }

    const GeometryType& r_geometry = this->GetGeometry();

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
    {
        const NodeType& r_node = r_geometry[i_node];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
    }

    return 0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansVMSMonolithicKBasedWallCondition<TDim, TNumNodes>::Initialize()
{
    KRATOS_TRY;

    if (RansCalculationUtilities::IsWall(*this))
    {
        const array_1d<double, 3>& rNormal = this->GetValue(NORMAL);
        KRATOS_ERROR_IF(norm_2(rNormal) == 0.0)
            << "NORMAL must be calculated before using this " << this->Info() << "\n";

        KRATOS_ERROR_IF(this->GetValue(NEIGHBOUR_ELEMENTS).size() == 0)
            << this->Info() << " cannot find parent element\n";

        mWallHeight = RansCalculationUtilities::CalculateWallHeight(*this, rNormal);

        KRATOS_ERROR_IF(mWallHeight == 0.0) << this->Info() << " has zero wall height.\n";
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string RansVMSMonolithicKBasedWallCondition<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "RansVMSMonolithicKBasedWallCondition" << TDim << "D";
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansVMSMonolithicKBasedWallCondition<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "RansVMSMonolithicKBasedWallCondition";
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansVMSMonolithicKBasedWallCondition<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansVMSMonolithicKBasedWallCondition<TDim, TNumNodes>::ApplyWallLaw(
    MatrixType& rLocalMatrix, VectorType& rLocalVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (RansCalculationUtilities::IsWall(*this))
    {
        const GeometryType& r_geometry = this->GetGeometry();
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        RansCalculationUtilities::CalculateConditionGeometryData(
            r_geometry, this->GetIntegrationMethod(), gauss_weights, shape_functions);
        const IndexType num_gauss_points = gauss_weights.size();

        const size_t block_size = TDim + 1;

        const double c_mu_25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
        const double kappa = rCurrentProcessInfo[WALL_VON_KARMAN];
        const double inv_kappa = 1.0 / kappa;
        const double beta = rCurrentProcessInfo[WALL_SMOOTHNESS_BETA];
        const double y_plus_limit = rCurrentProcessInfo[RANS_Y_PLUS_LIMIT];

        const double eps = std::numeric_limits<double>::epsilon();

        double condition_y_plus{0.0};
        array_1d<double, 3> condition_u_tau = ZeroVector(3);

        for (size_t g = 0; g < num_gauss_points; ++g)
        {
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const array_1d<double, 3>& r_wall_velocity =
                RansCalculationUtilities::EvaluateInPoint(
                    r_geometry, VELOCITY, gauss_shape_functions);
            const double wall_velocity_magnitude = norm_2(r_wall_velocity);

            const double tke = RansCalculationUtilities::EvaluateInPoint(
                r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
            const double rho = RansCalculationUtilities::EvaluateInPoint(
                r_geometry, DENSITY, gauss_shape_functions);
            const double nu = RansCalculationUtilities::EvaluateInPoint(
                r_geometry, KINEMATIC_VISCOSITY, gauss_shape_functions);

            double y_plus{0.0}, u_tau{0.0};
            RansCalculationUtilities::CalculateYPlusAndUtau(
                y_plus, u_tau, wall_velocity_magnitude, mWallHeight, nu, kappa, beta);
            y_plus = std::max(y_plus, y_plus_limit);

            condition_y_plus += y_plus;

            if (wall_velocity_magnitude > eps)
            {
                const double u_tau = RansCalculationUtilities::SoftMax(
                    c_mu_25 * std::sqrt(std::max(tke, 0.0)),
                    wall_velocity_magnitude / (inv_kappa * std::log(y_plus) + beta));

                noalias(condition_u_tau) -= r_wall_velocity * u_tau / wall_velocity_magnitude;

                const double value = rho * std::pow(u_tau, 2) *
                                     gauss_weights[g] / wall_velocity_magnitude;

                for (size_t a = 0; a < r_geometry.PointsNumber(); ++a)
                {
                    for (size_t dim = 0; dim < TDim; ++dim)
                    {
                        for (size_t b = 0; b < r_geometry.PointsNumber(); ++b)
                        {
                            rLocalMatrix(a * block_size + dim, b * block_size + dim) +=
                                gauss_shape_functions[a] * gauss_shape_functions[b] * value;
                        }
                        rLocalVector[a * block_size + dim] -=
                            gauss_shape_functions[a] * value * r_wall_velocity[dim];
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
void RansVMSMonolithicKBasedWallCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansVMSMonolithicKBasedWallCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

// template instantiations

template class RansVMSMonolithicKBasedWallCondition<2, 2>;
template class RansVMSMonolithicKBasedWallCondition<3, 3>;

} // namespace Kratos.
