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
#include <iostream>
#include <limits>
#include <string>
#include <cmath>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_evm_k_epsilon_vms_monolithic_wall.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>& RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>::operator=(
    RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes> const& rOther)
{
    Condition::operator=(rOther);

    return *this;
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<RansEvmKEpsilonVmsMonolithicWall>(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<RansEvmKEpsilonVmsMonolithicWall>(NewId, pGeom, pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>::Clone(
    IndexType NewId, NodesArrayType const& rThisNodes) const
{
    Condition::Pointer pNewCondition = Create(
        NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

    pNewCondition->SetData(this->GetData());
    pNewCondition->SetFlags(this->GetFlags());

    return pNewCondition;
}

template <unsigned int TDim, unsigned int TNumNodes>
int RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
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

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_Y_PLUS, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
    }

    return 0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "RansEvmKEpsilonVmsMonolithicWall" << TDim << "D";
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "RansEvmKEpsilonVmsMonolithicWall";
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>::ApplyLogarithmicWallLaw(
    MatrixType& rLocalMatrix, VectorType& rLocalVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeometryType& rGeometry = this->GetGeometry();
    const size_t BlockSize = TDim + 1;
    const double NodalFactor = 1.0 / double(TDim);

    double area = NodalFactor * rGeometry.DomainSize();
    // DomainSize() is the way to ask the geometry's length/area/volume (whatever is relevant for its dimension) without asking for the number of spatial dimensions first

    for (size_t itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
    {
        const NodeType& rConstNode = rGeometry[itNode];
        const double y = rConstNode.FastGetSolutionStepValue(DISTANCE); // wall distance to use in stress calculation
        if (y > 0.0 && rConstNode.Is(SLIP))
        {
            array_1d<double, 3> Vel = rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3>& VelMesh =
                rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
            Vel -= VelMesh;
            const double Ikappa = 1.0 / 0.41; // inverse of Von Karman's kappa
            const double B = 5.2;
            const double limit_yplus = 10.9931899; // limit between linear and log regions

            const double rho = rGeometry[itNode].FastGetSolutionStepValue(DENSITY);
            const double nu = rGeometry[itNode].FastGetSolutionStepValue(VISCOSITY);

            double wall_vel = 0.0;
            for (size_t d = 0; d < TDim; ++d)
            {
                wall_vel += Vel[d] * Vel[d];
            }
            wall_vel = sqrt(wall_vel);

            if (wall_vel > 1e-12) // do not bother if velocity is zero
            {
                // linear region
                double utau = sqrt(wall_vel * nu / y);
                double yplus = y * utau / nu;

                // log region
                if (yplus > limit_yplus)
                {
                    // wall_vel / utau = 1/kappa * log(yplus) + B
                    // this requires solving a nonlinear problem:
                    // f(utau) = utau*(1/kappa * log(y*utau/nu) + B) - wall_vel = 0
                    // note that f'(utau) = 1/kappa * log(y*utau/nu) + B + 1/kappa

                    IndexType iter = 0;
                    double dx = 1e10;
                    const double tol = 1e-6;
                    double uplus = Ikappa * log(yplus) + B;

                    while (iter < 100 && fabs(dx) > tol * utau)
                    {
                        // Newton-Raphson iteration
                        double f = utau * uplus - wall_vel;
                        double df = uplus + Ikappa;
                        dx = f / df;

                        // Update variables
                        utau -= dx;
                        yplus = y * utau / nu;
                        uplus = Ikappa * log(yplus) + B;
                        ++iter;
                    }
                    if (iter == 100)
                    {
                        std::cout
                            << "Warning: wall condition Newton-Raphson did "
                               "not converge. Residual is "
                            << dx << std::endl;
                    }
                }
                const double Tmp = area * utau * utau * rho / wall_vel;
                for (size_t d = 0; d < TDim; ++d)
                {
                    size_t k = itNode * BlockSize + d;
                    rLocalVector[k] -= Vel[d] * Tmp;
                    rLocalMatrix(k, k) += Tmp;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>::ApplyRansBasedWallLaw(
    MatrixType& rLocalMatrix, VectorType& rLocalVector, ProcessInfo& rCurrentProcessInfo)

{
    KRATOS_TRY

    GeometryType& r_geometry = this->GetGeometry();

    const GeometryType::IntegrationPointsArrayType& integration_points =
        r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const std::size_t number_of_gauss_points = integration_points.size();
    MatrixType shape_functions = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

    array_1d<double, 3> normal;
    this->CalculateNormal(normal); // this already contains the area
    double A = norm_2(normal);

    // CAUTION: "Jacobian" is 2.0*A for triangles but 0.5*A for lines
    double J = (TDim == 2) ? 0.5 * A : 2.0 * A;

    const size_t block_size = TDim + 1;

    const double c_mu_25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
    const double inv_von_karman = 1.0 / rCurrentProcessInfo[WALL_VON_KARMAN];
    const double beta = rCurrentProcessInfo[WALL_SMOOTHNESS_BETA];

    const double eps = std::numeric_limits<double>::epsilon();

    for (size_t g = 0; g < number_of_gauss_points; ++g)
    {
        const Vector& gauss_shape_functions = row(shape_functions, g);
        const double weight = J * integration_points[g].Weight();

        const array_1d<double, 3>& r_wall_velocity = RansCalculationUtilities::EvaluateInPoint(
            r_geometry, VELOCITY, gauss_shape_functions);
        const double wall_velocity_magnitude = norm_2(r_wall_velocity);

        const double tke = RansCalculationUtilities::EvaluateInPoint(
            r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
        const double y_plus = RansCalculationUtilities::EvaluateInPoint(
            r_geometry, RANS_Y_PLUS, gauss_shape_functions);
        const double rho = RansCalculationUtilities::EvaluateInPoint(
            r_geometry, DENSITY, gauss_shape_functions);

        if (wall_velocity_magnitude > eps)
        {
            const double u_tau = RansCalculationUtilities::SoftMax(
                c_mu_25 * std::sqrt(std::max(tke, 0.0)),
                wall_velocity_magnitude / (inv_von_karman * std::log(y_plus) + beta));
            const double value = rho * std::pow(u_tau, 2) * weight / wall_velocity_magnitude;

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

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>::ApplyWallLaw(
    MatrixType& rLocalMatrix, VectorType& rLocalVector, ProcessInfo& rCurrentProcessInfo)
{
    if (this->Is(SLIP))
    {
        if (rCurrentProcessInfo[IS_CO_SOLVING_PROCESS_ACTIVE])
        {
            this->ApplyRansBasedWallLaw(rLocalMatrix, rLocalVector, rCurrentProcessInfo);
        }
        else
        {
            this->ApplyLogarithmicWallLaw(rLocalMatrix, rLocalVector, rCurrentProcessInfo);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmKEpsilonVmsMonolithicWall<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

// template instantiations

template class RansEvmKEpsilonVmsMonolithicWall<2, 2>;
template class RansEvmKEpsilonVmsMonolithicWall<3, 3>;

} // namespace Kratos.
