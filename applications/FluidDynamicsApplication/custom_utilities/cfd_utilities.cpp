//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes
#include <functional>

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "utilities/variable_utils.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

// Include base h
#include "cfd_utilities.h"

namespace Kratos
{
void CFDUtilities::ReactionBasedYPlus::CalculateData(ModelPart::ConditionType& rCondition)
{
}

void CFDUtilities::ReactionBasedYPlus::CalculateYPlusAndUTau(
    double& rYPlus, array_1d<double, 3>& rUTau, const array_1d<double, 3>& rNormal)
{
    rUTau = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < rGeometry.PointsNumber(); ++i_node)
    {
        const NodeType& r_node = rGeometry[i_node];
        noalias(rUTau) += r_node.FastGetSolutionStepValue(rReactionVariable) /
                          r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS);
    }
    // calculate shear stress
    noalias(rUTau) = rUTau / norm_2(rNormal);
    const double reaction_magnitude = norm_2(rUTau);

    noalias(rUTau) = rUTau / std::sqrt(reaction_magnitude * Density);
    rYPlus = norm_2(rUTau) * WallHeight / KinematicViscosity;
}

void CFDUtilities::CalculateConditionGeometryData(const GeometryType& rGeometry,
                                                  const GeometryData::IntegrationMethod& rIntegrationMethod,
                                                  Vector& rGaussWeights,
                                                  Matrix& rNContainer)
{
    const GeometryType::IntegrationPointsArrayType& integration_points =
        rGeometry.IntegrationPoints(rIntegrationMethod);

    const std::size_t number_of_integration_points = integration_points.size();
    const int dimension = rGeometry.WorkingSpaceDimension();
    const double domain_size = rGeometry.DomainSize();

    if (rGaussWeights.size() != number_of_integration_points)
    {
        rGaussWeights.resize(number_of_integration_points, false);
    }

    rNContainer = rGeometry.ShapeFunctionsValues(rIntegrationMethod);

    // CAUTION: "Jacobian" is 2.0*A for triangles but 0.5*A for lines
    double det_J = (dimension == 2) ? 0.5 * domain_size : 2.0 * domain_size;

    for (unsigned int g = 0; g < number_of_integration_points; g++)
    {
        rGaussWeights[g] = det_J * integration_points[g].Weight();
    }
}

template <>
void CFDUtilities::CalculateNormal<2>(array_1d<double, 3>& rNormal, const ConditionType& rCondition)
{
    const GeometryType& pGeometry = rCondition.GetGeometry();

    rNormal[0] = pGeometry[1].Y() - pGeometry[0].Y();
    rNormal[1] = -(pGeometry[1].X() - pGeometry[0].X());
    rNormal[2] = 0.00;
}

template <>
void CFDUtilities::CalculateNormal<3>(array_1d<double, 3>& rNormal, const ConditionType& rCondition)
{
    const GeometryType& pGeometry = rCondition.GetGeometry();

    array_1d<double, 3> v1, v2;
    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(rNormal, v1, v2);
    rNormal *= 0.5;
}

double CFDUtilities::CalculateWallHeight(const ConditionType& rCondition,
                                         const array_1d<double, 3>& rNormal)
{
    KRATOS_TRY

    array_1d<double, 3> normal = rNormal / norm_2(rNormal);

    const ElementType& r_parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0];

    const GeometryType& r_parent_geometry = r_parent_element.GetGeometry();
    const GeometryType& r_condition_geometry = rCondition.GetGeometry();

    auto calculate_cell_center = [](const GeometryType& rGeometry) -> array_1d<double, 3> {
        const int number_of_nodes = rGeometry.PointsNumber();
        array_1d<double, 3> cell_center = ZeroVector(3);
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            noalias(cell_center) =
                cell_center + rGeometry[i_node].Coordinates() *
                                  (1.0 / static_cast<double>(number_of_nodes));
        }

        return cell_center;
    };

    const array_1d<double, 3>& parent_center = calculate_cell_center(r_parent_geometry);

    const array_1d<double, 3>& condition_center =
        calculate_cell_center(r_condition_geometry);

    return inner_prod(condition_center - parent_center, normal);

    KRATOS_CATCH("");
}

template <>
double CFDUtilities::EvaluateInPoint<double>(const GeometryType& rGeometry,
                                             const Variable<double>& rVariable,
                                             const Vector& rShapeFunction,
                                             const int Step)
{
    const unsigned int number_of_nodes = rGeometry.PointsNumber();
    double value = 0.0;
    for (unsigned int c = 0; c < number_of_nodes; ++c)
    {
        value += rShapeFunction[c] * rGeometry[c].FastGetSolutionStepValue(rVariable, Step);
    }

    return value;
}

template <>
array_1d<double, 3> CFDUtilities::EvaluateInPoint<array_1d<double, 3>>(
    const GeometryType& rGeometry,
    const Variable<array_1d<double, 3>>& rVariable,
    const Vector& rShapeFunction,
    const int Step)
{
    const unsigned int number_of_nodes = rGeometry.PointsNumber();
    array_1d<double, 3> value = ZeroVector(3);
    for (unsigned int c = 0; c < number_of_nodes; ++c)
    {
        value += rShapeFunction[c] * rGeometry[c].FastGetSolutionStepValue(rVariable, Step);
    }

    return value;
}

void CFDUtilities::CalculateNumberOfNeighbourConditions(ModelPart& rModelPart)
{
    ModelPart::NodesContainerType& r_nodes = rModelPart.Nodes();

    // reseting all nodal variables which needs to be calculated
    VariableUtils().SetNonHistoricalVariableToZero(NUMBER_OF_NEIGHBOUR_CONDITIONS, r_nodes);

    // updating
    const int number_of_conditions = rModelPart.NumberOfConditions();
#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
    {
        ConditionType& r_cond = *(rModelPart.ConditionsBegin() + i_cond);
        ConditionType::GeometryType& r_geometry = r_cond.GetGeometry();
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
        {
            NodeType& r_node = r_geometry[i_node];
            r_node.SetLock();
            r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS) += 1;
            r_node.UnSetLock();
        }
    }

    rModelPart.GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_CONDITIONS);
}

void CFDUtilities::CalculateReactionBasedUTau(double& YPlus,
                                              array_1d<double, 3>& rUTau,
                                              const GeometryType& rGeometry,
                                              const double Density,
                                              const double KinematicViscosity,
                                              const double WallHeight,
                                              const double Area,
                                              const Variable<array_1d<double, 3>>& rReactionVariable)
{
}

void CFDUtilities::CalculateWallFunctionBasedYPlusUTau(double& YPlus,
                                                       array_1d<double, 3>& rUTau,
                                                       const GeometryType& rGeometry,
                                                       const double Density,
                                                       const double KinematicViscosity,
                                                       const double WallHeight,
                                                       const double VonKarman,
                                                       const double Smoothness,
                                                       const Vector& GaussPointShapeFunctions)
{
    const array_1d<double, 3>& r_wall_velocity =
        CFDUtilities::EvaluateInPoint(rGeometry, VELOCITY, GaussPointShapeFunctions);
    const array_1d<double, 3>& r_mesh_velocity = CFDUtilities::EvaluateInPoint(
        rGeometry, MESH_VELOCITY, GaussPointShapeFunctions);
    const array_1d<double, 3>& r_effective_velocity = r_wall_velocity - r_mesh_velocity;
    const double wall_velocity_magnitude = norm_2(r_effective_velocity);

    double u_tau{0.0};
    CFDUtilities::CalculateYPlusAndUtau(YPlus, u_tau, wall_velocity_magnitude, WallHeight,
                                        KinematicViscosity, VonKarman, Smoothness);

    if (wall_velocity_magnitude > 0.0)
    {
        noalias(rUTau) = r_effective_velocity * u_tau / wall_velocity_magnitude;
    }
    else
    {
        rUTau = ZeroVector(3);
    }
}

void CFDUtilities::CalculateReactionBasedYPlus(ModelPart& rModelPart,
                                               const Variable<double>& rKinematicViscosityVariable,
                                               const Variable<array_1d<double, 3>>& rReactionVariable)
{
    ModelPart::NodesContainerType& r_nodes = rModelPart.Nodes();

    // reseting all nodal variables which needs to be calculated
    VariableUtils().SetNonHistoricalVariableToZero(Y_PLUS, r_nodes);
    VariableUtils().SetNonHistoricalVariableToZero(FRICTION_VELOCITY, r_nodes);

    // calculate number of neighour conditions
    CFDUtilities::CalculateNumberOfNeighbourConditions(rModelPart);

    const int domain_size = rModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const std::function<void(array_1d<double, 3>&, const ConditionType&)> normal_calculation_method =
        (domain_size == 2) ? CFDUtilities::CalculateNormal<2>
                           : CFDUtilities::CalculateNormal<3>;

    const int number_of_conditions = rModelPart.NumberOfConditions();
#pragma omp parallel
    {
        array_1d<double, 3> normal;
        array_1d<double, 3> reaction;
        array_1d<double, 3> u_tau;
#pragma omp for
        for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
        {
            ConditionType& r_condition = *(rModelPart.ConditionsBegin() + i_cond);
            GeometryType& r_geometry = r_condition.GetGeometry();

            Vector gauss_weights;
            Matrix shape_functions;
            CFDUtilities::CalculateConditionGeometryData(
                r_condition.GetGeometry(), GeometryData::GI_GAUSS_1,
                gauss_weights, shape_functions);

            const Vector& gauss_shape_functions = row(shape_functions, 0);

            // calculate normal for the condition
            normal_calculation_method(normal, r_condition);
            const double normal_area = norm_2(normal);

            const double y_wall = CFDUtilities::CalculateWallHeight(r_condition, normal);
            const double nu = CFDUtilities::EvaluateInPoint(
                r_geometry, rKinematicViscosityVariable, gauss_shape_functions);
            const double rho = CFDUtilities::EvaluateInPoint(
                r_geometry, DENSITY, gauss_shape_functions);

            // calculate appropriate reaction force for condition
            reaction = ZeroVector(3);
            for (unsigned int i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
            {
                const NodeType& r_node = r_geometry[i_node];
                noalias(reaction) += r_node.FastGetSolutionStepValue(rReactionVariable) /
                                     r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS);
            }
            // calculate shear stress
            noalias(reaction) = reaction / normal_area;
            const double reaction_magnitude = norm_2(reaction);

            noalias(u_tau) = reaction / std::sqrt(reaction_magnitude * rho);
            const double y_plus = norm_2(u_tau) * y_wall / nu;

            r_condition.SetValue(FRICTION_VELOCITY, u_tau);
            r_condition.SetValue(Y_PLUS, y_plus);

            // distribute calculated values to nodes of the condition
            for (unsigned int i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
            {
                NodeType& r_node = r_geometry[i_node];

                const double nodal_y_plus = r_node.GetValue(Y_PLUS);
                const array_1d<double, 3>& nodal_u_tau = r_node.GetValue(FRICTION_VELOCITY);

                r_node.SetLock();
                r_node.SetValue(Y_PLUS, nodal_y_plus + y_plus);
                r_node.SetValue(FRICTION_VELOCITY, nodal_u_tau + u_tau);
                r_node.UnSetLock();
            }
        }
    }

    rModelPart.GetCommunicator().AssembleNonHistoricalData(Y_PLUS);
    rModelPart.GetCommunicator().AssembleNonHistoricalData(FRICTION_VELOCITY);

    const int number_of_nodes = r_nodes.size();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(rModelPart.NodesBegin() + i_node);
        const double inv_number_of_neighbour_conditions =
            1.0 / static_cast<double>(r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS));

        const double nodal_y_plus = r_node.GetValue(Y_PLUS);
        const array_1d<double, 3>& nodal_u_tau = r_node.GetValue(FRICTION_VELOCITY);

        r_node.SetValue(Y_PLUS, nodal_y_plus * inv_number_of_neighbour_conditions);
        r_node.SetValue(FRICTION_VELOCITY, nodal_u_tau * inv_number_of_neighbour_conditions);
    }
}

double CFDUtilities::CalculateLogarithmicYPlusLimit(const double VonKarman,
                                                    const double Smoothness,
                                                    const int MaxIterations,
                                                    const double Tolerance)
{
    double y_plus = 11.06;
    const double inv_kappa = 1.0 / VonKarman;
    double dx = 0.0;
    for (int i = 0; i < MaxIterations; ++i)
    {
        const double value = inv_kappa * std::log(y_plus) + Smoothness;
        dx = value - y_plus;

        if (std::abs(dx) < Tolerance)
            return y_plus;

        y_plus = value;
    }

    KRATOS_WARNING("LogarithmicYPlusLimit")
        << "Logarithmic y_plus limit reached max iterations with dx > "
           "Tolerance [ "
        << dx << " > " << Tolerance << ", MaxIterations = " << MaxIterations << " ].\n";
    return y_plus;
}

void CFDUtilities::CalculateYPlusAndUtau(double& rYPlus,
                                         double& rUTau,
                                         const double WallVelocity,
                                         const double WallHeight,
                                         const double KinematicViscosity,
                                         const double VonKarman,
                                         const double Smoothness,
                                         const int MaxIterations,
                                         const double Tolerance)
{
    const double limit_y_plus = CalculateLogarithmicYPlusLimit(
        VonKarman, Smoothness, MaxIterations, Tolerance);

    // linear region
    rUTau = std::sqrt(WallVelocity * KinematicViscosity / WallHeight);
    rYPlus = rUTau * WallHeight / KinematicViscosity;
    const double inv_kappa = 1.0 / VonKarman;

    // log region
    if (rYPlus > limit_y_plus)
    {
        int iter = 0;
        double dx = 1e10;
        double u_plus = inv_kappa * std::log(rYPlus) + Smoothness;

        while (iter < MaxIterations && std::fabs(dx) > Tolerance * rUTau)
        {
            // Newton-Raphson iteration
            double f = rUTau * u_plus - WallVelocity;
            double df = u_plus + inv_kappa;
            dx = f / df;

            // Update variables
            rUTau -= dx;
            rYPlus = WallHeight * rUTau / KinematicViscosity;
            u_plus = inv_kappa * std::log(rYPlus) + Smoothness;
            ++iter;
        }
        if (iter == MaxIterations)
        {
            std::cout << "Warning: wall condition Newton-Raphson did not "
                         "converge. Residual is "
                      << dx << std::endl;
        }
    }
}

void CFDUtilities::CalculateWallFunctionBasedYPlus(ModelPart& rModelPart,
                                                   const double VonKarman = 0.41,
                                                   const double Smoothness = 5.2)
{
}

} // namespace Kratos.
