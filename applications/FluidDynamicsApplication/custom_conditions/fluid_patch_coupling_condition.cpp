//     |  /           |
//     ' /   __| _` | __|  _ \   __|
//     . \  |   (   | |   (   |\__ `
//    _|\_\_|  \__,_|\__|\___/ ____/
//                    Multi-Physics
//
//   Fluid Dynamics Application
//
//   License:         BSD License
//                    Kratos default license: kratos/license.txt
//
//   Main author:     Nicolò Antonelli
//

// System includes

// Project includes
#include "custom_conditions/fluid_patch_coupling_condition.h"
#include "includes/checks.h"
#include "includes/element.h"
#include "integration/line_gauss_legendre_integration_points.h"
#include "utilities/element_size_calculator.h"

namespace Kratos
{
namespace
{
Point interface_point(const Matrix& rCoordinates, const double Xi)
{
    return Point(
        0.5 * ((1.0 - Xi) * rCoordinates(0, 0) + (1.0 + Xi) * rCoordinates(1, 0)),
        0.5 * ((1.0 - Xi) * rCoordinates(0, 1) + (1.0 + Xi) * rCoordinates(1, 1)), 0.0);
}
} // unnamed namespace

Condition::Pointer FluidPatchCouplingCondition2D::Create(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<FluidPatchCouplingCondition2D>(NewId, pGeometry, pProperties);
}

Condition::Pointer FluidPatchCouplingCondition2D::Create(
    IndexType NewId, const NodesArrayType& rNodes, PropertiesType::Pointer pProperties) const
{
    return Create(NewId, Kratos::make_shared<Geometry<Node>>(rNodes), pProperties);
}

void FluidPatchCouplingCondition2D::CalculateCouplingMatrix(
    MatrixType& rLHS, const ProcessInfo& rProcessInfo) const
{
    KRATOS_TRY
    Check(rProcessInfo);
    const IndexType parent_nodes = GetValue(NEIGHBOUR_ELEMENTS)[0].GetGeometry().PointsNumber();
    const IndexType local_size = 3 * GetGeometry().PointsNumber();
    rLHS = ZeroMatrix(local_size, local_size);
    const auto& r_parents = GetValue(NEIGHBOUR_ELEMENTS);
    const auto& r_coordinates = GetValue(COUPLING_INTERFACE_COORDINATES);
    const double dx = r_coordinates(1, 0) - r_coordinates(0, 0);
    const double dy = r_coordinates(1, 1) - r_coordinates(0, 1);
    const double length = std::hypot(dx, dy);
    const array_1d<double, 2> normal{dy / length, -dx / length};
    const double coefficient = GetProperties()[PENALTY_COEFFICIENT];
    const double theta = coefficient > 0.0 ? 1.0 : -1.0;

    std::array<BoundedMatrix<double, 3, 3>, 2> constitutive_matrices;
    double penalty_scale = 0.0;
    for (IndexType side = 0; side < 2; ++side) {
        const auto& r_parent = r_parents[side];
        const double mu = r_parent.GetProperties()[DYNAMIC_VISCOSITY];
        auto& r_c = constitutive_matrices[side];
        noalias(r_c) = ZeroMatrix(3, 3);
        // Plane representation of Newtonian2DLaw's 3D deviatoric stress.
        r_c(0, 0) = r_c(1, 1) = 4.0 * mu / 3.0;
        r_c(0, 1) = r_c(1, 0) = -2.0 * mu / 3.0;
        r_c(2, 2) = mu;
        const double h = parent_nodes == 3
            ? ElementSizeCalculator<2, 3>::MinimumElementSize(r_parent.GetGeometry())
            : ElementSizeCalculator<2, 4>::MinimumElementSize(r_parent.GetGeometry());
        penalty_scale = std::max(penalty_scale, mu / h);
    }
    const double alpha = std::max(coefficient, 0.0) * penalty_scale; // P1 and Q1 have order one

    // Two-point line quadrature integrates trace products on each common
    // subsegment. It is exact for P1 and affine Q1 parents; distorted Q1
    // parents use the physical gradients evaluated at each Gauss point.
    for (const auto& r_point : LineGaussLegendreIntegrationPoints2::IntegrationPoints()) {
        const Point global_point = interface_point(r_coordinates, r_point.X());
        const double weight = 0.5 * length * r_point.Weight();
        Matrix jump = ZeroMatrix(2, local_size);
        Matrix traction = ZeroMatrix(2, local_size);
        for (IndexType side = 0; side < 2; ++side) {
            const auto& r_geometry = r_parents[side].GetGeometry();
            GeometryType::CoordinatesArrayType local_coordinates = ZeroVector(3);
            r_geometry.PointLocalCoordinates(local_coordinates, global_point);
            Vector shape_values;
            r_geometry.ShapeFunctionsValues(shape_values, local_coordinates);
            // Q1 gradients vary in physical space on a general quadrilateral.
            Matrix local_gradients, inverse_jacobian;
            r_geometry.ShapeFunctionsLocalGradients(local_gradients, local_coordinates);
            KRATOS_ERROR_IF(r_geometry.DeterminantOfJacobian(local_coordinates) <= 0.0)
                << "Parent elements must have positive Jacobians at interface Gauss points." << std::endl;
            r_geometry.InverseOfJacobian(inverse_jacobian, local_coordinates);
            const Matrix gradients = prod(local_gradients, inverse_jacobian);
            const auto& r_c = constitutive_matrices[side];
            const double sign = side == 0 ? 1.0 : -1.0;
            for (IndexType i = 0; i < parent_nodes; ++i) {
                const IndexType offset = 3 * (parent_nodes * side + i);
                const double gx = gradients(i, 0);
                const double gy = gradients(i, 1);
                jump(0, offset) = jump(1, offset + 1) = sign * shape_values[i];
                // Each traction column is (C B_i)n/2; pressure contributes
                // -N_i n/2. Both stencils include all off-interface nodes.
                for (IndexType component = 0; component < 2; ++component) {
                    const array_1d<double, 3> strain = component == 0
                        ? array_1d<double, 3>{gx, 0.0, gy}
                        : array_1d<double, 3>{0.0, gy, gx};
                    const array_1d<double, 3> stress = prod(r_c, strain);
                    traction(0, offset + component) = 0.5 * (stress[0] * normal[0] + stress[2] * normal[1]);
                    traction(1, offset + component) = 0.5 * (stress[2] * normal[0] + stress[1] * normal[1]);
                    traction(component, offset + 2) = -0.5 * shape_values[i] * normal[component];
                }
            }
        }
        // -J^T T - theta T^T J + alpha J^T J. Keeping the pressure columns
        // in T also supplies the q/[u] terms; pressure is not penalized.
        for (IndexType i = 0; i < local_size; ++i) {
            for (IndexType j = 0; j < local_size; ++j) {
                // The new Stokes gap path uses q div(u) in the bulk. Its
                // pressure test must retain -{q}[u].n for either velocity
                // adjoint choice; otherwise pressure/velocity work doubles.
                const double adjoint_theta = UsesBackgroundExtension() && i % 3 == 2 ? -1.0 : theta;
                for (IndexType d = 0; d < 2; ++d) {
                    rLHS(i, j) += weight * (-jump(d, i) * traction(d, j)
                        - adjoint_theta * traction(d, i) * jump(d, j)
                        + alpha * jump(d, i) * jump(d, j));
                }
            }
        }
    }
    KRATOS_CATCH("")
}

void FluidPatchCouplingCondition2D::CalculateResidual(const MatrixType& rLHS, VectorType& rRHS) const
{
    Vector values;
    GetFirstDerivativesVector(values);
    rRHS.resize(3 * GetGeometry().PointsNumber(), false);
    noalias(rRHS) = -prod(rLHS, values);
}

void FluidPatchCouplingCondition2D::CalculateLocalSystem(
    MatrixType& rLHS, VectorType& rRHS, const ProcessInfo& rProcessInfo)
{
    // Match the fluid application's split between internally integrated BDF
    // elements and externally integrated (e.g. Bossak) damping contributions.
    if (rProcessInfo.Has(BDF_COEFFICIENTS)) {
        CalculateCouplingMatrix(rLHS, rProcessInfo);
        CalculateResidual(rLHS, rRHS);
    } else {
        rLHS = ZeroMatrix(3 * GetGeometry().PointsNumber(), 3 * GetGeometry().PointsNumber());
        rRHS = ZeroVector(3 * GetGeometry().PointsNumber());
    }
}

void FluidPatchCouplingCondition2D::CalculateLeftHandSide(MatrixType& rLHS, const ProcessInfo& rProcessInfo)
{
    Vector rhs;
    CalculateLocalSystem(rLHS, rhs, rProcessInfo);
}

void FluidPatchCouplingCondition2D::CalculateRightHandSide(VectorType& rRHS, const ProcessInfo& rProcessInfo)
{
    Matrix lhs;
    CalculateLocalSystem(lhs, rRHS, rProcessInfo);
}

void FluidPatchCouplingCondition2D::CalculateLocalVelocityContribution(
    MatrixType& rLHS, VectorType& rRHS, const ProcessInfo& rProcessInfo)
{
    if (rProcessInfo.Has(BDF_COEFFICIENTS)) {
        rLHS = ZeroMatrix(3 * GetGeometry().PointsNumber(), 3 * GetGeometry().PointsNumber());
        // Preserve the residual already supplied by CalculateLocalSystem.
        if (rRHS.size() != 3 * GetGeometry().PointsNumber()) {
            rRHS = ZeroVector(3 * GetGeometry().PointsNumber());
        }
    } else {
        CalculateCouplingMatrix(rLHS, rProcessInfo);
        CalculateResidual(rLHS, rRHS);
    }
}

void FluidPatchCouplingCondition2D::CalculateDampingMatrix(MatrixType& rLHS, const ProcessInfo& rProcessInfo)
{
    if (rProcessInfo.Has(BDF_COEFFICIENTS)) {
        rLHS = ZeroMatrix(3 * GetGeometry().PointsNumber(), 3 * GetGeometry().PointsNumber());
    } else {
        CalculateCouplingMatrix(rLHS, rProcessInfo);
    }
}

void FluidPatchCouplingCondition2D::EquationIdVector(
    EquationIdVectorType& rResult, const ProcessInfo& rProcessInfo) const
{
    rResult.resize(3 * GetGeometry().PointsNumber(), false);
    for (IndexType i = 0; i < GetGeometry().PointsNumber(); ++i) {
        rResult[3 * i] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
        rResult[3 * i + 1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
        rResult[3 * i + 2] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
    }
}

void FluidPatchCouplingCondition2D::GetDofList(DofsVectorType& rDofs, const ProcessInfo& rProcessInfo) const
{
    rDofs.resize(3 * GetGeometry().PointsNumber());
    for (IndexType i = 0; i < GetGeometry().PointsNumber(); ++i) {
        rDofs[3 * i] = GetGeometry()[i].pGetDof(VELOCITY_X);
        rDofs[3 * i + 1] = GetGeometry()[i].pGetDof(VELOCITY_Y);
        rDofs[3 * i + 2] = GetGeometry()[i].pGetDof(PRESSURE);
    }
}

void FluidPatchCouplingCondition2D::GetFirstDerivativesVector(Vector& rValues, const int Step) const
{
    rValues.resize(3 * GetGeometry().PointsNumber(), false);
    for (IndexType i = 0; i < GetGeometry().PointsNumber(); ++i) {
        const auto& r_node = GetGeometry()[i];
        const auto& r_velocity = r_node.FastGetSolutionStepValue(VELOCITY, Step);
        rValues[3 * i] = r_velocity[0];
        rValues[3 * i + 1] = r_velocity[1];
        rValues[3 * i + 2] = r_node.FastGetSolutionStepValue(PRESSURE, Step);
    }
}

void FluidPatchCouplingCondition2D::GetSecondDerivativesVector(Vector& rValues, const int Step) const
{
    rValues.resize(3 * GetGeometry().PointsNumber(), false);
    noalias(rValues) = ZeroVector(3 * GetGeometry().PointsNumber());
    for (IndexType i = 0; i < GetGeometry().PointsNumber(); ++i) {
        const auto& r_acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        rValues[3 * i] = r_acceleration[0];
        rValues[3 * i + 1] = r_acceleration[1];
    }
}

int FluidPatchCouplingCondition2D::Check(const ProcessInfo& rProcessInfo) const
{
    KRATOS_TRY
    KRATOS_ERROR_IF_NOT(Has(NEIGHBOUR_ELEMENTS)) << Info() << " requires NEIGHBOUR_ELEMENTS." << std::endl;
    const auto& r_parents = GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_ERROR_IF(r_parents.size() != 2) << Info() << " requires two parent elements." << std::endl;
    const IndexType parent_nodes = r_parents[0].GetGeometry().PointsNumber();
    KRATOS_ERROR_IF(UsesBackgroundExtension() && parent_nodes != 4)
        << "Gap-SBM coupling requires Q1 parents." << std::endl;
    KRATOS_ERROR_IF(GetGeometry().size() != 2 * parent_nodes ||
                    r_parents[1].GetGeometry().GetGeometryType() != r_parents[0].GetGeometry().GetGeometryType())
        << Info() << " requires two parents of the same type and their independent nodes." << std::endl;
    KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_COEFFICIENT))
        << Info() << " requires PENALTY_COEFFICIENT in its Properties." << std::endl;
    const double coefficient = GetProperties()[PENALTY_COEFFICIENT];
    KRATOS_ERROR_IF(!std::isfinite(coefficient) || (coefficient < 0.0 && coefficient != -1.0))
        << "PENALTY_COEFFICIENT must be non-negative or -1 (penalty-free)." << std::endl;
    KRATOS_ERROR_IF_NOT(Has(COUPLING_INTERFACE_COORDINATES))
        << Info() << " requires COUPLING_INTERFACE_COORDINATES." << std::endl;
    const auto& r_coordinates = GetValue(COUPLING_INTERFACE_COORDINATES);
    KRATOS_ERROR_IF(r_coordinates.size1() != 2 || r_coordinates.size2() != 3)
        << "COUPLING_INTERFACE_COORDINATES must be a 2x3 matrix." << std::endl;
    for (IndexType i = 0; i < 2; ++i) {
        for (IndexType d = 0; d < 3; ++d) {
            KRATOS_ERROR_IF(!std::isfinite(r_coordinates(i, d))) << "Non-finite interface coordinate." << std::endl;
        }
        KRATOS_ERROR_IF(std::abs(r_coordinates(i, 2)) > 1e-10) << "Only the z=0 plane is supported." << std::endl;
    }
    const double dx = r_coordinates(1, 0) - r_coordinates(0, 0);
    const double dy = r_coordinates(1, 1) - r_coordinates(0, 1);
    const double length = std::hypot(dx, dy);
    KRATOS_ERROR_IF(length <= 0.0) << "Zero-length coupling segment." << std::endl;
    const array_1d<double, 2> normal{dy / length, -dx / length};
    const Point midpoint = interface_point(r_coordinates, 0.0);
    std::set<IndexType> node_ids;
    for (IndexType side = 0; side < 2; ++side) {
        const auto& r_parent = r_parents[side];
        const auto& r_geometry = r_parent.GetGeometry();
        KRATOS_ERROR_IF(r_geometry.GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Triangle2D3 &&
                        r_geometry.GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
            << Info() << " requires Triangle2D3 (P1) or Quadrilateral2D4 (Q1) parents." << std::endl;
        KRATOS_ERROR_IF(r_geometry.Area() <= 0.0) << "Parent elements must have positive area." << std::endl;
        KRATOS_ERROR_IF_NOT(r_parent.GetProperties().Has(DYNAMIC_VISCOSITY))
            << "Each parent requires DYNAMIC_VISCOSITY." << std::endl;
        const double mu = r_parent.GetProperties()[DYNAMIC_VISCOSITY];
        KRATOS_ERROR_IF(!std::isfinite(mu) || mu <= 0.0) << "DYNAMIC_VISCOSITY must be positive." << std::endl;
        const auto center = r_geometry.Center();
        const double normal_distance = normal[0] * (center.X() - midpoint.X())
                                     + normal[1] * (center.Y() - midpoint.Y());
        KRATOS_ERROR_IF((side == 0 && !UsesBackgroundExtension() && normal_distance >= 0.0) || (side == 1 && normal_distance <= 0.0))
            << "The interface normal must point from parent A to parent B." << std::endl;
        std::array<Vector, 2> endpoint_values;
        for (IndexType endpoint = 0; endpoint < 2; ++endpoint) {
            GeometryType::CoordinatesArrayType local = ZeroVector(3);
            r_geometry.PointLocalCoordinates(local, Point(
                r_coordinates(endpoint, 0), r_coordinates(endpoint, 1), 0.0));
            r_geometry.ShapeFunctionsValues(endpoint_values[endpoint], local);
            for (const double value : endpoint_values[endpoint]) {
                KRATOS_ERROR_IF(!(side == 0 && UsesBackgroundExtension()) &&
                                (value < -1e-9 || value > 1.0 + 1e-9))
                    << "Coupling endpoints must lie on both parent edges; gaps are not supported." << std::endl;
            }
        }
        // Test membership in the same physical edge. A single zero shape
        // function is sufficient for triangles, but not for quadrilaterals.
        bool common_edge = false;
        for (const auto& r_edge : r_geometry.GenerateEdges()) {
            const double ex = r_edge[1].X() - r_edge[0].X();
            const double ey = r_edge[1].Y() - r_edge[0].Y();
            const double edge_length = std::hypot(ex, ey);
            bool contains_segment = true;
            for (IndexType endpoint = 0; endpoint < 2; ++endpoint) {
                const double px = r_coordinates(endpoint, 0) - r_edge[0].X();
                const double py = r_coordinates(endpoint, 1) - r_edge[0].Y();
                const double projection = (px * ex + py * ey) / edge_length;
                contains_segment = contains_segment &&
                    std::abs(px * ey - py * ex) <= 1e-9 * edge_length &&
                    projection >= -1e-9 && projection <= edge_length + 1e-9;
            }
            common_edge = common_edge || contains_segment;
        }
        for (IndexType i = 0; i < parent_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            KRATOS_ERROR_IF(GetGeometry()[parent_nodes * side + i].Id() != r_node.Id())
                << "Condition nodes must follow the two parent node orders." << std::endl;
            KRATOS_ERROR_IF_NOT(node_ids.insert(r_node.Id()).second)
                << "Coupled parents must have independent nodes." << std::endl;
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION, r_node);
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X, r_node);
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y, r_node);
            KRATOS_CHECK_DOF_IN_NODE(PRESSURE, r_node);
        }
        KRATOS_ERROR_IF(!common_edge && !(side == 0 && UsesBackgroundExtension())) << "Interface segment must belong to a single edge of each parent." << std::endl;
    }
    return 0;
    KRATOS_CATCH("")
}

} // namespace Kratos
