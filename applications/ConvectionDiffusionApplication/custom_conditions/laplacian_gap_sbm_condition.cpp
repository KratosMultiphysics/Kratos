// BSD License, Kratos default license: kratos/license.txt
#include <algorithm>
#include <cmath>
#include <limits>
#include "custom_conditions/laplacian_gap_sbm_condition.h"
#include "convection_diffusion_application_variables.h"
#include "includes/convection_diffusion_settings.h"
#include "includes/element.h"
#include "includes/checks.h"
#include "integration/line_gauss_legendre_integration_points.h"
#include "utilities/math_utils.h"

namespace Kratos {

Condition::Pointer LaplacianGapDirichletCondition2D::Create(IndexType Id, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const
{ return Kratos::make_intrusive<LaplacianGapDirichletCondition2D>(Id, pGeometry, pProperties); }
Condition::Pointer LaplacianGapDirichletCondition2D::Create(IndexType Id, const NodesArrayType& rNodes, PropertiesType::Pointer pProperties) const
{ return Create(Id, Kratos::make_shared<Geometry<Node>>(rNodes), pProperties); }
Condition::Pointer LaplacianGapConnectorCondition2D::Create(IndexType Id, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const
{ return Kratos::make_intrusive<LaplacianGapConnectorCondition2D>(Id, pGeometry, pProperties); }
Condition::Pointer LaplacianGapConnectorCondition2D::Create(IndexType Id, const NodesArrayType& rNodes, PropertiesType::Pointer pProperties) const
{ return Create(Id, Kratos::make_shared<Geometry<Node>>(rNodes), pProperties); }

void LaplacianGapDirichletCondition2D::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rValues, const ProcessInfo& rInfo)
{
    if (rVariable != INTEGRATION_COORDINATES) {
        Condition::CalculateOnIntegrationPoints(rVariable, rValues, rInfo);
        return;
    }
    const auto& r_segment = GetValue(IsConnector() ? SURROGATE_BOUNDARY_FACE_COORDINATES : GAP_INTERFACE_COORDINATES);
    const auto& r_rule = LineGaussLegendreIntegrationPoints3::IntegrationPoints();
    rValues.resize(r_rule.size());
    for (std::size_t g = 0; g < r_rule.size(); ++g) {
        const double xi = r_rule[g].X();
        for (std::size_t d = 0; d < 3; ++d)
            rValues[g][d] = 0.5 * ((1.0 - xi) * r_segment(0,d) + (1.0 + xi) * r_segment(1,d));
    }
}

void LaplacianGapDirichletCondition2D::GetDofList(DofsVectorType& rDofs, const ProcessInfo& rInfo) const
{
    const auto& r_unknown = rInfo[CONVECTION_DIFFUSION_SETTINGS]->GetUnknownVariable();
    rDofs.resize(GetGeometry().size());
    for (IndexType i = 0; i < GetGeometry().size(); ++i) rDofs[i] = GetGeometry()[i].pGetDof(r_unknown);
}

void LaplacianGapDirichletCondition2D::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rInfo) const
{
    const auto& r_unknown = rInfo[CONVECTION_DIFFUSION_SETTINGS]->GetUnknownVariable();
    rResult.resize(GetGeometry().size(), false);
    for (IndexType i = 0; i < GetGeometry().size(); ++i) rResult[i] = GetGeometry()[i].GetDof(r_unknown).EquationId();
}

int LaplacianGapDirichletCondition2D::Check(const ProcessInfo& rInfo) const
{
    KRATOS_ERROR_IF_NOT(rInfo.Has(CONVECTION_DIFFUSION_SETTINGS)) << "Missing convection-diffusion settings." << std::endl;
    KRATOS_ERROR_IF_NOT(Has(NEIGHBOUR_ELEMENTS)) << "Scalar gap condition requires parent elements." << std::endl;
    KRATOS_ERROR_IF(GetValue(NEIGHBOUR_ELEMENTS).size() != (IsConnector() ? 2 : 1)) << "Invalid scalar gap parents." << std::endl;
    KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_COEFFICIENT) && GetProperties()[PENALTY_COEFFICIENT] > 0.0)
        << "Scalar gap condition requires a positive penalty." << std::endl;
    const auto& r_unknown = rInfo[CONVECTION_DIFFUSION_SETTINGS]->GetUnknownVariable();
    for (const auto& r_node : GetGeometry()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_unknown, r_node);
        KRATOS_CHECK_DOF_IN_NODE(r_unknown, r_node);
    }
    return 0;
}

void LaplacianGapDirichletCondition2D::CalculateLocalSystem(MatrixType& rLHS, VectorType& rRHS, const ProcessInfo& rInfo)
{
    KRATOS_TRY
    const auto& r_parents = GetValue(NEIGHBOUR_ELEMENTS);
    const Matrix& r_segment = GetValue(IsConnector() ? SURROGATE_BOUNDARY_FACE_COORDINATES : GAP_INTERFACE_COORDINATES);
    const auto& r_rule = LineGaussLegendreIntegrationPoints3::IntegrationPoints();
    const auto& r_unknown = rInfo[CONVECTION_DIFFUSION_SETTINGS]->GetUnknownVariable();
    const auto& r_nodes = GetGeometry();
    const std::size_t size = r_nodes.size();
    rLHS = ZeroMatrix(size, size);
    rRHS = ZeroVector(size);
    const double tx = r_segment(1, 0) - r_segment(0, 0);
    const double ty = r_segment(1, 1) - r_segment(0, 1);
    const double length = std::hypot(tx, ty);
    KRATOS_ERROR_IF(length <= 0.0) << "Scalar gap segment has zero length." << std::endl;
    const double nx = ty / length, ny = -tx / length;
    double h = std::numeric_limits<double>::max();
    for (std::size_t s = 0; s < r_parents.size(); ++s) h = std::min(h, std::sqrt(r_parents[s].GetGeometry().Area()));
    KRATOS_ERROR_IF(h <= 0.0) << "Scalar gap parent has zero area." << std::endl;
    const double penalty = GetProperties()[PENALTY_COEFFICIENT] / h;
    if (!IsConnector()) {
        KRATOS_ERROR_IF_NOT(Has(GAP_DIRICHLET_VALUES) && GetValue(GAP_DIRICHLET_VALUES).size() == r_rule.size())
            << "Scalar gap Dirichlet condition needs one value per integration point." << std::endl;
    }
    for (std::size_t g = 0; g < r_rule.size(); ++g) {
        const double xi = r_rule[g].X();
        const Point point(0.5 * ((1.0-xi)*r_segment(0,0)+(1.0+xi)*r_segment(1,0)),
                          0.5 * ((1.0-xi)*r_segment(0,1)+(1.0+xi)*r_segment(1,1)), 0.0);
        Vector jump = ZeroVector(size), flux = ZeroVector(size);
        for (std::size_t s = 0; s < r_parents.size(); ++s) {
            const auto& r_parent = r_parents[s].GetGeometry();
            GeometryType::CoordinatesArrayType local = ZeroVector(3);
            r_parent.PointLocalCoordinates(local, point);
            Vector shape;
            Matrix local_gradients, jacobian, inverse;
            double determinant = 0.0;
            r_parent.ShapeFunctionsValues(shape, local);
            r_parent.ShapeFunctionsLocalGradients(local_gradients, local);
            r_parent.Jacobian(jacobian, local);
            MathUtils<double>::InvertMatrix(jacobian, inverse, determinant);
            KRATOS_ERROR_IF(std::abs(determinant) < 1e-20) << "Singular scalar gap parent." << std::endl;
            const Matrix gradients = prod(local_gradients, inverse);
            for (std::size_t i = 0; i < r_parent.size(); ++i) {
                std::size_t k = 0;
                while (k < size && r_nodes[k].Id() != r_parent[i].Id()) ++k;
                KRATOS_ERROR_IF(k == size) << "Scalar gap parent node is absent from condition geometry." << std::endl;
                jump[k] += (s == 0 ? 1.0 : -1.0) * shape[i];
                flux[k] += (IsConnector() ? 0.5 : 1.0) * (gradients(i,0)*nx + gradients(i,1)*ny);
            }
        }
        const double weight = 0.5 * length * r_rule[g].Weight();
        const double boundary_value = IsConnector() ? 0.0 : GetValue(GAP_DIRICHLET_VALUES)[g];
        for (std::size_t i = 0; i < size; ++i) {
            rRHS[i] += weight * (penalty * jump[i] - flux[i]) * boundary_value;
            for (std::size_t j = 0; j < size; ++j)
                rLHS(i,j) += weight * (-jump[i]*flux[j] - flux[i]*jump[j] + penalty*jump[i]*jump[j]);
        }
    }
    for (std::size_t i = 0; i < size; ++i)
        for (std::size_t j = 0; j < size; ++j)
            rRHS[i] -= rLHS(i,j) * r_nodes[j].FastGetSolutionStepValue(r_unknown);
    KRATOS_CATCH("")
}
}
