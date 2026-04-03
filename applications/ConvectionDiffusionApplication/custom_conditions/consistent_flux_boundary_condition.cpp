// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    GitHub Copilot
//

#include "consistent_flux_boundary_condition.h"

#include "includes/checks.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/integration_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos
{

ConsistentFluxBoundaryCondition::ConsistentFluxBoundaryCondition(
    IndexType NewId,
    Geometry<Node>::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
}

ConsistentFluxBoundaryCondition::ConsistentFluxBoundaryCondition(
    IndexType NewId,
    Geometry<Node>::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

ConsistentFluxBoundaryCondition::~ConsistentFluxBoundaryCondition() = default;

Condition::Pointer ConsistentFluxBoundaryCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConsistentFluxBoundaryCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer ConsistentFluxBoundaryCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConsistentFluxBoundaryCondition>(NewId, pGeom, pProperties);
}

void ConsistentFluxBoundaryCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    CalculateConditionSystem(&rLeftHandSideMatrix, &rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("");
}

void ConsistentFluxBoundaryCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    CalculateConditionSystem(&rLeftHandSideMatrix, nullptr, rCurrentProcessInfo);

    KRATOS_CATCH("");
}

void ConsistentFluxBoundaryCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    CalculateConditionSystem(nullptr, &rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("");
}

void ConsistentFluxBoundaryCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const auto& r_process_info = rCurrentProcessInfo;
    const auto& r_settings = *r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_parent_geometry = GetParentElement().GetGeometry();
    const auto parent_num_nodes = r_parent_geometry.PointsNumber();

    if (rResult.size() != parent_num_nodes) {
        rResult.resize(parent_num_nodes, false);
    }

    for (IndexType i = 0; i < parent_num_nodes; ++i) {
        rResult[i] = r_parent_geometry[i].GetDof(r_unknown_var).EquationId();
    }

    KRATOS_CATCH("");
}

void ConsistentFluxBoundaryCondition::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const auto& r_process_info = rCurrentProcessInfo;
    const auto& r_settings = *r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_parent_geometry = GetParentElement().GetGeometry();
    const auto parent_num_nodes = r_parent_geometry.PointsNumber();

    if (rConditionalDofList.size() != parent_num_nodes) {
        rConditionalDofList.resize(parent_num_nodes);
    }

    for (IndexType i = 0; i < parent_num_nodes; ++i) {
        rConditionalDofList[i] = r_parent_geometry[i].pGetDof(r_unknown_var);
    }

    KRATOS_CATCH("");
}

GeometryData::IntegrationMethod ConsistentFluxBoundaryCondition::GetIntegrationMethod() const
{
    return IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(this->GetGeometry());
}

int ConsistentFluxBoundaryCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    int check = Condition::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS))
        << "No CONVECTION_DIFFUSION_SETTINGS defined in ProcessInfo." << std::endl;

    const auto& r_settings = *rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];

    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedUnknownVariable())
        << "No unknown variable defined in CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedDiffusionVariable())
        << "No diffusion variable defined in CONVECTION_DIFFUSION_SETTINGS." << std::endl;

    KRATOS_ERROR_IF_NOT(this->Has(NEIGHBOUR_ELEMENTS))
        << "NEIGHBOUR_ELEMENTS was not assigned to condition " << this->Id() << '.' << std::endl;

    const auto& r_neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_ERROR_IF(r_neighbours.size() != 1)
        << "Condition " << this->Id() << " requires exactly one parent element. Got " << r_neighbours.size() << '.' << std::endl;

    const auto& r_parent_geometry = r_neighbours[0].GetGeometry();
    KRATOS_ERROR_IF(r_parent_geometry.WorkingSpaceDimension() != r_parent_geometry.LocalSpaceDimension())
        << "Condition " << this->Id() << " requires a parent geometry with invertible Jacobian." << std::endl;

    const auto& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_diffusion_var = r_settings.GetDiffusionVariable();
    for (IndexType i = 0; i < r_parent_geometry.PointsNumber(); ++i) {
        const auto& r_node = r_parent_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_unknown_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_diffusion_var, r_node);
        KRATOS_CHECK_DOF_IN_NODE(r_unknown_var, r_node);
    }

    return check;

    KRATOS_CATCH("");
}

std::string ConsistentFluxBoundaryCondition::Info() const
{
    std::stringstream buffer;
    buffer << "ConsistentFluxBoundaryCondition #" << Id();
    return buffer.str();
}

void ConsistentFluxBoundaryCondition::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "ConsistentFluxBoundaryCondition #" << Id();
}

void ConsistentFluxBoundaryCondition::PrintData(std::ostream& rOStream) const
{
    rOStream << "ConsistentFluxBoundaryCondition #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}

Element& ConsistentFluxBoundaryCondition::GetParentElement()
{
    KRATOS_ERROR_IF_NOT(this->Has(NEIGHBOUR_ELEMENTS))
        << "NEIGHBOUR_ELEMENTS was not assigned to condition " << this->Id() << '.' << std::endl;

    auto& r_neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_ERROR_IF(r_neighbours.size() != 1)
        << "Condition " << this->Id() << " requires exactly one parent element. Got " << r_neighbours.size() << '.' << std::endl;

    return r_neighbours[0];
}

const Element& ConsistentFluxBoundaryCondition::GetParentElement() const
{
    KRATOS_ERROR_IF_NOT(this->Has(NEIGHBOUR_ELEMENTS))
        << "NEIGHBOUR_ELEMENTS was not assigned to condition " << this->Id() << '.' << std::endl;

    const auto& r_neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_ERROR_IF(r_neighbours.size() != 1)
        << "Condition " << this->Id() << " requires exactly one parent element. Got " << r_neighbours.size() << '.' << std::endl;

    return r_neighbours[0];
}

void ConsistentFluxBoundaryCondition::CalculateConditionSystem(
    MatrixType* pLeftHandSideMatrix,
    VectorType* pRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_condition_geometry = this->GetGeometry();
    const auto& r_parent_geometry = GetParentElement().GetGeometry();
    const auto parent_num_nodes = r_parent_geometry.PointsNumber();

    if (pLeftHandSideMatrix != nullptr) {
        if (pLeftHandSideMatrix->size1() != parent_num_nodes || pLeftHandSideMatrix->size2() != parent_num_nodes) {
            pLeftHandSideMatrix->resize(parent_num_nodes, parent_num_nodes, false);
        }
        noalias(*pLeftHandSideMatrix) = ZeroMatrix(parent_num_nodes, parent_num_nodes);
    }

    if (pRightHandSideVector != nullptr) {
        if (pRightHandSideVector->size() != parent_num_nodes) {
            pRightHandSideVector->resize(parent_num_nodes, false);
        }
        noalias(*pRightHandSideVector) = ZeroVector(parent_num_nodes);
    }

    const auto& r_settings = *rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_diffusion_var = r_settings.GetDiffusionVariable();

    const auto integration_method = this->GetIntegrationMethod();
    const auto& r_integration_points = r_condition_geometry.IntegrationPoints(integration_method);

    Vector det_j = ZeroVector(r_integration_points.size());
    r_condition_geometry.DeterminantOfJacobian(det_j, integration_method);

    MatrixType lhs_contribution(parent_num_nodes, parent_num_nodes);
    noalias(lhs_contribution) = ZeroMatrix(parent_num_nodes, parent_num_nodes);
    IntegrationPointData data;
    for (IndexType g = 0; g < r_integration_points.size(); ++g) {
        ComputeIntegrationPointData(
            r_condition_geometry,
            r_parent_geometry,
            g,
            r_integration_points,
            det_j,
            r_diffusion_var,
            data);

        for (IndexType i = 0; i < parent_num_nodes; ++i) {
            for (IndexType j = 0; j < parent_num_nodes; ++j) {
                double flux_projection = 0.0;
                for (IndexType d = 0; d < r_parent_geometry.WorkingSpaceDimension(); ++d) {
                    flux_projection += data.UnitNormal[d] * data.ParentDNDX(j, d);
                }
                lhs_contribution(i, j) += data.Weight * data.Conductivity * data.ParentN[i] * flux_projection;
            }
        }
    }

    if (pLeftHandSideMatrix != nullptr) {
        noalias(*pLeftHandSideMatrix) = lhs_contribution;
    }

    if (pRightHandSideVector != nullptr) {
        const auto parent_unknown_values = GetParentUnknownValues(r_unknown_var);
        noalias(*pRightHandSideVector) -= prod(lhs_contribution, parent_unknown_values);
    }
}

void ConsistentFluxBoundaryCondition::ComputeIntegrationPointData(
    const GeometryType& rConditionGeometry,
    const GeometryType& rParentGeometry,
    const IndexType IntegrationPointIndex,
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const Vector& rDetJ,
    const Variable<double>& rDiffusionVariable,
    IntegrationPointData& rData) const
{
    rData.Weight = rDetJ[IntegrationPointIndex] * rIntegrationPoints[IntegrationPointIndex].Weight();
    rData.UnitNormal = rConditionGeometry.UnitNormal(rIntegrationPoints[IntegrationPointIndex].Coordinates());

    CalculateParentGeometryValues(
        rConditionGeometry,
        rParentGeometry,
        rIntegrationPoints[IntegrationPointIndex].Coordinates(),
        rData.ParentN,
        rData.ParentDNDX);

    rData.Conductivity = 0.0;
    for (IndexType i = 0; i < rParentGeometry.PointsNumber(); ++i) {
        rData.Conductivity += rData.ParentN[i] * rParentGeometry[i].FastGetSolutionStepValue(rDiffusionVariable);
    }
}

void ConsistentFluxBoundaryCondition::CalculateParentGeometryValues(
    const GeometryType& rConditionGeometry,
    const GeometryType& rParentGeometry,
    const GeometryType::CoordinatesArrayType& rConditionLocalCoordinates,
    Vector& rParentN,
    Matrix& rParentDNDX) const
{
    Point global_coordinates;
    rConditionGeometry.GlobalCoordinates(global_coordinates, rConditionLocalCoordinates);

    Point parent_local_coordinates;
    rParentGeometry.PointLocalCoordinates(parent_local_coordinates, global_coordinates);

    rParentN.resize(rParentGeometry.PointsNumber(), false);
    rParentGeometry.ShapeFunctionsValues(rParentN, parent_local_coordinates);

    Matrix local_gradients;
    rParentGeometry.ShapeFunctionsLocalGradients(local_gradients, parent_local_coordinates);

    Matrix jacobian;
    rParentGeometry.Jacobian(jacobian, parent_local_coordinates);

    Matrix jacobian_inverse;
    double det_j = 0.0;
    MathUtils<double>::InvertMatrix(jacobian, jacobian_inverse, det_j);

    if (rParentDNDX.size1() != local_gradients.size1() || rParentDNDX.size2() != jacobian_inverse.size2()) {
        rParentDNDX.resize(local_gradients.size1(), jacobian_inverse.size2(), false);
    }
    noalias(rParentDNDX) = prod(local_gradients, jacobian_inverse);
}

Vector ConsistentFluxBoundaryCondition::GetParentUnknownValues(const Variable<double>& rUnknownVariable) const
{
    const auto& r_parent_geometry = GetParentElement().GetGeometry();
    Vector values(r_parent_geometry.PointsNumber());
    for (IndexType i = 0; i < r_parent_geometry.PointsNumber(); ++i) {
        values[i] = r_parent_geometry[i].FastGetSolutionStepValue(rUnknownVariable);
    }
    return values;
}

ConsistentFluxBoundaryCondition::ConsistentFluxBoundaryCondition()
    : Condition()
{
}

void ConsistentFluxBoundaryCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

void ConsistentFluxBoundaryCondition::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

} // namespace Kratos