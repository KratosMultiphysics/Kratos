//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ricky Aristio
//

// System includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>

// External includes

// Project includes
#include "custom_utilities/iga_sbm_mls_extension_operator_utility.h"
#include "custom_conditions/coupling_sbm_extension_operator_6p_condition.h"
#include "utilities/mls_shape_functions_utility.h"
#include "iga_application_variables.h"

namespace Kratos
{

std::vector<IgaSbmMlsExtensionOperatorUtility::IndexType> IgaSbmMlsExtensionOperatorUtility::FindNClosestPoints(
    const std::vector<array_1d<double, 3>>& rCandidatePoints,
    const array_1d<double, 3>& rEvalPoint,
    const SizeType NClosest)
{
    KRATOS_TRY

    const SizeType n_candidates = rCandidatePoints.size();
    const SizeType n_closest = std::min(NClosest, n_candidates);

    std::vector<IndexType> indices(n_candidates);
    std::iota(indices.begin(), indices.end(), 0);

    std::vector<double> distances(n_candidates);
    for (IndexType i = 0; i < n_candidates; ++i) {
        distances[i] = norm_2(rCandidatePoints[i] - rEvalPoint);
    }

    std::partial_sort(
        indices.begin(), indices.begin() + n_closest, indices.end(),
        [&distances](IndexType a, IndexType b) { return distances[a] < distances[b]; });

    indices.resize(n_closest);
    return indices;

    KRATOS_CATCH("")
}

void IgaSbmMlsExtensionOperatorUtility::ComputeExtensionOperator(
    const Matrix& rCloudPointsParametric,
    const array_1d<double, 3>& rEvalPointParametric,
    const SizeType MLSOrder,
    Vector& rN,
    Matrix& rDN_DLocal)
{
    KRATOS_TRY

    const SizeType n_cloud = rCloudPointsParametric.size1();
    KRATOS_ERROR_IF(n_cloud == 0) << "IgaSbmMlsExtensionOperatorUtility: empty cloud." << std::endl;

    // Kernel radius = 1.3 * MAX distance in the cloud 
    double max_dist = 0.0;
    for (IndexType i = 0; i < n_cloud; ++i) {
        array_1d<double, 3> cloud_point = ZeroVector(3);
        cloud_point[0] = rCloudPointsParametric(i, 0);
        cloud_point[1] = rCloudPointsParametric(i, 1);
        cloud_point[2] = rCloudPointsParametric(i, 2);
        const double dist = norm_2(cloud_point - rEvalPointParametric);
        max_dist = std::max(max_dist, dist);
    }
    KRATOS_ERROR_IF(max_dist < 1e-14)
        << "IgaSbmMlsExtensionOperatorUtility: degenerate cloud, all points coincide with the evaluation point." << std::endl;
    const double h = 1.3 * max_dist;

    switch (MLSOrder) {
        case 1:
            MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<2, 1>(
                rCloudPointsParametric, rEvalPointParametric, h, rN, rDN_DLocal);
            break;
        case 2:
            MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<2, 2>(
                rCloudPointsParametric, rEvalPointParametric, h, rN, rDN_DLocal);
            break;
        default:
            KRATOS_ERROR << "IgaSbmMlsExtensionOperatorUtility: unsupported MLS order "
                         << MLSOrder << ". Only 1 and 2 are supported." << std::endl;
    }

    KRATOS_CATCH("")
}

std::vector<Element::Pointer> IgaSbmMlsExtensionOperatorUtility::ClassifyActiveElements(
    const std::vector<Element::Pointer>& rDomainElements,
    const array_1d<double, 3>& rEvalPointParametric,
    const SizeType KBuffer)
{
    KRATOS_TRY

    const SizeType n_elements = rDomainElements.size();
    std::vector<Element::Pointer> elements;
    elements.reserve(n_elements);
    std::vector<double> distances;
    distances.reserve(n_elements);

    for (const Element::Pointer& p_elem : rDomainElements) {
        const auto& r_geometry = p_elem->GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();
        const SizeType n_gp = r_integration_points.size();
        KRATOS_ERROR_IF(n_gp == 0) << "IgaSbmMlsExtensionOperatorUtility: element "
            << p_elem->Id() << " has no integration points." << std::endl;

        array_1d<double, 3> center = ZeroVector(3);
        for (IndexType gp = 0; gp < n_gp; ++gp) {
            center[0] += r_integration_points[gp][0];
            center[1] += r_integration_points[gp][1];
        }
        center /= static_cast<double>(n_gp);

        elements.push_back(p_elem);
        distances.push_back(norm_2(center - rEvalPointParametric));
    }

    const SizeType k_buffer = std::min(KBuffer, n_elements);
    std::vector<IndexType> order(n_elements);
    std::iota(order.begin(), order.end(), 0);
    std::partial_sort(
        order.begin(), order.begin() + k_buffer, order.end(),
        [&distances](IndexType a, IndexType b) { return distances[a] < distances[b]; });

    std::vector<bool> excluded(n_elements, false);
    for (IndexType i = 0; i < k_buffer; ++i) {
        excluded[order[i]] = true;
    }

    std::vector<Element::Pointer> active_elements;
    active_elements.reserve(n_elements - k_buffer);
    for (IndexType i = 0; i < n_elements; ++i) {
        if (!excluded[i]) {
            active_elements.push_back(elements[i]);
        }
    }
    return active_elements;

    KRATOS_CATCH("")
}

IgaSbmMlsExtensionOperatorUtility::ExtensionOperatorResult IgaSbmMlsExtensionOperatorUtility::ComputeAssembledExtensionOperator(
    const std::vector<Element::Pointer>& rActiveElements,
    const array_1d<double, 3>& rEvalPointParametric,
    const SizeType NClosest,
    const SizeType MLSOrder)
{
    KRATOS_TRY

    // Build the full cloud of quadrature points from the active elements
    std::vector<array_1d<double, 3>> cloud_positions;
    std::vector<std::vector<IndexType>> cloud_node_ids;
    std::vector<Vector> cloud_N_values;

    for (const auto& p_elem : rActiveElements) {
        const auto& r_geometry = p_elem->GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();
        const Matrix& r_N_all = r_geometry.ShapeFunctionsValues();
        const SizeType n_gp = r_integration_points.size();
        const SizeType n_cp = r_geometry.PointsNumber();

        std::vector<IndexType> node_ids(n_cp);
        for (IndexType i = 0; i < n_cp; ++i) {
            node_ids[i] = r_geometry[i].Id();
        }

        for (IndexType gp = 0; gp < n_gp; ++gp) {
            array_1d<double, 3> pos = ZeroVector(3);
            pos[0] = r_integration_points[gp][0];
            pos[1] = r_integration_points[gp][1];
            cloud_positions.push_back(pos);
            cloud_node_ids.push_back(node_ids);
            cloud_N_values.push_back(row(r_N_all, gp));
        }
    }

    KRATOS_ERROR_IF(cloud_positions.empty())
        << "IgaSbmMlsExtensionOperatorUtility: no active elements to build a cloud from." << std::endl;

    // Adaptive cloud growth
    constexpr double pou_tolerance = 1e-6;
    constexpr double max_weight_tolerance = 3.0;
    constexpr double linear_repro_tolerance = 1e-6;
    const SizeType n_available = cloud_positions.size();

    std::vector<IndexType> closest;
    Matrix cloud_coords_selected;
    Vector rN;
    Matrix rDN_DLocal;

    SizeType n_closest_try = std::min(NClosest, n_available);
    while (true) {
        closest = FindNClosestPoints(cloud_positions, rEvalPointParametric, n_closest_try);
        const SizeType n_closest = closest.size();

        cloud_coords_selected.resize(n_closest, 3, false);
        for (IndexType row_idx = 0; row_idx < n_closest; ++row_idx) {
            const auto& pos = cloud_positions[closest[row_idx]];
            cloud_coords_selected(row_idx, 0) = pos[0];
            cloud_coords_selected(row_idx, 1) = pos[1];
            cloud_coords_selected(row_idx, 2) = pos[2];
        }

        ComputeExtensionOperator(cloud_coords_selected, rEvalPointParametric, MLSOrder, rN, rDN_DLocal);

        double pou = 0.0;
        double max_abs_w = 0.0;
        double repro_xi = 0.0, repro_eta = 0.0;
        double max_dist = 0.0;
        for (IndexType row_idx = 0; row_idx < n_closest; ++row_idx) {
            pou += rN[row_idx];
            max_abs_w = std::max(max_abs_w, std::abs(rN[row_idx]));
            repro_xi += rN[row_idx] * cloud_coords_selected(row_idx, 0);
            repro_eta += rN[row_idx] * cloud_coords_selected(row_idx, 1);
            const double dx = cloud_coords_selected(row_idx, 0) - rEvalPointParametric[0];
            const double dy = cloud_coords_selected(row_idx, 1) - rEvalPointParametric[1];
            max_dist = std::max(max_dist, std::sqrt(dx * dx + dy * dy));
        }
        const double linear_scale = std::max(max_dist, 1e-12);
        const double linear_err_xi = std::abs(repro_xi - rEvalPointParametric[0]) / linear_scale;
        const double linear_err_eta = std::abs(repro_eta - rEvalPointParametric[1]) / linear_scale;

        const bool pou_ok = std::abs(pou - 1.0) <= pou_tolerance;
        const bool weight_ok = max_abs_w <= max_weight_tolerance;
        const bool linear_ok = (linear_err_xi <= linear_repro_tolerance) && (linear_err_eta <= linear_repro_tolerance);

        if (pou_ok && weight_ok && linear_ok) {
            if (n_closest_try > NClosest) {
                KRATOS_INFO("IgaSbmMlsExtensionOperatorUtility")
                    << "DIAG: cloud grown at eval=(" << rEvalPointParametric[0] << "," << rEvalPointParametric[1]
                    << ") from " << NClosest << " to " << n_closest_try << " (final POU=" << pou
                    << ", max|w|=" << max_abs_w << ", linear_err=(" << linear_err_xi << "," << linear_err_eta
                    << "))" << std::endl;
            }
            break;
        }
        if (n_closest_try >= n_available) {
            KRATOS_ERROR
                << "IgaSbmMlsExtensionOperatorUtility: MLS cloud degenerate even using all "
                << n_available << " available points (partition of unity = " << pou
                << ", expected 1.0; max|weight| = " << max_abs_w << ", expected <= "
                << max_weight_tolerance << "; linear reproduction error = (" << linear_err_xi << ","
                << linear_err_eta << "), expected <= " << linear_repro_tolerance << ") at MLSOrder "
                << MLSOrder << ". The point cloud does "
                << "not have enough spread in both directions to support this MLS order -- "
                << "try a lower MLSOrder or provide more/less graded active elements." << std::endl;
        }
        n_closest_try = std::min(n_closest_try * 2, n_available);
    }

    const SizeType n_closest = closest.size();

    // Assemble per-control-point weights
    std::unordered_map<IndexType, double> weight_map;
    std::unordered_map<IndexType, array_1d<double, 2>> gradient_map;

    for (IndexType row_idx = 0; row_idx < n_closest; ++row_idx) {
        const IndexType cloud_idx = closest[row_idx];
        const auto& node_ids = cloud_node_ids[cloud_idx];
        const Vector& N_values = cloud_N_values[cloud_idx];
        const double mls_weight = rN[row_idx];
        const double mls_grad_xi = rDN_DLocal(row_idx, 0);
        const double mls_grad_eta = rDN_DLocal(row_idx, 1);

        for (IndexType j = 0; j < node_ids.size(); ++j) {
            const IndexType gid = node_ids[j];
            const double Nj = N_values[j];

            weight_map[gid] += mls_weight * Nj;

            if (gradient_map.find(gid) == gradient_map.end()) {
                gradient_map[gid] = ZeroVector(2);
            }
            gradient_map[gid][0] += mls_grad_xi * Nj;
            gradient_map[gid][1] += mls_grad_eta * Nj;
        }
    }

    ExtensionOperatorResult result;
    result.NodeIds.reserve(weight_map.size());
    result.Weights.resize(weight_map.size());
    result.GradientWeights.resize(weight_map.size(), 2);

    IndexType row_idx = 0;
    for (const auto& entry : weight_map) {
        result.NodeIds.push_back(entry.first);
        result.Weights[row_idx] = entry.second;
        const auto& grad = gradient_map[entry.first];
        result.GradientWeights(row_idx, 0) = grad[0];
        result.GradientWeights(row_idx, 1) = grad[1];
        ++row_idx;
    }

    return result;

    KRATOS_CATCH("")
}

namespace {

// Builds an id -> Node::Pointer lookup covering every node in rDomainElements'
// own geometries
std::unordered_map<std::size_t, Node::Pointer> BuildNodeIdToPointerMap(
    const std::vector<Element::Pointer>& rDomainElements)
{
    std::unordered_map<std::size_t, Node::Pointer> id_to_node;
    for (const auto& p_elem : rDomainElements) {
        auto& r_geometry = p_elem->GetGeometry();
        for (std::size_t i = 0; i < r_geometry.PointsNumber(); ++i) {
            id_to_node.emplace(r_geometry[i].Id(), r_geometry.pGetPoint(i));
        }
    }
    return id_to_node;
}

} 

namespace {

std::vector<Element::Pointer> FilterElementsByExcludedNodes(
    const std::vector<Element::Pointer>& rElements,
    const std::vector<IgaSbmMlsExtensionOperatorUtility::IndexType>& rExcludedNodeIds)
{
    if (rExcludedNodeIds.empty()) {
        return rElements;
    }
    std::unordered_set<std::size_t> excluded(rExcludedNodeIds.begin(), rExcludedNodeIds.end());
    std::vector<Element::Pointer> filtered;
    filtered.reserve(rElements.size());
    for (const auto& p_elem : rElements) {
        const auto& r_geometry = p_elem->GetGeometry();
        bool touches_excluded = false;
        for (std::size_t i = 0; i < r_geometry.PointsNumber(); ++i) {
            if (excluded.count(r_geometry[i].Id()) > 0) {
                touches_excluded = true;
                break;
            }
        }
        if (!touches_excluded) {
            filtered.push_back(p_elem);
        }
    }
    return filtered;
}

}

void IgaSbmMlsExtensionOperatorUtility::PrecomputeAndStoreCouplingExtensionOperators(
    ModelPart& rCouplingModelPart,
    const std::vector<Element::Pointer>& rMasterDomainElements,
    const std::vector<Element::Pointer>& rSlaveDomainElements,
    const SizeType KBuffer,
    const SizeType NClosest,
    const SizeType MLSOrder,
    const std::vector<IndexType>& rExcludedNodeIds)
{
    KRATOS_TRY

    const std::vector<Element::Pointer> filtered_master_elements =
        FilterElementsByExcludedNodes(rMasterDomainElements, rExcludedNodeIds);
    const std::vector<Element::Pointer> filtered_slave_elements =
        FilterElementsByExcludedNodes(rSlaveDomainElements, rExcludedNodeIds);
    KRATOS_ERROR_IF(filtered_master_elements.empty() || filtered_slave_elements.empty())
        << "IgaSbmMlsExtensionOperatorUtility::PrecomputeAndStoreCouplingExtensionOperators: "
        << "rExcludedNodeIds filtered out ALL master or slave domain elements -- "
        << "the exclusion set is too aggressive relative to K_BUFFER/mesh density." << std::endl;

    for (auto& r_condition : rCouplingModelPart.Conditions()) {
        const auto& r_geometry_master = r_condition.GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = r_condition.GetGeometry().GetGeometryPart(1);

        const auto& r_integration_points_master = r_geometry_master.IntegrationPoints();
        const auto& r_integration_points_slave = r_geometry_slave.IntegrationPoints();
        const SizeType n_gp = r_integration_points_master.size();

        KRATOS_ERROR_IF(r_integration_points_slave.size() != n_gp)
            << "IgaSbmMlsExtensionOperatorUtility::PrecomputeAndStoreCouplingExtensionOperators: "
            << "condition " << r_condition.Id() << " has mismatched master/slave quadrature point counts." << std::endl;

        const Matrix& r_N_master = r_geometry_master.ShapeFunctionsValues();
        const Matrix& r_N_slave = r_geometry_slave.ShapeFunctionsValues();
        const SizeType n_cp_master = r_geometry_master.PointsNumber();
        const SizeType n_cp_slave = r_geometry_slave.PointsNumber();

        std::vector<ExtensionOperatorResult> master_results(n_gp);
        std::vector<ExtensionOperatorResult> slave_results(n_gp);

        for (IndexType gp = 0; gp < n_gp; ++gp) {
            array_1d<double, 3> eval_point_master = ZeroVector(3);
            for (IndexType i = 0; i < n_cp_master; ++i) {
                eval_point_master[0] += r_N_master(gp, i) * r_geometry_master[i].X0();
                eval_point_master[1] += r_N_master(gp, i) * r_geometry_master[i].Y0();
            }

            array_1d<double, 3> eval_point_slave = ZeroVector(3);
            for (IndexType i = 0; i < n_cp_slave; ++i) {
                eval_point_slave[0] += r_N_slave(gp, i) * r_geometry_slave[i].X0();
                eval_point_slave[1] += r_N_slave(gp, i) * r_geometry_slave[i].Y0();
            }

            const auto active_master = ClassifyActiveElements(filtered_master_elements, eval_point_master, KBuffer);
            const auto active_slave = ClassifyActiveElements(filtered_slave_elements, eval_point_slave, KBuffer);

            master_results[gp] = ComputeAssembledExtensionOperator(active_master, eval_point_master, NClosest, MLSOrder);
            slave_results[gp] = ComputeAssembledExtensionOperator(active_slave, eval_point_slave, NClosest, MLSOrder);
        }

        std::vector<Node::Pointer> all_dof_nodes;
        all_dof_nodes.reserve(n_cp_master + n_cp_slave);
        std::unordered_map<IndexType, IndexType> id_to_full_index;

        for (IndexType i = 0; i < n_cp_master; ++i) {
            const auto p_node = r_geometry_master.pGetPoint(i);
            id_to_full_index.emplace(p_node->Id(), all_dof_nodes.size());
            all_dof_nodes.push_back(p_node);
        }
        for (IndexType i = 0; i < n_cp_slave; ++i) {
            const auto p_node = r_geometry_slave.pGetPoint(i);
            if (id_to_full_index.find(p_node->Id()) == id_to_full_index.end()) {
                id_to_full_index.emplace(p_node->Id(), all_dof_nodes.size());
                all_dof_nodes.push_back(p_node);
            }
        }

        const auto master_id_to_node = BuildNodeIdToPointerMap(rMasterDomainElements);
        const auto slave_id_to_node = BuildNodeIdToPointerMap(rSlaveDomainElements);

        for (const auto& r_result : master_results) {
            for (const IndexType id : r_result.NodeIds) {
                if (id_to_full_index.find(id) == id_to_full_index.end()) {
                    id_to_full_index.emplace(id, all_dof_nodes.size());
                    all_dof_nodes.push_back(master_id_to_node.at(id));
                }
            }
        }
        for (const auto& r_result : slave_results) {
            for (const IndexType id : r_result.NodeIds) {
                if (id_to_full_index.find(id) == id_to_full_index.end()) {
                    id_to_full_index.emplace(id, all_dof_nodes.size());
                    all_dof_nodes.push_back(slave_id_to_node.at(id));
                }
            }
        }

        const SizeType n_dof_full = all_dof_nodes.size();
        Matrix master_weights_full = ZeroMatrix(n_gp, n_dof_full);
        Matrix slave_weights_full = ZeroMatrix(n_gp, n_dof_full);
        Matrix master_grad_xi_full = ZeroMatrix(n_gp, n_dof_full);
        Matrix master_grad_eta_full = ZeroMatrix(n_gp, n_dof_full);
        Matrix slave_grad_xi_full = ZeroMatrix(n_gp, n_dof_full);
        Matrix slave_grad_eta_full = ZeroMatrix(n_gp, n_dof_full);

        for (IndexType gp = 0; gp < n_gp; ++gp) {
            const auto& r_result_master = master_results[gp];
            for (IndexType i = 0; i < r_result_master.NodeIds.size(); ++i) {
                const IndexType col = id_to_full_index.at(r_result_master.NodeIds[i]);
                master_weights_full(gp, col) = r_result_master.Weights[i];
                master_grad_xi_full(gp, col) = r_result_master.GradientWeights(i, 0);
                master_grad_eta_full(gp, col) = r_result_master.GradientWeights(i, 1);
            }
            const auto& r_result_slave = slave_results[gp];
            for (IndexType i = 0; i < r_result_slave.NodeIds.size(); ++i) {
                const IndexType col = id_to_full_index.at(r_result_slave.NodeIds[i]);
                slave_weights_full(gp, col) = r_result_slave.Weights[i];
                slave_grad_xi_full(gp, col) = r_result_slave.GradientWeights(i, 0);
                slave_grad_eta_full(gp, col) = r_result_slave.GradientWeights(i, 1);
            }
        }

        r_condition.SetValue(SBM_MLS_ALL_DOF_NODES, all_dof_nodes);
        r_condition.SetValue(SBM_MLS_MASTER_WEIGHTS, master_weights_full);
        r_condition.SetValue(SBM_MLS_SLAVE_WEIGHTS, slave_weights_full);
        r_condition.SetValue(SBM_MLS_MASTER_GRADIENT_WEIGHTS_XI, master_grad_xi_full);
        r_condition.SetValue(SBM_MLS_MASTER_GRADIENT_WEIGHTS_ETA, master_grad_eta_full);
        r_condition.SetValue(SBM_MLS_SLAVE_GRADIENT_WEIGHTS_XI, slave_grad_xi_full);
        r_condition.SetValue(SBM_MLS_SLAVE_GRADIENT_WEIGHTS_ETA, slave_grad_eta_full);
    }

    KRATOS_CATCH("")
}

std::vector<IgaSbmMlsExtensionOperatorUtility::IndexType> IgaSbmMlsExtensionOperatorUtility::GetStoredDofNodeIds(
    const Condition& rCondition)
{
    KRATOS_TRY

    const auto& r_nodes = rCondition.GetValue(SBM_MLS_ALL_DOF_NODES);
    std::vector<IndexType> ids;
    ids.reserve(r_nodes.size());
    for (const auto& p_node : r_nodes) {
        ids.push_back(p_node->Id());
    }
    return ids;

    KRATOS_CATCH("")
}

namespace {

// Shape-function-weighted position of rGeometry's own single
// evaluation point 
array_1d<double, 3> PhysicalCenterOfGeometry(const Geometry<Node>& rGeometry)
{
    const Matrix& r_N = rGeometry.ShapeFunctionsValues();
    const std::size_t n = rGeometry.PointsNumber();
    array_1d<double, 3> center = ZeroVector(3);
    for (std::size_t i = 0; i < n; ++i) {
        center[0] += r_N(0, i) * rGeometry[i].X0();
        center[1] += r_N(0, i) * rGeometry[i].Y0();
        center[2] += r_N(0, i) * rGeometry[i].Z0();
    }
    return center;
}

// Single closest condition among rCandidates.
Condition::Pointer FindClosestConditionByGeometryImpl(
    const std::vector<Condition::Pointer>& rCandidates,
    const array_1d<double, 3>& rEvalPoint)
{
    Condition::Pointer best = nullptr;
    double best_dist = std::numeric_limits<double>::max();
    for (const auto& p_condition : rCandidates) {
        const double dist = norm_2(PhysicalCenterOfGeometry(p_condition->GetGeometry()) - rEvalPoint);
        if (dist < best_dist) {
            best_dist = dist;
            best = p_condition;
        }
    }
    KRATOS_ERROR_IF(best == nullptr)
        << "IgaSbmMlsExtensionOperatorUtility: no candidate conditions to search." << std::endl;
    return best;
}

} 

void IgaSbmMlsExtensionOperatorUtility::PrecomputeAndStoreCouplingShiftSources(
    ModelPart& rCouplingModelPart,
    const std::vector<Condition::Pointer>& rMasterSurrogateConditions,
    const std::vector<Condition::Pointer>& rSlaveSurrogateConditions)
{
    KRATOS_TRY

    for (auto& r_condition : rCouplingModelPart.Conditions()) {
        auto p_typed = dynamic_cast<CouplingSbmExtensionOperator6pCondition*>(&r_condition);
        KRATOS_ERROR_IF(p_typed == nullptr)
            << "IgaSbmMlsExtensionOperatorUtility::PrecomputeAndStoreCouplingShiftSources: condition "
            << r_condition.Id() << " is not a CouplingSbmExtensionOperator6pCondition." << std::endl;

        const array_1d<double, 3> eval_point_master =
            PhysicalCenterOfGeometry(r_condition.GetGeometry().GetGeometryPart(0));
        const array_1d<double, 3> eval_point_slave =
            PhysicalCenterOfGeometry(r_condition.GetGeometry().GetGeometryPart(1));

        Condition::Pointer p_master_source = FindClosestConditionByGeometryImpl(rMasterSurrogateConditions, eval_point_master);
        Condition::Pointer p_slave_source = FindClosestConditionByGeometryImpl(rSlaveSurrogateConditions, eval_point_slave);

        auto all_dof_nodes = r_condition.GetValue(SBM_MLS_ALL_DOF_NODES);
        std::unordered_map<std::size_t, std::size_t> id_to_index;
        id_to_index.reserve(all_dof_nodes.size());
        for (std::size_t k = 0; k < all_dof_nodes.size(); ++k) {
            id_to_index.emplace(all_dof_nodes[k]->Id(), k);
        }

        const std::size_t n_dof_full_old = all_dof_nodes.size();
        for (const Condition::Pointer& p_source : {p_master_source, p_slave_source}) {
            auto& r_source_geometry = p_source->GetGeometry();
            for (std::size_t i = 0; i < r_source_geometry.PointsNumber(); ++i) {
                const auto p_node = r_source_geometry.pGetPoint(i);
                if (id_to_index.find(p_node->Id()) == id_to_index.end()) {
                    id_to_index.emplace(p_node->Id(), all_dof_nodes.size());
                    all_dof_nodes.push_back(p_node);
                }
            }
        }

        if (all_dof_nodes.size() != n_dof_full_old) {
            const std::size_t n_dof_full_new = all_dof_nodes.size();
            for (const Variable<Matrix>* p_var : {
                     &SBM_MLS_MASTER_WEIGHTS, &SBM_MLS_SLAVE_WEIGHTS,
                     &SBM_MLS_MASTER_GRADIENT_WEIGHTS_XI, &SBM_MLS_MASTER_GRADIENT_WEIGHTS_ETA,
                     &SBM_MLS_SLAVE_GRADIENT_WEIGHTS_XI, &SBM_MLS_SLAVE_GRADIENT_WEIGHTS_ETA}) {
                const Matrix& r_old = r_condition.GetValue(*p_var);
                Matrix grown = ZeroMatrix(r_old.size1(), n_dof_full_new);
                for (std::size_t r = 0; r < r_old.size1(); ++r) {
                    for (std::size_t c = 0; c < r_old.size2(); ++c) {
                        grown(r, c) = r_old(r, c);
                    }
                }
                r_condition.SetValue(*p_var, grown);
            }
            r_condition.SetValue(SBM_MLS_ALL_DOF_NODES, all_dof_nodes);
        }

        p_typed->SetShiftSources(p_master_source, p_slave_source);
    }

    KRATOS_CATCH("")
}

double IgaSbmMlsExtensionOperatorUtility::GetDifferentialArea(const Element& rElement)
{
    KRATOS_TRY

    const auto& r_geometry = rElement.GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints();

    KRATOS_ERROR_IF(r_integration_points.size() != 1)
        << "IgaSbmMlsExtensionOperatorUtility::GetDifferentialArea: element " << rElement.Id()
        << " has " << r_integration_points.size() << " integration points, expected exactly 1." << std::endl;

    Vector det_jacobian;
    r_geometry.DeterminantOfJacobian(det_jacobian);

    return r_integration_points[0].Weight() * det_jacobian[0];

    KRATOS_CATCH("")
}

}  // namespace Kratos.
