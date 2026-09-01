// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ricky Aristio
//

// System includes
#include <unordered_map>
#include <map>
#include <utility>
#include <vector>
#include <algorithm>

// External includes

// Project includes
#include "initial_flattening_utility.h"

namespace Kratos {

namespace {

typedef std::size_t SizeType;
typedef array_1d<double, 3> Vector3;

void NormalizeVector(Vector3& rVector)
{
    const double vec_norm = norm_2(rVector);
    KRATOS_ERROR_IF(vec_norm < 1e-12) << "Vector has length of zero!" << std::endl;
    rVector /= vec_norm;
}

void ValidateParameters(Parameters ThisParameters)
{
    Parameters default_parameters = Parameters(R"(
        {
            "projection_type"  : "planar_mean_normal",
            "global_direction" : [0.0, 0.0, 1.0],
            "sub_model_parts"  : [],
            "length_fit_iterations"  : 1000,
            "length_fit_tolerance"   : 1e-10,
            "length_fit_relaxation"  : 0.8,
            "echo_level"       : 0
        })" );

    ThisParameters.ValidateAndAssignDefaults(default_parameters);
}

void CollectEdges(const ModelPart& rModelPart, std::vector<std::pair<SizeType, SizeType>>& rEdges, std::vector<double>& rTargetLengths)
{
    std::map<std::pair<SizeType, SizeType>, double> unique_edges;

    for (const auto& r_element : rModelPart.Elements()) {
        const auto& r_geom = r_element.GetGeometry();
        if (r_geom.LocalSpaceDimension() != 2) continue;

        const SizeType number_of_nodes = r_geom.size();
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node_a = r_geom[i];
            const auto& r_node_b = r_geom[(i + 1) % number_of_nodes];

            const SizeType id_a = r_node_a.Id();
            const SizeType id_b = r_node_b.Id();

            const Vector3& r_pos_a = r_node_a.GetInitialPosition();
            const Vector3& r_pos_b = r_node_b.GetInitialPosition();
            const Vector3 delta = r_pos_a - r_pos_b;

            unique_edges[(id_a < id_b) ? std::make_pair(id_a, id_b)
                                       : std::make_pair(id_b, id_a)] = norm_2(delta);
        }
    }

    KRATOS_ERROR_IF(unique_edges.empty())
        << "InitialFlatteningUtility: no 2D elements found to build edges from!" << std::endl;

    rEdges.clear();
    rTargetLengths.clear();
    rEdges.reserve(unique_edges.size());
    rTargetLengths.reserve(unique_edges.size());
    for (const auto& r_entry : unique_edges) {
        rEdges.push_back(r_entry.first);
        rTargetLengths.push_back(r_entry.second);
    }
}
} // helpers namespace


void InitialFlatteningUtility::Execute(ModelPart& rModelPart, Parameters ThisParameters)
{
    ValidateParameters(ThisParameters);
    const int echo_level = ThisParameters["echo_level"].GetInt();

    const std::string& r_projection_type = ThisParameters["projection_type"].GetString();
    const std::vector<std::string> sub_model_part_names = ThisParameters["sub_model_parts"].GetStringArray();

    if (r_projection_type == "planar_fixed_direction") {
        const Vector3 normal = GetFixedDirection(ThisParameters);
        ProjectNodesOntoPlane(rModelPart, normal, echo_level);
    } else if (r_projection_type == "planar_fixed_direction_lsq") {
        const Vector3 normal = GetFixedDirection(ThisParameters);
        ProjectNodesOntoPlane(rModelPart, normal, echo_level);
        FitEdgeLengths(rModelPart, normal,
                       ThisParameters["length_fit_iterations"].GetInt(),
                       ThisParameters["length_fit_tolerance"].GetDouble(),
                       ThisParameters["length_fit_relaxation"].GetDouble(),
                       echo_level);
    } else if (r_projection_type == "planar_mean_normal") {
        if (sub_model_part_names.empty()) {
            const Vector3 normal = ComputeMeanSurfaceNormal(rModelPart);
            ProjectNodesOntoPlane(rModelPart, normal, echo_level);
        } else {
            ProjectNodesOntoPlanePerPart(rModelPart, sub_model_part_names, echo_level);
        }
    } else {
        KRATOS_ERROR << "projection type: " << r_projection_type
            << " not available, please use planar_fixed_direction, planar_fixed_direction_lsq, planar_mean_normal" << std::endl;
    }
}

InitialFlatteningUtility::Vector3 InitialFlatteningUtility::GetFixedDirection(Parameters ThisParameters)
{
    const Vector direction_vec = ThisParameters["global_direction"].GetVector();
    KRATOS_ERROR_IF_NOT(direction_vec.size() == 3) << "\"global_direction\" must be of size 3!" << std::endl;

    Vector3 direction;
    direction[0] = direction_vec[0];
    direction[1] = direction_vec[1];
    direction[2] = direction_vec[2];

    NormalizeVector(direction);
    return direction;
}

InitialFlatteningUtility::Vector3 InitialFlatteningUtility::ComputeMeanSurfaceNormal(const ModelPart& rModelPart)
{
    std::unordered_map<SizeType, Vector3> nodal_directors;

    for (const auto& r_element : rModelPart.Elements()) {
        const auto& r_geom = r_element.GetGeometry();
        if (r_geom.LocalSpaceDimension() != 2) continue;

        const auto integration_method = r_geom.GetDefaultIntegrationMethod();
        const auto& r_integration_points = r_geom.IntegrationPoints(integration_method);

        Vector3 element_normal = ZeroVector(3);
        for (SizeType gp = 0; gp < r_integration_points.size(); ++gp) {
            element_normal += r_integration_points[gp].Weight() * r_geom.Normal(gp, integration_method);
        }
        NormalizeVector(element_normal);

        for (const auto& r_node : r_geom) {
            auto it = nodal_directors.find(r_node.Id());
            if (it == nodal_directors.end()) {
                nodal_directors.emplace(r_node.Id(), element_normal);
            } else {
                it->second += element_normal;
            }
        }
    }

    KRATOS_ERROR_IF(nodal_directors.empty()) << "InitialFlatteningUtility: no 2D elements found to compute a mean surface normal from!" << std::endl;

    Vector3 mean_normal = ZeroVector(3);
    for (auto& r_entry : nodal_directors) {
        NormalizeVector(r_entry.second);
        mean_normal += r_entry.second;
    }
    NormalizeVector(mean_normal);

    return mean_normal;
}

void InitialFlatteningUtility::ProjectNodesOntoPlane(ModelPart& rModelPart, const Vector3& rNormal, const int EchoLevel)
{
    // Projection of the reference position onto the plane through the origin with normal rNormal.
    const SizeType buffer_size = rModelPart.GetBufferSize();

    for (auto& r_node : rModelPart.Nodes()) {
        const Vector3& r_ref_position = r_node.GetInitialPosition();
        const Vector3 proj_vec = -inner_prod(rNormal, r_ref_position) * rNormal;

        for (SizeType i = 0; i < buffer_size; ++i) {
            r_node.FastGetSolutionStepValue(DISPLACEMENT, i) = proj_vec;
        }
    }

    KRATOS_INFO_IF("InitialFlatteningUtility", EchoLevel > 0)
        << "Initial flattening projected " << rModelPart.NumberOfNodes()
        << " nodes onto plane with normal " << rNormal << std::endl;
}

void InitialFlatteningUtility::ProjectNodesOntoPlanePerPart(ModelPart& rModelPart, const std::vector<std::string>& rSubModelPartNames, const int EchoLevel)
{
    // One mean normal per named sub model part 
    std::vector<Vector3> part_normals;
    part_normals.reserve(rSubModelPartNames.size());
    for (const auto& r_name : rSubModelPartNames) {
        const Vector3 part_normal = ComputeMeanSurfaceNormal(rModelPart.GetSubModelPart(r_name));
        part_normals.push_back(part_normal);
        KRATOS_INFO_IF("InitialFlatteningUtility", EchoLevel > 0)
            << "Mean surface normal for \"" << r_name << "\": " << part_normal << std::endl;
    }

    std::unordered_map<SizeType, Vector3> node_normal_sum;
    for (SizeType i = 0; i < rSubModelPartNames.size(); ++i) {
        for (const auto& r_node : rModelPart.GetSubModelPart(rSubModelPartNames[i]).Nodes()) {
            auto it = node_normal_sum.find(r_node.Id());
            if (it == node_normal_sum.end()) {
                node_normal_sum.emplace(r_node.Id(), part_normals[i]);
            } else {
                it->second += part_normals[i];
            }
        }
    }

    const SizeType buffer_size = rModelPart.GetBufferSize();
    SizeType num_projected = 0;
    for (auto& r_node : rModelPart.Nodes()) {
        auto it = node_normal_sum.find(r_node.Id());
        KRATOS_ERROR_IF(it == node_normal_sum.end())
            << "InitialFlatteningUtility: node " << r_node.Id()
            << " does not belong to any of the given \"sub_model_parts\"!" << std::endl;

        Vector3 normal = it->second;
        NormalizeVector(normal);

        const Vector3& r_ref_position = r_node.GetInitialPosition();
        const Vector3 proj_vec = -inner_prod(normal, r_ref_position) * normal;

        for (SizeType i = 0; i < buffer_size; ++i) {
            r_node.FastGetSolutionStepValue(DISPLACEMENT, i) = proj_vec;
        }
        ++num_projected;
    }

    KRATOS_INFO_IF("InitialFlatteningUtility", EchoLevel > 0)
        << "Initial flattening (per-part) projected " << num_projected << " nodes across "
        << rSubModelPartNames.size() << " sub model parts." << std::endl;
}

void InitialFlatteningUtility::FitEdgeLengths(
    ModelPart& rModelPart,
    const Vector3& rNormal,
    const int MaxIterations,
    const double Tolerance,
    const double RelaxationFactor,
    const int EchoLevel)
{
    std::vector<std::pair<SizeType, SizeType>> edges;
    std::vector<double> target_lengths;
    CollectEdges(rModelPart, edges, target_lengths);

    // current in-plane positions = X0 + u
    std::unordered_map<SizeType, Vector3> pos;
    for (auto& r_node : rModelPart.Nodes()) {
        const Vector3& r_ref_position = r_node.GetInitialPosition();
        pos[r_node.Id()] = r_ref_position + r_node.FastGetSolutionStepValue(DISPLACEMENT);
    }

    std::unordered_map<SizeType, Vector3> correction;
    std::unordered_map<SizeType, double> count;

    double max_step = 0.0;
    int it = 0;
    for (; it < MaxIterations; ++it) {
        for (auto& r_entry : pos) {
            correction[r_entry.first] = ZeroVector(3);
            count[r_entry.first] = 0.0;
        }

        for (SizeType m = 0; m < edges.size(); ++m) {
            const Vector3 d = pos[edges[m].first] - pos[edges[m].second];
            const double length = norm_2(d);
            if (length < 1e-14) continue;

            const Vector3 shift = (0.5 * (target_lengths[m] - length) / length) * d;
            correction[edges[m].first]  += shift;
            correction[edges[m].second] -= shift;
            count[edges[m].first]  += 1.0;
            count[edges[m].second] += 1.0;
        }

        max_step = 0.0;
        for (auto& r_entry : pos) {
            if (count[r_entry.first] < 1.0) continue;
            Vector3 step = (RelaxationFactor / count[r_entry.first]) * correction[r_entry.first];
            step -= inner_prod(step, rNormal) * rNormal;   // stay in the flattening plane
            r_entry.second += step;
            max_step = std::max(max_step, norm_2(step));
        }
        if (max_step < Tolerance) break;
    }

    const SizeType buffer_size = rModelPart.GetBufferSize();
    for (auto& r_node : rModelPart.Nodes()) {
        const Vector3& r_ref_position = r_node.GetInitialPosition();
        const Vector3 disp = pos[r_node.Id()] - r_ref_position;
        for (SizeType i = 0; i < buffer_size; ++i) {
            r_node.FastGetSolutionStepValue(DISPLACEMENT, i) = disp;
        }
    }

    double lambda_min = 1.0e30;
    double lambda_max = 0.0;
    for (SizeType m = 0; m < edges.size(); ++m) {
        const double lambda = norm_2(pos[edges[m].first] - pos[edges[m].second]) / target_lengths[m];
        lambda_min = std::min(lambda_min, lambda);
        lambda_max = std::max(lambda_max, lambda);
    }

    KRATOS_INFO_IF("InitialFlatteningUtility", EchoLevel > 0)
        << "Least-squares length fit: " << it << " iterations, max step " << max_step
        << ", edge stretch in [" << lambda_min << ", " << lambda_max << "]" << std::endl;

    KRATOS_WARNING_IF("InitialFlatteningUtility", lambda_min < 0.75)
        << "Initial pattern still has edges compressed to " << lambda_min
        << " of their reference length. Green-Lagrange strain is near the -0.5 floor "
        << "and the relaxation phase is likely to fold elements." << std::endl;
}

} // namespace Kratos.
