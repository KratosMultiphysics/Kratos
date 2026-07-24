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
            "echo_level"       : 0
        })" );

    ThisParameters.ValidateAndAssignDefaults(default_parameters);
}

} // helpers namespace


void InitialFlatteningUtility::Execute(ModelPart& rModelPart, Parameters ThisParameters)
{
    ValidateParameters(ThisParameters);
    const int echo_level = ThisParameters["echo_level"].GetInt();

    const std::string& r_projection_type = ThisParameters["projection_type"].GetString();

    Vector3 normal;
    if (r_projection_type == "planar_fixed_direction") {
        normal = GetFixedDirection(ThisParameters);
    } else if (r_projection_type == "planar_mean_normal") {
        normal = ComputeMeanSurfaceNormal(rModelPart);
    } else {
        KRATOS_ERROR << "projection type: " << r_projection_type
            << " not available, please use planar_fixed_direction, planar_mean_normal" << std::endl;
    }

    ProjectNodesOntoPlane(rModelPart, normal, echo_level);
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

} // namespace Kratos.
