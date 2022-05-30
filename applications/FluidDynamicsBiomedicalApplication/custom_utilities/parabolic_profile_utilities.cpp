//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduardo Soudah
//                   Ruben Zorilla
//

// System includes

// External includes

// Project includes
#include "processes/parallel_distance_calculation_process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "fluid_dynamics_biomedical_application_variables.h"
#include "parabolic_profile_utilities.h"

namespace Kratos {

double ParabolicProfileUtilities::CalculateInletArea(const ModelPart& rModelPart)
{
    double inlet_area = 0.0;
    if (rModelPart.GetCommunicator().LocalMesh().NumberOfConditions() != 0) {
        inlet_area = block_for_each<SumReduction<double>>(rModelPart.Conditions(), [](const Condition& rCondition){
            return rCondition.GetGeometry().DomainSize();
        });
    }
    inlet_area = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(inlet_area);
    return inlet_area;
}

void ParabolicProfileUtilities::CalculateWallParallelDistance(
    ModelPart& rWallModelPart,
    ModelPart& rFluidModelPart,
    const std::size_t WallDistanceLevels)
{
    // Initialize the distance field to the geometry bounding box characteristic length
    // Also initialize the BOUNDARY flag that will be used to mark the wall
    const double char_length = CalculateBoundingBoxCharacteristicLength(rFluidModelPart);
    block_for_each(rFluidModelPart.Nodes(), [char_length](NodeType& rNode){
        rNode.Set(BOUNDARY, false);
        rNode.SetValue(WALL_DISTANCE, char_length);
    });

    // Set a small negative distance value in all wall nodes
    const double wall_d = -char_length * 1.0e-6;
    KRATOS_ERROR_IF(std::abs(wall_d) < std::numeric_limits<double>::epsilon())
        << "Auxiliary negative wall distance is close to zero. Check the domain dimensions." << std::endl;
    block_for_each(rWallModelPart.Nodes(), [wall_d](NodeType& rNode){
        rNode.Set(BOUNDARY, true);
        rNode.GetValue(WALL_DISTANCE) = wall_d;
    });
    rWallModelPart.GetCommunicator().SynchronizeOrNodalFlags(BOUNDARY);
    rWallModelPart.GetCommunicator().SynchronizeNonHistoricalDataToMin(WALL_DISTANCE);

    // Calculate the distance in the first layer of nodes
    // Note that these will be the values to be preserved by the parallel distance calculation below
    auto aux_tls = std::make_pair(array_1d<double,3>(), array_1d<double,3>());
    KRATOS_ERROR_IF(rWallModelPart.GetCommunicator().GlobalNumberOfConditions() == 0) << "Wall model part has no conditions." << std::endl;
    block_for_each(rWallModelPart.Conditions(), aux_tls, [&rWallModelPart](
        Condition& rCondition,
        std::pair<array_1d<double,3>, array_1d<double,3>>& rTLS)
    {
        // Get TLS
        auto& r_normal = std::get<0>(rTLS);
        auto& r_wall_vect = std::get<1>(rTLS);

        // Get condition geometry and calculate normal
        auto& r_geom = rCondition.GetGeometry();
        r_normal = r_geom.UnitNormal(0,GeometryData::IntegrationMethod::GI_GAUSS_1);

        // Get parent geometry to search for the first layer of nodes
        const auto& r_neighbours = rCondition.GetValue(NEIGHBOUR_ELEMENTS);
        KRATOS_ERROR_IF(r_neighbours.size() == 0) << "No parent found in condition with id " << rCondition.Id() << " in " <<  rWallModelPart.FullName() << "." << std::endl;
        auto& r_parent_geom = r_neighbours[0].GetGeometry();

        // Calculate the minimum distance from the wall to the parent nodes
        for (auto& r_parent_node : r_parent_geom) {
            if (!r_parent_node.Is(BOUNDARY)) {
                r_wall_vect = r_geom.Center() - r_parent_node.Coordinates();
                const double distance = std::abs(inner_prod(r_wall_vect, r_normal));
                r_parent_node.SetLock();
                double& r_parent_node_d = r_parent_node.GetValue(WALL_DISTANCE);
                r_parent_node_d = std::min(distance, r_parent_node_d);
                r_parent_node.UnSetLock();
            }
        }
    });
    rFluidModelPart.GetCommunicator().SynchronizeNonHistoricalDataToMin(WALL_DISTANCE);

    // Set the parallel distance settings
    // Note that the non-historical database and WALL_DISTANCE variables are used
    Parameters parallel_distance_settings(R"({
        "max_levels" : 1,
        "max_distance" : 1.0,
        "distance_variable" : "WALL_DISTANCE",
        "distance_database" : "nodal_non_historical"
    })");
    parallel_distance_settings["max_distance"].SetDouble(char_length);
    const std::size_t domain_size = rFluidModelPart.GetProcessInfo()[DOMAIN_SIZE];

    // Calculate parallel distance
    // Note that we check that the initial maximum levels are sufficient
    double current_max_dist = char_length;
    std::size_t current_max_levels = WallDistanceLevels / 2;
    while (std::abs(current_max_dist - char_length) < 1.0e-12) {
        // Calculate the parallel distance with the new number of levels
        current_max_levels *= 2;
        parallel_distance_settings["max_levels"].SetInt(current_max_levels);
        if (domain_size == 2) {
            ParallelDistanceCalculationProcess<2>(rFluidModelPart, parallel_distance_settings).Execute();
        } else if (domain_size == 3) {
            ParallelDistanceCalculationProcess<3>(rFluidModelPart, parallel_distance_settings).Execute();
        } else {
            KRATOS_ERROR << "Unsupported domain size in " << rFluidModelPart.Name() << ". Found DOMAIN_SIZE " << domain_size << " in ProcessInfo container." << std::endl;
        }

        // Get the new maximum distance
        current_max_dist = block_for_each<MaxReduction<double>>(rFluidModelPart.Nodes(), [](const NodeType& rNode){
            return rNode.GetValue(WALL_DISTANCE);
        });
        rFluidModelPart.GetCommunicator().GetDataCommunicator().MaxAll(current_max_dist);
    }

    // Reset boundary values to zero
    block_for_each(rWallModelPart.Nodes(), [](NodeType& rNode){
        rNode.GetValue(WALL_DISTANCE) = 0.0;
    });
}

void ParabolicProfileUtilities::ImposeParabolicInlet(
    ModelPart& rModelPart,
    const double MaxParabolaValue,
    const double MaxValueFactor)
{
    // Impose the parabolic inlet profile
    ImposeParabolicProfile(rModelPart, MaxParabolaValue, MaxValueFactor);
}

void ParabolicProfileUtilities::ImposeParabolicInlet(
    ModelPart& rModelPart,
    const GenericFunctionUtility::Pointer MaxParabolaValue,
    const double MaxValueFactor)
{
    // Impose the parabolic inlet profile
    ImposeParabolicProfile(rModelPart, MaxParabolaValue, MaxValueFactor);
}

template<class TInputType>
void ParabolicProfileUtilities::ImposeParabolicProfile(
    ModelPart &rModelPart,
    const TInputType& rMaxParabolaValue,
    const double MaxValueFactor)
{
    // Get time value from model part ProcessInfo
    const auto& r_process_info = rModelPart.GetProcessInfo();
    KRATOS_ERROR_IF_NOT(r_process_info.Has(TIME)) << "TIME is not present in '" << rModelPart.FullName() << "' ProcessInfo." << std::endl;
    const double time = r_process_info[TIME];

    // Get the maximum distance to the wall
    const double max_dist = block_for_each<MaxReduction<double>>(rModelPart.Nodes(), [](NodeType& rNode){
        return rNode.GetValue(WALL_DISTANCE);
    });
    KRATOS_ERROR_IF(std::abs(max_dist) < 1.0e-12) << "WALL_DISTANCE is close to zero." << std::endl;

    // Impose the parabolic profile values
    block_for_each(rModelPart.Nodes(), array_1d<double,3>(), [time, max_dist, MaxValueFactor, &rMaxParabolaValue](NodeType& rNode, array_1d<double,3>& rTLSValue){
        //Calculate distance for each node
        const double wall_dist = rNode.GetValue(WALL_DISTANCE) < 0.0 ? 0.0 : rNode.GetValue(WALL_DISTANCE);

        // Calculate the inlet direction
        rTLSValue = rNode.GetValue(INLET_NORMAL);
        const double n_norm = norm_2(rTLSValue);
        if (n_norm > 1.0e-12) {
            rTLSValue /= -n_norm;
        } else {
            KRATOS_WARNING("ImposeParabolicProfile") << "Node " << rNode.Id() << " INLET_NORMAL is close to zero." << std::endl;
            rTLSValue /= -1.0;
        }

        // Calculate the inlet value module
        const double max_value = MaxValueFactor * GetMaxParabolaValue(time, rNode, rMaxParabolaValue);
        const double value_in = max_value * (1.0-(std::pow(max_dist-wall_dist,2)/std::pow(max_dist,2)));
        rTLSValue *= value_in;

        // Set and fix the VELOCITY DOFs
        rNode.FastGetSolutionStepValue(VELOCITY) = rTLSValue;
        rNode.Fix(VELOCITY_X);
        rNode.Fix(VELOCITY_Y);
        rNode.Fix(VELOCITY_Z);
    });
}

void ParabolicProfileUtilities::FreeParabolicInlet(ModelPart& rModelPart)
{
    block_for_each(rModelPart.Nodes(), [](NodeType& rNode){
        rNode.Free(VELOCITY_X);
        rNode.Free(VELOCITY_Y);
        rNode.Free(VELOCITY_Z);
    });
}

template<>
double ParabolicProfileUtilities::GetMaxParabolaValue<double>(
    const double Time,
    const NodeType& rNode,
    const double& rMaxParabolaValue)
{
    return rMaxParabolaValue;
}

template<>
double ParabolicProfileUtilities::GetMaxParabolaValue<GenericFunctionUtility::Pointer>(
    const double Time,
    const NodeType& rNode,
    const GenericFunctionUtility::Pointer& rMaxParabolaValue)
{
    if(!rMaxParabolaValue->UseLocalSystem()) {
        return rMaxParabolaValue->CallFunction(rNode.X(), rNode.Y(), rNode.Z(), Time, rNode.X0(), rNode.Y0(), rNode.Z0());
    } else {
        return rMaxParabolaValue->RotateAndCallFunction(rNode.X(), rNode.Y(), rNode.Z(), Time, rNode.X0(), rNode.Y0(), rNode.Z0());
    }
}

double ParabolicProfileUtilities::CalculateBoundingBoxCharacteristicLength(const ModelPart& rModelPart)
{
    // Get the background mesh model part
    KRATOS_ERROR_IF(rModelPart.GetCommunicator().GlobalNumberOfNodes() == 0) << "Model part has no nodes." << std::endl;

    // Compute the domain characteristic length
    typedef CombinedReduction<MaxReduction<double>,MaxReduction<double>,MaxReduction<double>,MinReduction<double>,MinReduction<double>,MinReduction<double>> CustomReduction;
    double max_x, max_y, max_z, min_x, min_y, min_z;
    std::tie(max_x,max_y,max_z,min_x,min_y,min_z) = block_for_each<CustomReduction>(rModelPart.Nodes(),[](const NodeType& rNode){
        return std::make_tuple(rNode[0],rNode[1],rNode[2],rNode[0],rNode[1],rNode[2]);}
    );
    auto max_vector = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(std::vector<double>{max_x, max_y, max_z});
    auto min_vector = rModelPart.GetCommunicator().GetDataCommunicator().MinAll(std::vector<double>{min_x, min_y, min_z});

    const double char_length = std::sqrt(std::pow(max_vector[0] - min_vector[0], 2) + std::pow(max_vector[1] - min_vector[1], 2) + std::pow(max_vector[2] - min_vector[2], 2));
    KRATOS_ERROR_IF(char_length < std::numeric_limits<double>::epsilon()) << "Domain characteristic length is close to zero." << std::endl;

    return char_length;
}

template KRATOS_API(FLUID_DYNAMICS_BIOMEDICAL_APPLICATION) void ParabolicProfileUtilities::ImposeParabolicProfile<double>(ModelPart&, const double&, const double);
template KRATOS_API(FLUID_DYNAMICS_BIOMEDICAL_APPLICATION) void ParabolicProfileUtilities::ImposeParabolicProfile<GenericFunctionUtility::Pointer>(ModelPart&, const GenericFunctionUtility::Pointer&, const double);

} //namespace Kratos
