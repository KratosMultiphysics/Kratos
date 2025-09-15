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
#include "utilities/variable_utils.h"
#include "../FluidDynamicsApplication/fluid_dynamics_application_variables.h"

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

ModelPart& ParabolicProfileUtilities::CreateAndFillInletAuxiliaryVolumeModelPart(ModelPart& rInletModelPart)
{
    // Create an auxiliary model part to store the elements attached to the inlet face
    // Note that we do this at the level of the parent in order to avoid adding elements to the current inlet model part
    auto& r_inlet_aux_model_part = rInletModelPart.GetParentModelPart().CreateSubModelPart(rInletModelPart.Name()+"_auxiliary");

    // Flag the first layer of elements attached to the inlet face
    // For this we flag the inlet nodes and then check the elements sharing these nodes
    std::set<std::size_t> nodes_ids;
    std::vector<std::size_t> elem_ids;
    VariableUtils().SetFlag(SELECTED, true, rInletModelPart.Nodes());
    for (auto& r_element : rInletModelPart.GetRootModelPart().Elements()) {
        bool is_selected = false;
        const auto& r_geom = r_element.GetGeometry();
        // Check if current element is attached to the inlet by checking if it owns an inlet node
        for (const auto& r_node : r_geom) {
            if (r_node.Is(SELECTED)) {
                is_selected = true;
                elem_ids.push_back(r_element.Id());
                break;
            }
        }
        // If current element is an inlet element, add also its nodes (required by the parallel distance)
        if (is_selected) {
            for (const auto& r_node : r_geom) {
                nodes_ids.insert(r_node.Id());
            }
        }
    }

    // Add nodes and elements to the auxiliary model part
    r_inlet_aux_model_part.AddNodes(std::vector<std::size_t>(nodes_ids.begin(), nodes_ids.end()));
    r_inlet_aux_model_part.AddElements(elem_ids);

    // Return the new model part
    return r_inlet_aux_model_part;
}

void ParabolicProfileUtilities::CalculateWallParallelDistance(
    ModelPart& rWallModelPart,
    ModelPart& rFluidModelPart,
    const std::size_t WallDistanceLevels)
{
    // Initialize the distance field to the geometry bounding box characteristic length
    // Also initialize the BOUNDARY flag that will be used to mark the wall
    // Note that we do this with the parent model part in case only the slice of elements attached to the inlet is passed
    const auto& r_parent_model_part = rFluidModelPart.GetParentModelPart();
    const double char_length = CalculateBoundingBoxCharacteristicLength(r_parent_model_part);
    block_for_each(r_parent_model_part.Nodes(), [char_length](NodeType& rNode){
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
    std::size_t dist_it = 0;
    const std::size_t max_dist_its = 10;
    double current_max_dist = char_length;
    std::size_t current_max_levels = WallDistanceLevels / 2;
    while (std::abs(current_max_dist - char_length) < 1.0e-12 && dist_it < max_dist_its) {
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

        // Update the parallel distance iterations counter
        ++dist_it;
        KRATOS_WARNING_IF("CalculateWallParallelDistance", dist_it == max_dist_its) << "Reached maximum wall distance calculation iterations. Check input geometry." << std::endl;
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
    ImposeParabolicProfile(rModelPart, *MaxParabolaValue, MaxValueFactor);
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

    // TLS with the Inlet normal and the rMaxParabolaValue (to prevent race conditions if its a function)
    struct TLSType{
        array_1d<double,3> mInlet;
        TInputType mParabolaValue;

        TLSType(const TInputType & parabolaValue) : mInlet(array_1d<double,3>()), mParabolaValue(parabolaValue) {};
    };

    TLSType tls(rMaxParabolaValue);

    // Impose the parabolic profile values
    block_for_each(rModelPart.Nodes(), tls, [time, max_dist, MaxValueFactor](NodeType& rNode, TLSType& rTLSValue){
        //Calculate distance for each node
        const double wall_dist = rNode.GetValue(WALL_DISTANCE) < 0.0 ? 0.0 : rNode.GetValue(WALL_DISTANCE);

        // Calculate the inlet direction
        rTLSValue.mInlet = rNode.GetValue(INLET_NORMAL);
        const double n_norm = norm_2(rTLSValue.mInlet);
        if (n_norm > 1.0e-12) {
            rTLSValue.mInlet /= -n_norm;
        } else {
            KRATOS_WARNING("ImposeParabolicProfile") << "Node " << rNode.Id() << " INLET_NORMAL is close to zero." << std::endl;
            rTLSValue.mInlet /= -1.0;
        }

        // Calculate the inlet value module
        const double max_value = MaxValueFactor * GetMaxParabolaValue(time, rNode, rTLSValue.mParabolaValue);
        const double value_in = max_value * (1.0-(std::pow(max_dist-wall_dist,2)/std::pow(max_dist,2)));
        rTLSValue.mInlet *= value_in;

        // Set and fix the VELOCITY DOFs
        rNode.FastGetSolutionStepValue(VELOCITY) = rTLSValue.mInlet;
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
    double& rMaxParabolaValue)
{
    return rMaxParabolaValue;
}

template<>
double ParabolicProfileUtilities::GetMaxParabolaValue<GenericFunctionUtility>(
    const double Time,
    const NodeType& rNode,
    GenericFunctionUtility& rMaxParabolaValue)
{
    if(!rMaxParabolaValue.UseLocalSystem()) {
        return rMaxParabolaValue.CallFunction(rNode.X(), rNode.Y(), rNode.Z(), Time, rNode.X0(), rNode.Y0(), rNode.Z0());
    } else {
        return rMaxParabolaValue.RotateAndCallFunction(rNode.X(), rNode.Y(), rNode.Z(), Time, rNode.X0(), rNode.Y0(), rNode.Z0());
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
template KRATOS_API(FLUID_DYNAMICS_BIOMEDICAL_APPLICATION) void ParabolicProfileUtilities::ImposeParabolicProfile<GenericFunctionUtility>(ModelPart&, const GenericFunctionUtility&, const double);

} //namespace Kratos
