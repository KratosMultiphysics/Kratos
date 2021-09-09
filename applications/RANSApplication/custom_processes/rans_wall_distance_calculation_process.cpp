//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <limits>

// External includes

// Project includes
#include "includes/communicator.h"
#include "includes/variables.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/parallel_levelset_distance_calculator.h"
#include "utilities/variable_utils.h"

// Application includes

// Include base h
#include "rans_wall_distance_calculation_process.h"

namespace Kratos
{
namespace WallDistanceCalculationUtilities
{
/**
 * @brief Calculates wall distance and updates nodal value
 *
 * This method calculates wall distance. Wall distance is the projection of distance vector
 * from node to wall point in the direction of normal
 *
 * @param rNode                 Node, where wall distance is required
 * @param rWallLocation         Nearest wall location
 * @param rUnitNormal           Wall normal (this should be always outward pointing normal)
 * @param rDistanceVariable     Distance variable to store nodal distance
 */
void CalculateAndUpdateNodalMinimumWallDistance(
    ModelPart::NodeType& rNode,
    const array_1d<double, 3>& rWallLocation,
    const array_1d<double, 3>& rUnitNormal,
    const Variable<double>& rDistanceVariable)
{
    // Even if the rUnitNormal is outward pointing, there are situations (such as trailing edge point)
    // where the nodal averaged unit normal will be outward pointing, but the nodal location of interest
    // may give negative wall distances. So std::abs is used in here to avoid that.
    const double wall_distance =
        std::abs(inner_prod(rWallLocation - rNode.Coordinates(), rUnitNormal));

    rNode.SetLock();
    double& current_distance = rNode.FastGetSolutionStepValue(rDistanceVariable);
    current_distance = std::min(current_distance, wall_distance);
    rNode.Set(VISITED, true);
    rNode.UnSetLock();
}
} // namespace WallDistanceCalculationUtilities

RansWallDistanceCalculationProcess::RansWallDistanceCalculationProcess(
    Model& rModel,
    Parameters rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mMaxLevels = rParameters["max_levels"].GetInt();
    mEchoLevel = rParameters["echo_level"].GetInt();
    mMaxDistance = rParameters["max_distance"].GetDouble();
    mMainModelPartName = rParameters["main_model_part_name"].GetString();
    mWallModelPartName = rParameters["wall_model_part_name"].GetString();
    mDistanceVariableName = rParameters["distance_variable_name"].GetString();
    mNodalAreaVariableName = rParameters["nodal_area_variable_name"].GetString();
    mRecalculateAtEachTimeStep = rParameters["re_calculate_at_each_time_step"].GetBool();

    KRATOS_CATCH("");
}

int RansWallDistanceCalculationProcess::Check()
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mMainModelPartName);

    const auto& r_distance_variable =
        KratosComponents<Variable<double>>::Get(mDistanceVariableName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(r_distance_variable))
        << r_distance_variable.Name() << " is not found in nodal solution step variables list of "
        << mMainModelPartName << " used to store nodal distances.";

    const auto& r_nodal_area_variable = KratosComponents<Variable<double>>::Get(mNodalAreaVariableName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(r_nodal_area_variable))
        << r_nodal_area_variable.Name() << " is not found in nodal solution step variables list of "
        << mMainModelPartName << " used to store nodal areas.";

    return RansFormulationProcess::Check();

    KRATOS_CATCH("");
}

void RansWallDistanceCalculationProcess::ExecuteInitialize()
{
    CalculateWallDistances();
}

void RansWallDistanceCalculationProcess::ExecuteInitializeSolutionStep()
{
    if (mRecalculateAtEachTimeStep) {
        CalculateWallDistances();
    }
}

std::string RansWallDistanceCalculationProcess::Info() const
{
    return std::string("RansWallDistanceCalculationProcess");
}

void RansWallDistanceCalculationProcess::CalculateWallDistances()
{
    KRATOS_TRY

    using namespace WallDistanceCalculationUtilities;

    auto& r_main_model_part = mrModel.GetModelPart(mMainModelPartName);
    auto& r_wall_model_part = mrModel.GetModelPart(mWallModelPartName);
    auto& r_communicator = r_main_model_part.GetCommunicator();

    const auto& r_distance_variable = KratosComponents<Variable<double>>::Get(mDistanceVariableName);
    const auto& r_nodal_area_variable = KratosComponents<Variable<double>>::Get(mNodalAreaVariableName);

    // variable initialization
    block_for_each(r_main_model_part.Nodes(), [&](ModelPart::NodeType& rNode){
        rNode.SetValue(NORMAL, NORMAL.Zero());
        rNode.Set(VISITED, false);
        rNode.FastGetSolutionStepValue(r_distance_variable) = mMaxDistance;
    });

    // identify all wall nodes to negative small distance values
    block_for_each(r_wall_model_part.Nodes(), [&](ModelPart::NodeType& rNode){
        rNode.FastGetSolutionStepValue(r_distance_variable) = -1e-12;
        rNode.Set(VISITED, true);
    });

    // calculate wall distances based on conditions
    block_for_each(r_wall_model_part.Conditions(), [&](ModelPart::ConditionType& rCondition) {
        array_1d<double, 3> normal = rCondition.GetValue(NORMAL);
        const double normal_magnitude = norm_2(normal);
        KRATOS_ERROR_IF(normal_magnitude == 0.0)
            << "NORMAL is not properly initialized in condition with id "
            << rCondition.Id() << " at " << rCondition.GetGeometry().Center() << ".\n";
        normal /= normal_magnitude;

        auto& parent_geometry = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0].GetGeometry();
        for (auto& r_parent_node : parent_geometry) {
            if (r_parent_node.FastGetSolutionStepValue(r_distance_variable) > 0.0) {
                CalculateAndUpdateNodalMinimumWallDistance(
                    r_parent_node, rCondition.GetGeometry().Center(), normal, r_distance_variable);
            }
        }

        // update nodal non-historical NORMAL to be used in case where
        // some nodes needs updating of nodal distances via Elements.
        auto& r_geometry = rCondition.GetGeometry();
        for (auto& r_node : r_geometry) {
            r_node.SetLock();
            r_node.GetValue(NORMAL) += normal;
            r_node.UnSetLock();
        }
    });

    // communication for all ranks
    r_communicator.AssembleNonHistoricalData(NORMAL);
    r_communicator.SynchronizeCurrentDataToMin(r_distance_variable);
    r_communicator.SynchronizeOrNodalFlags(VISITED);

    // calculate distances based on elements (on wall adjacent nodes which are not covered by conditions)
    block_for_each(r_main_model_part.Elements(), [&](ModelPart::ElementType& rElement) {
        auto& r_geometry = rElement.GetGeometry();

        array_1d<double, 3> normal = ZeroVector(3);
        array_1d<double, 3> normal_center = ZeroVector(3);
        double normals_count = 0.0;
        std::vector<int> nodal_indices_to_update;

        // identify nodes' distances which are not calculated by conditions, and are adjacent to wall nodes.
        // compute average normal and average wall location from nodal normals and nodal coordinates which
        // lies on the wall.
        for (std::size_t i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {
            const auto& r_node = r_geometry[i_node];
            if (r_node.FastGetSolutionStepValue(r_distance_variable) < 0.0) {
                noalias(normal) += r_node.GetValue(NORMAL);
                noalias(normal_center) += r_node.Coordinates();
                normals_count += 1.0;
            } else if (!r_node.Is(VISITED)) {
                nodal_indices_to_update.push_back(i_node);
            }
        }

        if (normals_count > 0.0 && nodal_indices_to_update.size() > 0) {
            const double normal_magnitude = norm_2(normal);
            KRATOS_ERROR_IF(normal_magnitude == 0.0)
                << "NORMAL is not properly initialized in adjacent wall nodes "
                   "in element with id "
                << rElement.Id() << " at " << rElement.GetGeometry().Center() << ".\n";

            normal /= normal_magnitude;
            normal_center /= normals_count;

            for (const auto i_node : nodal_indices_to_update) {
                auto& r_node = r_geometry[i_node];
                CalculateAndUpdateNodalMinimumWallDistance(
                    r_node, normal_center, normal, r_distance_variable);
            }
        }
    });

    // communication for all ranks
    r_communicator.SynchronizeCurrentDataToMin(r_distance_variable);
    r_communicator.SynchronizeOrNodalFlags(VISITED);

    // update rest of the domain
    VariableUtils().SetVariable(r_distance_variable, 0.0, r_main_model_part.Nodes(),
                                VISITED, false);

    const int domain_size = r_main_model_part.GetProcessInfo()[DOMAIN_SIZE];
    if (domain_size == 2) {
        ParallelDistanceCalculator<2>().CalculateDistances(
            r_main_model_part, r_distance_variable, r_nodal_area_variable, mMaxLevels, mMaxDistance);
    } else if (domain_size == 3) {
        ParallelDistanceCalculator<3>().CalculateDistances(
            r_main_model_part, r_distance_variable, r_nodal_area_variable, mMaxLevels, mMaxDistance);
    } else {
        KRATOS_ERROR << "Unsupported domain size in " << r_main_model_part.Name()
                     << ". [ DOMAIN_SIZE = " << domain_size << " ]\n";
    }

    // revert boundary negative distances to zero
    VariableUtils().SetVariable(r_distance_variable, 0.0, r_wall_model_part.Nodes());

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Wall distances calculated in " << mMainModelPartName << ".\n";

    KRATOS_CATCH("");
}

const Parameters RansWallDistanceCalculationProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "main_model_part_name"             : "PLEASE_SPECIFY_MAIN_MODEL_PART_NAME",
            "wall_model_part_name"             : "PLEASE_SPECIFY_WALL_MODEL_PART_NAME",
            "max_levels"                       : 100,
            "max_distance"                     : 1e+30,
            "echo_level"                       : 0,
            "distance_variable_name"           : "DISTANCE",
            "nodal_area_variable_name"         : "NODAL_AREA",
            "re_calculate_at_each_time_step"   : false
        })");

    return default_parameters;
}

} // namespace Kratos.
