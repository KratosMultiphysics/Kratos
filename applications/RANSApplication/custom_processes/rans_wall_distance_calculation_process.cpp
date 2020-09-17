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
#include <unordered_map>
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
    // rUnitNormal is assumed to be outward pointing, hence wall_distance will be always positive.
    const double wall_distance =
        inner_prod(rWallLocation - rNode.Coordinates(), rUnitNormal);

    rNode.SetLock();
    double& current_distance = rNode.FastGetSolutionStepValue(rDistanceVariable);
    if (current_distance > wall_distance) {
        current_distance = wall_distance;
    }
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

    if (rParameters.Has("max_distance")) {
        if (rParameters["max_distance"].IsString()) {
            if (rParameters["max_distance"].GetString() == "max") {
                rParameters["max_distance"].SetDouble(std::numeric_limits<double>::max());
            }
            else {
                KRATOS_ERROR << "Unsupported \"max_distance\" value. Supports "
                                "only \"max\". [ max_distance = "
                             << rParameters["max_distance"].GetString() << " ]\n";
            }
        }
    }

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mMaxLevels = rParameters["max_levels"].GetInt();
    mEchoLevel = rParameters["echo_level"].GetInt();
    mMaxDistance = rParameters["max_distance"].GetDouble();
    mModelPartName = rParameters["model_part_name"].GetString();
    mWallFlagVariableName = rParameters["wall_flag_variable_name"].GetString();
    mWallFlagVariableValue = rParameters["wall_flag_variable_value"].GetBool();
    mDistanceVariableName = rParameters["distance_variable_name"].GetString();
    mNodalAreaVariableName = rParameters["nodal_area_variable_name"].GetString();
    mRecalculateAtEachTimeStep = rParameters["re_calculate_at_each_time_step"].GetBool();
    mCorrectDistancesUsingNeighbors =
        rParameters["correct_distances_using_neighbors"].GetBool();

    KRATOS_CATCH("");
}

int RansWallDistanceCalculationProcess::Check()
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    const auto& r_distance_variable =
        KratosComponents<Variable<double>>::Get(mDistanceVariableName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(r_distance_variable))
        << r_distance_variable.Name() << " is not found in nodal solution step variables list of "
        << mModelPartName << " used to store nodal distances.";

    const auto& r_nodal_area_variable = KratosComponents<Variable<double>>::Get(mNodalAreaVariableName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(r_nodal_area_variable))
        << r_nodal_area_variable.Name() << " is not found in nodal solution step variables list of "
        << mModelPartName << " used to store nodal areas.";

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

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    auto& r_communicator = r_model_part.GetCommunicator();

    const auto& r_wall_flag = KratosComponents<Flags>::Get(mWallFlagVariableName);
    const auto& r_distance_variable = KratosComponents<Variable<double>>::Get(mDistanceVariableName);
    const auto& r_nodal_area_variable = KratosComponents<Variable<double>>::Get(mNodalAreaVariableName);

    // variable initialization
    block_for_each(r_model_part.Nodes(), [&](ModelPart::NodeType& rNode){
        rNode.SetValue(NORMAL, NORMAL.Zero());
        rNode.Set(VISITED, false);
        rNode.FastGetSolutionStepValue(r_distance_variable) = mMaxDistance;
    });

    // calculate wall distances based on conditions
    block_for_each(r_model_part.Conditions(), [&](ModelPart::ConditionType& rCondition) {
        if (rCondition.Is(r_wall_flag) == mWallFlagVariableValue) {
            array_1d<double, 3> normal = rCondition.GetValue(NORMAL);
            const double normal_magnitude = norm_2(normal);
            KRATOS_ERROR_IF(normal_magnitude == 0.0)
                << "NORMAL is not properly initialied in condition with id "
                << rCondition.Id() << " at "
                << rCondition.GetGeometry().Center() << ".\n";
            normal /= normal_magnitude;

            auto& parent_geometry =
                rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0].GetGeometry();
            for (auto& r_parent_node : parent_geometry) {
                if (r_parent_node.Is(r_wall_flag) != mWallFlagVariableValue) {
                    CalculateAndUpdateNodalMinimumWallDistance(
                        r_parent_node, rCondition.GetGeometry().Center(),
                        normal, r_distance_variable);
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
        }
    });

    // communication for all ranks
    r_communicator.AssembleNonHistoricalData(NORMAL);
    r_communicator.SynchronizeCurrentDataToMin(r_distance_variable);
    r_communicator.SynchronizeOrNodalFlags(VISITED);

    // calculate distances based on elements (on wall adjacent nodes which are not covered by conditions)
    block_for_each(r_model_part.Elements(), [&](ModelPart::ElementType& rElement) {
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
            if (r_node.Is(r_wall_flag) == mWallFlagVariableValue) {
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
                << "NORMAL is not properly initialied in adjacent wall nodes "
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
    block_for_each(r_model_part.Nodes(), [&](ModelPart::NodeType& rNode) {
        if (rNode.Is(r_wall_flag) == mWallFlagVariableValue) {
            rNode.FastGetSolutionStepValue(r_distance_variable) = -1e-12;
        } else if (!rNode.Is(VISITED)) {
            rNode.FastGetSolutionStepValue(r_distance_variable) = 0.0;
        }
    });

    const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
    if (domain_size == 2) {
        auto p_distance_smoother = Kratos::make_shared<ParallelDistanceCalculator<2>>();
        p_distance_smoother->CalculateDistances(
            r_model_part, r_distance_variable, r_nodal_area_variable, mMaxLevels, mMaxDistance);
    } else if (domain_size == 3) {
        auto p_distance_smoother = Kratos::make_shared<ParallelDistanceCalculator<3>>();
        p_distance_smoother->CalculateDistances(
            r_model_part, r_distance_variable, r_nodal_area_variable, mMaxLevels, mMaxDistance);
    } else {
        KRATOS_ERROR << "Unsupported domain size in " << r_model_part.Name()
                     << ". [ DOMAIN_SIZE = " << domain_size << " ]\n";
    }

    // revert boundary negative distances to zero
    VariableUtils().SetVariable(r_distance_variable, 0.0, r_model_part.Nodes(),
                                r_wall_flag, mWallFlagVariableValue);

    if (mCorrectDistancesUsingNeighbors) {
        CorrectWallDistances();
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Wall distances calculated in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

void RansWallDistanceCalculationProcess::CorrectWallDistances()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    auto& r_communicator = r_model_part.GetCommunicator();
    auto& r_data_communicator = r_communicator.GetDataCommunicator();

    auto& r_nodes = r_communicator.LocalMesh().Nodes();

    const auto& r_distance_variable = KratosComponents<Variable<double>>::Get(mDistanceVariableName);

    VariableUtils().SetNonHistoricalVariableToZero(r_distance_variable, r_model_part.Nodes());

    FindGlobalNodalNeighboursProcess find_nodal_neighbours_process(
        r_data_communicator, r_model_part);
    find_nodal_neighbours_process.Execute();

    class GlobalPointerAdder
    {
    public:
        typedef GlobalPointersVector<ModelPart::NodeType> value_type;
        value_type gp_vector;

        value_type GetValue()
        {
            gp_vector.Unique();
            return gp_vector;
        }

        void LocalReduce(const value_type& rGPVector)
        {
            for (auto& r_gp : rGPVector.GetContainer()) {
                this->gp_vector.push_back(r_gp);
            }
        }
        void ThreadSafeReduce(GlobalPointerAdder& rOther)
        {
#pragma omp critical
            {
                for (auto& r_gp : rOther.gp_vector.GetContainer()) {
                    this->gp_vector.push_back(r_gp);
                }
            }
        }
    };

    auto all_global_pointers =
        block_for_each<GlobalPointerAdder>(r_nodes, [](ModelPart::NodeType& rNode) {
            return rNode.GetValue(NEIGHBOUR_NODES);
        });

    GlobalPointerCommunicator<ModelPart::NodeType> pointer_comm(
        r_data_communicator, all_global_pointers);

    auto distance_proxy =
        pointer_comm.Apply([&](const GlobalPointer<ModelPart::NodeType>& gp) -> double {
            return gp->FastGetSolutionStepValue(r_distance_variable);
        });

    int number_of_modified_nodes = block_for_each<SumReduction<int>>(
        r_nodes, [&](ModelPart::NodeType& rNode) -> int {
            if (rNode.FastGetSolutionStepValue(r_distance_variable) < 0.0) {
                const auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES);

                int count = 0;
                double average_value = 0.0;
                for (int j_node = 0;
                     j_node < static_cast<int>(r_neighbours.size()); ++j_node) {
                    const double j_distance = distance_proxy.Get(r_neighbours(j_node));
                    if (j_distance > 0.0) {
                        average_value += j_distance;
                        count++;
                    }
                }

                if (count > 0) {
                    rNode.SetValue(r_distance_variable,
                                   average_value / static_cast<double>(count));
                } else {
                    KRATOS_ERROR
                        << "Node " << rNode.Id() << " at " << rNode.Coordinates()
                        << " didn't find any neighbour(" << r_neighbours.size()
                        << ") with positive " << r_distance_variable.Name()
                        << ". Please recheck " << mModelPartName << ".\n";
                }

                return 1;
            }
            return 0;
        });

    block_for_each(r_nodes, [&](ModelPart::NodeType& rNode) {
        double& r_distance = rNode.FastGetSolutionStepValue(r_distance_variable);
        const double avg_distance = rNode.GetValue(r_distance_variable);
        if (r_distance < 0.0) {
            r_distance = avg_distance;
        }
    });

    r_communicator.SynchronizeVariable(r_distance_variable);
    number_of_modified_nodes =
        r_communicator.GetDataCommunicator().SumAll(number_of_modified_nodes);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0 && number_of_modified_nodes > 0)
        << "Corrected " << number_of_modified_nodes
        << " nodal wall distances in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

const Parameters RansWallDistanceCalculationProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name"                  : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "max_levels"                       : 100,
            "max_distance"                     : 1e+30,
            "echo_level"                       : 0,
            "distance_variable_name"           : "DISTANCE",
            "nodal_area_variable_name"         : "NODAL_AREA",
            "wall_flag_variable_name"          : "STRUCTURE",
            "wall_flag_variable_value"         : true,
            "re_calculate_at_each_time_step"   : false,
            "correct_distances_using_neighbors": true
        })");

    return default_parameters;
}

} // namespace Kratos.
