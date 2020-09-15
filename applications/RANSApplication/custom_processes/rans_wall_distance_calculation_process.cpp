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
void UpdateToMinimumDistances(
    std::vector<int>& rIndices,
    std::vector<double>& rDistances,
    const int NodeId,
    const double WallDistance)
{
    const auto itr = std::find(rIndices.begin(), rIndices.end(), NodeId);
    if (itr == rIndices.end()) {
        rIndices.push_back(NodeId);
        rDistances.push_back(WallDistance);
    } else {
        const auto index = std::distance(rIndices.begin(), itr);
        double& current_distance = rDistances[index];
        if (current_distance > WallDistance) {
            current_distance = WallDistance;
        }
    }
}

void UpdateModelPartDistancesToMinimum(
    ModelPart& rModelPart,
    std::vector<int>& rLocalIndices,
    std::vector<double>& rLocalDistances,
    const Variable<double>& rDistanceVariable)
{
    KRATOS_TRY

    const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();

    KRATOS_ERROR_IF(rLocalIndices.size() != rLocalDistances.size())
        << "Local indices list and distances list size mismatch. [ "
        << rLocalIndices.size() << " != " << rLocalDistances.size() << " ].\n";

    const auto& global_indices_vector = r_data_communicator.Gatherv(rLocalIndices, 0);
    const auto& global_distances_vector =
        r_data_communicator.Gatherv(rLocalDistances, 0);

    rLocalIndices.clear();
    rLocalDistances.clear();

    if (r_data_communicator.Rank() == 0) {
        for (std::size_t rank = 0; rank < global_indices_vector.size(); ++rank) {
            const auto& current_indices = global_indices_vector[rank];
            const auto& current_distances = global_distances_vector[rank];

            for (std::size_t i = 0; i < current_indices.size(); ++i) {
                UpdateToMinimumDistances(rLocalIndices, rLocalDistances,
                                         current_indices[i], current_distances[i]);
            }
        }
    }

    int number_of_indices = rLocalIndices.size();
    r_data_communicator.Broadcast(number_of_indices, 0);
    rLocalIndices.resize(number_of_indices);
    rLocalDistances.resize(number_of_indices);

    r_data_communicator.Broadcast(rLocalIndices, 0);
    r_data_communicator.Broadcast(rLocalDistances, 0);

    IndexPartition<int>(number_of_indices).for_each([&](const int i) {
        const int node_id = rLocalIndices[i];

        if (rModelPart.GetMesh().HasNode(node_id)) {
            const double distance = rLocalDistances[i];
            auto& r_node = rModelPart.GetNode(node_id);
            r_node.SetLock();
            r_node.FastGetSolutionStepValue(rDistanceVariable) = distance;
            r_node.Set(VISITED, true);
            r_node.UnSetLock();
        }
    });

    // no VISITED flag synchronization is required since, all ranks iterate
    // through same global indices thus setting VISITED flag in each rank

    KRATOS_CATCH("");
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

    return 0.0;

    KRATOS_CATCH("");
}

void RansWallDistanceCalculationProcess::ExecuteInitialize()
{
    KRATOS_TRY

    CalculateWallDistances();

    KRATOS_CATCH("");
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

    class NodalMinimumDistanceReduction
    {
    public:
        typedef std::tuple<std::vector<int>, std::vector<double>> value_type;
        std::vector<int> indices;
        std::vector<double> distances;

        /// access to reduced value
        value_type GetValue() const
        {
            return std::make_tuple(indices, distances);
        }

        /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
        void LocalReduce(const value_type& value)
        {
            const auto& value_indices = std::get<0>(value);
            const auto& value_distances = std::get<1>(value);

            for (std::size_t i = 0; i < value_indices.size(); ++i) {
                UpdateToMinimumDistances(indices, distances, value_indices[i],
                                         value_distances[i]);
            }
        }

        /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
        void ThreadSafeReduce(const NodalMinimumDistanceReduction& rOther)
        {
#pragma omp critical
            {
                const int number_of_indices = rOther.indices.size();
                for (int i = 0; i < number_of_indices; ++i) {
                    UpdateToMinimumDistances(indices, distances, rOther.indices[i],
                                             rOther.distances[i]);
                }
            }
        }
    };

    // first try setting distances for nodes based on condition normals
    std::vector<int> local_indices;
    std::vector<double> local_distances;

    const auto number_of_nodes = r_communicator.GlobalNumberOfNodes();
    local_indices.reserve(number_of_nodes);
    local_distances.reserve(number_of_nodes);

    std::tie(local_indices, local_distances) = block_for_each<NodalMinimumDistanceReduction>(
        r_model_part.Conditions(), [&](ModelPart::ConditionType& rCondition) {
            std::vector<int> indices;
            std::vector<double> distances;

            if (rCondition.Is(r_wall_flag) == mWallFlagVariableValue) {
                array_1d<double, 3> normal = rCondition.GetValue(NORMAL);
                normal /= norm_2(normal);

                auto& parent_geometry =
                    rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0].GetGeometry();
                for (auto& r_parent_node : parent_geometry) {
                    const int parent_node_id = r_parent_node.Id();
                    if (r_parent_node.Is(r_wall_flag) != mWallFlagVariableValue) {
                        const double wall_distance = inner_prod(
                            rCondition.GetGeometry().Center() - r_parent_node.Coordinates(),
                            normal);
                        UpdateToMinimumDistances(indices, distances,
                                                 parent_node_id, wall_distance);
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

            return std::make_tuple(indices, distances);
        });

    // communication for all ranks
    r_communicator.AssembleNonHistoricalData(NORMAL);

    // get condition based distances from all ranks
    UpdateModelPartDistancesToMinimum(r_model_part, local_indices, local_distances, r_distance_variable);

    VariableUtils().SetFlag(VISITED, false, r_model_part.Elements());

    // identify elements which has nodes needs distance calculated analytically.
    block_for_each(r_model_part.Elements(), [&](ModelPart::ElementType& rElement) {
        const auto& r_geometry = rElement.GetGeometry();
        bool found_non_wall_node_without_distance = false;

        array_1d<double, 3> normal = ZeroVector(3);
        array_1d<double, 3> normal_center = ZeroVector(3);
        double normals_count = 0.0;

        for (const auto& r_node : r_geometry) {
            if (r_node.Is(r_wall_flag) == mWallFlagVariableValue) {
                noalias(normal) += r_node.GetValue(NORMAL);
                noalias(normal_center) += r_node.Coordinates();
                normals_count += 1.0;
            } else {
                found_non_wall_node_without_distance = found_non_wall_node_without_distance
                                                           ? found_non_wall_node_without_distance
                                                           : !r_node.Is(VISITED);
            }
        }

        if (normals_count > 0.0 && found_non_wall_node_without_distance) {
            normal /= norm_2(normal);
            normal_center /= normals_count;
            rElement.SetValue(NORMAL, normal);
            rElement.SetValue(NODAL_VAUX, normal_center);
            rElement.Set(VISITED, true);
        }
    });

    // calculate distances for identified nodes in elements, which are not covered by conditions
    std::tie(local_indices, local_distances) = block_for_each<NodalMinimumDistanceReduction>(
        r_model_part.Elements(), [&](ModelPart::ElementType& rElement) {
            std::vector<int> indices;
            std::vector<double> distances;

            if (rElement.Is(VISITED)) {
                const auto& r_normal = rElement.GetValue(NORMAL);
                const auto& r_center = rElement.GetValue(NODAL_VAUX);
                const auto& r_geometry = rElement.GetGeometry();

                for (const auto& r_node : r_geometry) {
                    const int node_id = r_node.Id();
                    if (r_node.Is(r_wall_flag) != mWallFlagVariableValue &&
                        !r_node.Is(VISITED)) {
                        // assumes r_normal is always outward pointing, therefore wall_distance > 0.0
                        const double wall_distance =
                            inner_prod(r_center - r_node.Coordinates(), r_normal);
                        UpdateToMinimumDistances(indices, distances, node_id, wall_distance);
                    }
                }
            }

            return std::make_tuple(indices, distances);
        });

    // communication for all ranks
    UpdateModelPartDistancesToMinimum(r_model_part, local_indices, local_distances, r_distance_variable);

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
            "max_distance"                     : 1000.0,
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
