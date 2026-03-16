//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Suneth Warnakulasuriya
//

// System includes
#include <iostream>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/global_pointer_variables.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "processes/process.h"
#include "utilities/communication_coloring_utilities.h"
#include "utilities/global_pointer_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/pointer_communicator.h"
#include "utilities/variable_utils.h"

// Include base h
#include "find_global_nodal_neighbours_for_entities_process.h"

namespace Kratos
{

template <class TContainerType>
TContainerType& FindNodalNeighboursForEntitiesProcess<TContainerType>::GetContainer()
{
    if constexpr (std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        return this->mrModelPart.Elements();
    } else if constexpr (std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        return this->mrModelPart.Conditions();
    } else {
        KRATOS_ERROR << "Unsupported container type" << std::endl;
    }
}

template <class TContainerType>
void FindNodalNeighboursForEntitiesProcess<TContainerType>::AddHangingNodeIds(
    std::unordered_map<int, std::unordered_map<int, std::vector<int>>>& rNeighbourIds) const
{
    // do nothing for elements here since mettis partitioner is based on elements, there
    // cannot be any hanging nodes. If metis partitioner is based on conditions, then
    // this process need not to be used with elements, so this method won't be called.

    if constexpr (std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        // if the metis partitioner is based on conditions, this will still work,
        // but with this additional cost of checking.

        // this loop cannot run in parallel, since std::unordered_maps adds the key if it
        // is not found, otherwise do nothing.
        for (const NodeType& rNode : mrModelPart.Nodes()) {
            const int i_owner_rank = rNode.FastGetSolutionStepValue(PARTITION_INDEX);
            const auto node_id = rNode.Id();
            rNeighbourIds[i_owner_rank][node_id];
        }
    } else {
        KRATOS_ERROR << "Unsupported container type" << std::endl;
    }
}

template <class TContainerType>
void FindNodalNeighboursForEntitiesProcess<TContainerType>::Execute()
{
    KRATOS_TRY

    auto& r_nodes = this->mrModelPart.Nodes();

    // first of all the neighbour nodes and elements array are initialized
    // to the guessed size and empties the old entries
    VariableUtils().SetNonHistoricalVariable(
        mrOutputVariable, GlobalPointersVector<NodeType>(), r_nodes);

    // adding the neighbouring nodes
    if (!this->mrDataCommunicator.IsDistributed()) {
        block_for_each(this->GetContainer(), [&](typename TContainerType::value_type& rEntity) {
            auto& r_geometry = rEntity.GetGeometry();
            for (unsigned int i = 0; i < r_geometry.size(); ++i) {
                for (unsigned int j = 0; j < r_geometry.size(); ++j) {
                    if (j != i) {
                        auto gp = GlobalPointer<NodeType>(r_geometry(j), 0);
                        auto& r_node = r_geometry[i];
                        r_node.SetLock();
                        AddUniqueGlobalPointer<NodeType>(
                            r_node.GetValue(mrOutputVariable), gp);
                        r_node.UnSetLock();
                    }
                }
            }
        });

        block_for_each(r_nodes, [&](NodeType& rNode) {
            auto& r_neighbours = rNode.GetValue(mrOutputVariable);
            r_neighbours.shrink_to_fit();
            std::sort(r_neighbours.ptr_begin(), r_neighbours.ptr_end(),
                      [](GlobalPointer<NodeType> const& gp1,
                         GlobalPointer<NodeType> const& gp2) {
                          return gp1->Id() < gp2->Id();
                      });
        });

    } else { // mpi case!
        const int current_rank = this->mrDataCommunicator.Rank();

        using map_of_sets = std::unordered_map<int, std::vector<int>>;
        std::unordered_map<int, map_of_sets> neighbours_ids;

        for (auto& r_entity : this->GetContainer()) {
            const auto& r_geometry = r_entity.GetGeometry();
            for (unsigned int i = 0; i < r_geometry.size(); ++i) {
                const int i_owner_rank =
                    r_geometry[i].FastGetSolutionStepValue(PARTITION_INDEX);
                auto& container = neighbours_ids[i_owner_rank][r_geometry[i].Id()];
                for (unsigned int j = 0; j < r_geometry.size(); ++j) {
                    if (j != i) {
                        AddUnique(container, r_geometry[j].Id());
                    }
                }
            }
        }

        // there are some isolated nodes when metis partitioner performs
        // partitioning specially in the case where TContainerType =
        // ModelPart::ConditionsContainerType therefore we add those ids to
        // neighbours_ids so that proper communication scheduling can be computed.
        AddHangingNodeIds(neighbours_ids);

        // here communicate non local data
        // compute communication plan
        std::vector<int> send_list;
        send_list.reserve(neighbours_ids.size());
        for (const auto& it : neighbours_ids) {
            if (it.first != current_rank) {
                send_list.push_back(it.first);
            }
        }

        std::sort(send_list.begin(), send_list.end());
        const auto colors = MPIColoringUtilities::ComputeCommunicationScheduling(
            send_list, this->mrDataCommunicator);

        // finalize computation of neighbour ids on owner nodes
        std::unordered_map<int, std::vector<int>> non_local_node_ids; // this will contain the id of the nodes that will
                                                                      // need communicaiton
        for (const int color : colors) {
            if (color >= 0) {
                auto tmp = this->mrDataCommunicator.SendRecv(
                    neighbours_ids[color], color, color);
                for (const auto& item : tmp) {
                    auto& ids = neighbours_ids[current_rank][item.first];
                    for (int neighbour_id : item.second)
                        AddUnique(ids, neighbour_id);
                    non_local_node_ids[color].push_back(
                        item.first); // this are the nodes (ids) for which
                                     // neihbours are needed
                }
            }
        }

        for (auto& owner : neighbours_ids) {
            for (auto& item : owner.second) {
                std::sort(item.second.begin(), item.second.end());
                auto last = std::unique(item.second.begin(), item.second.end());
                item.second.erase(last, item.second.end());
            }
        }

        // obtain all global pointers needed
        std::vector<int> all_ids;
        for (const auto& item : neighbours_ids[current_rank]) {
            all_ids.push_back(item.first);
            for (int id : item.second) {
                all_ids.push_back(id);
            }
        }
        std::sort(all_ids.begin(), all_ids.end());
        auto last = std::unique(all_ids.begin(), all_ids.end());
        all_ids.erase(last, all_ids.end());

        auto all_gps_map = GlobalPointerUtilities::RetrieveGlobalIndexedPointersMap(
            this->mrModelPart.Nodes(), all_ids, this->mrDataCommunicator);

        // now construct the list of GlobalPointers - here neighbours are ok for
        // locally owned nodes
        for (const auto& item : neighbours_ids[current_rank]) {
            auto node_id = item.first;
            auto& r_node = this->mrModelPart.Nodes()[node_id];
            auto& neighbours = r_node.GetValue(mrOutputVariable);
            neighbours.reserve(item.second.size());

            for (int id : item.second) {
                auto found = all_gps_map.find(id);
                KRATOS_DEBUG_ERROR_IF(found == all_gps_map.end())
                    << "id " << id << " not found in all_gps_map" << std::endl;

                neighbours.push_back(found->second);
            }
            neighbours.shrink_to_fit();
        }

        // finalize computation by obtaining the neighbours for non-local nodes
        for (const int color : colors) {
            if (color >= 0) {
                std::unordered_map<int, GlobalPointersVector<NodeType>> neighbours_to_send;

                for (auto id : non_local_node_ids[color]) {
                    neighbours_to_send[id] =
                        this->mrModelPart.Nodes()[id].GetValue(mrOutputVariable);
                }
                auto received_neighbours = this->mrDataCommunicator.SendRecv(
                    neighbours_to_send, color, color);
                for (auto& item : received_neighbours) {
                    auto& r_node = this->mrModelPart.Nodes()[item.first];
                    r_node.SetValue(mrOutputVariable, item.second);
                }
            }
        }
    }

    KRATOS_CATCH("");
}

template <class TContainerType>
void FindNodalNeighboursForEntitiesProcess<TContainerType>::ClearNeighbours()
{
    block_for_each(this->mrModelPart.Nodes(), [&](NodeType& rNode) {
        if(rNode.Has(mrOutputVariable)){
            auto& r_gp_neighbour_nodes_vector = rNode.GetValue(mrOutputVariable);
            r_gp_neighbour_nodes_vector.erase(r_gp_neighbour_nodes_vector.begin(),
                                            r_gp_neighbour_nodes_vector.end());
            r_gp_neighbour_nodes_vector.shrink_to_fit();
        } else {
            rNode.SetValue(mrOutputVariable, GlobalPointersVector<NodeType>());
        }
    });
}

template <class TContainerType>
std::unordered_map<int, std::vector<int>> FindNodalNeighboursForEntitiesProcess<TContainerType>::GetNodalNeighbourIdsMap(
    ModelPart::NodesContainerType& rNodes,
    const DataCommunicator& rDataCommunicator,
    const Variable<GlobalPointersVector<NodeType>>& rNodalGlobalPointerVariable)
{
    class GlobalPointerAdder
    {
    public:
        typedef GlobalPointersVector<NodeType> value_type;
        typedef GlobalPointersVector<NodeType> return_type;

        return_type gp_vector;
        return_type GetValue()
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

    GlobalPointersVector<NodeType> all_global_pointers =
        block_for_each<GlobalPointerAdder>(rNodes, [&](NodeType& rNode) {
            return rNode.GetValue(rNodalGlobalPointerVariable);
        });

    GlobalPointerCommunicator<NodeType> pointer_comm(rDataCommunicator, all_global_pointers);

    auto node_id_proxy = pointer_comm.Apply(
        [](const GlobalPointer<NodeType>& rGP) { return rGP->Id(); });

    std::unordered_map<int, std::vector<int>> output;

    for (const auto& r_node : rNodes) {
        const auto& r_neighbours = r_node.GetValue(rNodalGlobalPointerVariable);
        std::vector<int> neighbour_id(r_neighbours.size());
        for (unsigned int i = 0; i < r_neighbours.size(); ++i) {
            neighbour_id[i] = node_id_proxy.Get(r_neighbours(i));
        }
        output[r_node.Id()] = neighbour_id;
    }

    return output;
}

template <class TContainerType>
void FindNodalNeighboursForEntitiesProcess<TContainerType>::AddUnique(
    std::vector<int>& rContainer,
    const int Item) const
{
    if (std::find(rContainer.begin(), rContainer.end(), Item) == rContainer.end()) {
        rContainer.push_back(Item);
    }
}

// template instantiations

template class KRATOS_API(KRATOS_CORE) FindNodalNeighboursForEntitiesProcess<ModelPart::ElementsContainerType>;
template class KRATOS_API(KRATOS_CORE) FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType>;

} // namespace Kratos.
