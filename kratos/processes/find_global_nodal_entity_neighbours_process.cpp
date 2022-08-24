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

// External includes

// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/element.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"
#include "utilities/compute_neighbour_list_functor.h"
#include "utilities/communication_coloring_utilities.h"
#include "utilities/pointer_communicator.h"
#include "utilities/global_pointer_utilities.h"

// Include base h
#include "processes/find_global_nodal_entity_neighbours_process.h"

namespace Kratos
{

template <>
ModelPart::ConditionsContainerType& FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType>::GetContainer()
{
    return this->mrModelPart.Conditions();
}

template <>
ModelPart::ElementsContainerType& FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType>::GetContainer()
{
    return this->mrModelPart.Elements();
}

template<>
FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType>::FindGlobalNodalEntityNeighboursProcess(
    ModelPart& rModelPart)
    : mrModelPart(rModelPart),
      mrOutputVariable(NEIGHBOUR_CONDITIONS)
{
}

template<>
FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType>::FindGlobalNodalEntityNeighboursProcess(
    ModelPart& rModelPart)
    : mrModelPart(rModelPart),
      mrOutputVariable(NEIGHBOUR_ELEMENTS)
{
}

template<class TContainerType>
FindGlobalNodalEntityNeighboursProcess<TContainerType>::FindGlobalNodalEntityNeighboursProcess(
    ModelPart& rModelPart,
    const Variable<GlobalEntityPointersVectorType>& rOutputVariable)
    : mrModelPart(rModelPart),
      mrOutputVariable(rOutputVariable)
{
}

template<class TContainerType>
void FindGlobalNodalEntityNeighboursProcess<TContainerType>::Execute()
{
    KRATOS_TRY

    // first of all the neighbour nodes and elements array are initialized to the guessed size
    // and empties the old entries
    auto& r_nodes = mrModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(mrOutputVariable, GlobalEntityPointersVectorType(), r_nodes);

    // compute the complete list of local neighbours
    const auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();
    const IndexType current_rank = r_data_communicator.Rank();

    block_for_each(GetContainer(), [&](EntityType& rEntity) {
        GlobalPointer<EntityType> gp_entity(&rEntity, current_rank);
        for(auto& r_node : rEntity.GetGeometry()) {
            r_node.SetLock();
            r_node.GetValue(mrOutputVariable).push_back(gp_entity);
            r_node.UnSetLock();
        }
    });

    if(r_data_communicator.IsDistributed()) {
        // construct the list of nodes that need to be sent
        NonLocalMapType non_local_map;

        for(const auto& r_node : mrModelPart.GetCommunicator().InterfaceMesh().Nodes()) {
            const int owner_rank = r_node.FastGetSolutionStepValue(PARTITION_INDEX);
            non_local_map[owner_rank][r_node.Id()] = r_node.GetValue(mrOutputVariable);
        }

        // here communicate non local data
        // compute communication plan
        std::vector<int> send_list;
        send_list.reserve(non_local_map.size());
        for(auto& it : non_local_map)
            send_list.push_back(it.first);

        std::sort(send_list.begin(), send_list.end());
        auto colors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list, r_data_communicator);

        // finalize computation of neighbour ids on owner nodes
        NonLocalMapType recv_map;
        for(int color : colors) {
            if (color >= 0) {
                // recev the global neighbours as computed on color
                recv_map[color] = r_data_communicator.SendRecv(non_local_map[color], color, color);

                for(auto& r_item : recv_map[color]) {
                    auto& recv_node = mrModelPart.GetNode(r_item.first);
                    auto& neighbours = recv_node.GetValue(mrOutputVariable);
                    for(auto& r_gp : r_item.second.GetContainer()) {
                        neighbours.push_back(r_gp);
                    }
                }
            }
        } // after this loop is finished neighbours are ok for the owner nodes

        // fill back the recv_map with the updated information
        for(int color : colors) {
            if (color >= 0) {
                for(auto& r_item : recv_map[color]) {
                    //r_item.first contains the id of the node
                    //r_item.second contains the list of neighbours
                    auto& r_recv_node = mrModelPart.GetNode(r_item.first);
                    r_item.second = r_recv_node.GetValue(mrOutputVariable);
                }

                //obtain the final list of neighbours for nodes owned on color
                auto final_gp_map = r_data_communicator.SendRecv(recv_map[color], color, color );

                //update the local database
                for(auto& r_item : final_gp_map) {
                    auto& r_recv_node = mrModelPart.GetNode(r_item.first);
                    r_recv_node.GetValue(mrOutputVariable) = r_item.second;
                }
            }
        }
    }

    auto constructor_functor =  ComputeNeighbourListFunctor<ModelPart::NodesContainerType, Variable<GlobalEntityPointersVectorType>>(r_nodes, mrOutputVariable);

    GlobalPointerCommunicator<EntityType> pointer_comm(r_data_communicator, constructor_functor);
    auto id_proxy = pointer_comm.Apply(
            [](GlobalPointer<EntityType> const& gp){return gp->Id();}
    );

    block_for_each(r_nodes, [&](NodeType& rNode){
        auto& neighbours = rNode.GetValue(mrOutputVariable);
        neighbours.shrink_to_fit();
        std::sort(neighbours.ptr_begin(), neighbours.ptr_end(),
            [&id_proxy](GlobalPointer<EntityType> const& gp1, GlobalPointer<EntityType> const& gp2)
            {
                return id_proxy.Get(gp1) < id_proxy.Get(gp2);
            }
        );
    });

    KRATOS_CATCH("");
}

template<class TContainerType>
void FindGlobalNodalEntityNeighboursProcess<TContainerType>::Clear()
{
    auto& rNodes = mrModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(mrOutputVariable, GlobalEntityPointersVectorType(), rNodes);
}

template<class TContainerType>
std::unordered_map<int, std::vector<int>> FindGlobalNodalEntityNeighboursProcess<TContainerType>::GetNeighbourIds(
    ModelPart::NodesContainerType& rNodes)
{
    KRATOS_TRY

    std::unordered_map<int, std::vector<int>> output;

    auto constructor_functor =  Kratos::ComputeNeighbourListFunctor<ModelPart::NodesContainerType, Variable<GlobalEntityPointersVectorType>>(rNodes, mrOutputVariable);

    GlobalPointerCommunicator<EntityType> pointer_comm(mrModelPart.GetCommunicator().GetDataCommunicator(), constructor_functor);

    auto result_proxy = pointer_comm.Apply(
        [](GlobalPointer<EntityType>& gp) {
            return gp->Id();
        }
    );

    block_for_each(rNodes, [&](NodeType& rNode){
        auto& r_neighbours = rNode.GetValue(mrOutputVariable);
        std::vector<int> tmp(r_neighbours.size());
        for(unsigned int i=0; i<r_neighbours.size(); ++i) {
            tmp[i] = result_proxy.Get(r_neighbours(i));
        }
        output[rNode.Id()] = tmp;
    });

    return output;

    KRATOS_CATCH("");
}

// template instantiations
template class FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType>;
template class FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType>;

}  // namespace Kratos.


