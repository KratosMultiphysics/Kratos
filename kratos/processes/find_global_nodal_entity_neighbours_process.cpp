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
#include "utilities/reduction_utilities.h"
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
ModelPart::ConditionsContainerType& FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType>::GetContainer(ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

template <>
ModelPart::ElementsContainerType& FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType>::GetContainer(ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template <>
const Variable<FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType>::GlobalEntityPointersVectorType>& FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType>::GetDefaultOutputVariable()
{
    return NEIGHBOUR_CONDITIONS;
}

template <>
const Variable<FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType>::GlobalEntityPointersVectorType>& FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType>::GetDefaultOutputVariable()
{
    return NEIGHBOUR_ELEMENTS;
}

template<class TContainerType>
const Parameters FindGlobalNodalEntityNeighboursProcess<TContainerType>::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name"     : "PLEASE_SPECIFY_MODEL_PART_NAME"
        })");

    return default_parameters;
}

template<class TContainerType>
FindGlobalNodalEntityNeighboursProcess<TContainerType>::FindGlobalNodalEntityNeighboursProcess(
    Model& rModel,
    Parameters Params)
    : mrModel(rModel),
      mrOutputVariable(GetDefaultOutputVariable())
{
    Params.ValidateAndAssignDefaults(GetDefaultParameters());
    mModelPartName = Params["model_part_name"].GetString();
}

template<>
FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType>::FindGlobalNodalEntityNeighboursProcess(
    ModelPart& rModelPart)
    : mrModel(rModelPart.GetModel()),
      mrOutputVariable(NEIGHBOUR_CONDITIONS)
{
    mModelPartName = rModelPart.FullName();
}

template<>
FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType>::FindGlobalNodalEntityNeighboursProcess(
    ModelPart& rModelPart)
    : mrModel(rModelPart.GetModel()),
      mrOutputVariable(NEIGHBOUR_ELEMENTS)
{
    mModelPartName = rModelPart.FullName();
}

template<class TContainerType>
FindGlobalNodalEntityNeighboursProcess<TContainerType>::FindGlobalNodalEntityNeighboursProcess(
    ModelPart& rModelPart,
    const Variable<GlobalEntityPointersVectorType>& rOutputVariable)
    : mrModel(rModelPart.GetModel()),
      mrOutputVariable(rOutputVariable)
{
    mModelPartName = rModelPart.FullName();
}

template<class TContainerType>
void FindGlobalNodalEntityNeighboursProcess<TContainerType>::Execute()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    // first of all the neighbour nodes and elements array are initialized to the guessed size
    // and empties the old entries
    auto& r_nodes = r_model_part.Nodes();
    VariableUtils().SetNonHistoricalVariable(mrOutputVariable, GlobalEntityPointersVectorType(), r_nodes);

    // compute the complete list of local neighbours
    const auto& r_data_communicator = r_model_part.GetCommunicator().GetDataCommunicator();
    const IndexType current_rank = r_data_communicator.Rank();

    block_for_each(GetContainer(r_model_part), [&](EntityType& rEntity) {
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

        for(const auto& r_node : r_model_part.GetCommunicator().InterfaceMesh().Nodes()) {
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
                    auto& recv_node = r_model_part.GetNode(r_item.first);
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
                    auto& r_recv_node = r_model_part.GetNode(r_item.first);
                    r_item.second = r_recv_node.GetValue(mrOutputVariable);
                }

                //obtain the final list of neighbours for nodes owned on color
                auto final_gp_map = r_data_communicator.SendRecv(recv_map[color], color, color );

                //update the local database
                for(auto& r_item : final_gp_map) {
                    auto& r_recv_node = r_model_part.GetNode(r_item.first);
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
    auto& rNodes = mrModel.GetModelPart(mModelPartName).Nodes();
    VariableUtils().SetNonHistoricalVariable(mrOutputVariable, GlobalEntityPointersVectorType(), rNodes);
}

template<class TContainerType>
typename FindGlobalNodalEntityNeighboursProcess<TContainerType>::IdMapType FindGlobalNodalEntityNeighboursProcess<TContainerType>::GetNeighbourIds(const ModelPart::NodesContainerType& rNodes)
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    auto constructor_functor =  Kratos::ComputeNeighbourListFunctor<ModelPart::NodesContainerType, Variable<GlobalEntityPointersVectorType>>(rNodes, mrOutputVariable);

    GlobalPointerCommunicator<EntityType> pointer_comm(r_model_part.GetCommunicator().GetDataCommunicator(), constructor_functor);

    auto result_proxy = pointer_comm.Apply(
        [](GlobalPointer<EntityType>& gp) {
            return gp->Id();
        }
    );

    return block_for_each<MapReduction<IdMapType>>(rNodes, [&](auto& rNode){
        auto& r_neighbours = rNode.GetValue(mrOutputVariable);
        std::vector<int> tmp(r_neighbours.size());
        for(unsigned int i=0; i<r_neighbours.size(); ++i) {
            tmp[i] = result_proxy.Get(r_neighbours(i));
        }
        return std::make_pair(rNode.Id(), tmp);
    });

    KRATOS_CATCH("");
}

// template instantiations
template class FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType>;
template class FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType>;

}  // namespace Kratos.
