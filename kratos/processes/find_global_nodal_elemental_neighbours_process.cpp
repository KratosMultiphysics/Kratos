//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "processes/find_global_nodal_elemental_neighbours_process.h"

namespace Kratos
{

    void FindGlobalNodalElementalNeighboursProcess::Execute()
    {
        unsigned int current_rank = mrComm.Rank();
        auto& rNodes = mr_model_part.Nodes();

        //first of all the neighbour nodes and elements array are initialized to the guessed size
        //and empties the old entries
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rNodes.size()); ++i)
        {
            auto in = rNodes.begin() + i;
            in->SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector< Element >());
        }

        //compute the complete list of local neighbours
        for(auto& relem : mr_model_part.Elements())
        {
            GlobalPointer<Element> gpelem(&relem, current_rank);
            for(auto& node : relem.GetGeometry())
            {    
                node.GetValue(NEIGHBOUR_ELEMENTS).push_back( gpelem );
            }
        }

        if(mrComm.IsDistributed())
        {
            typedef std::unordered_map<int, GlobalPointersVector<Element>> neighbours_map_type; //contains id vs vector_of_neighbours 
            typedef std::unordered_map< int, neighbours_map_type > non_local_map_type;

            //construct the list of nodes that need to be sent
            non_local_map_type non_local_map;

            for(const auto& node : mr_model_part.GetCommunicator().InterfaceMesh().Nodes())
            {
                int owner_rank = node.FastGetSolutionStepValue(PARTITION_INDEX);
                non_local_map[owner_rank][node.Id()] = node.GetValue(NEIGHBOUR_ELEMENTS);
            }

            //here communicate non local data
            //compute communication plan
            std::vector<int> send_list;
            send_list.reserve( non_local_map.size() );
            for(auto& it : non_local_map)
                send_list.push_back( it.first );

            std::sort(send_list.begin(), send_list.end());
            auto colors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list, mrComm);

            //finalize computation of neighbour ids on owner nodes
            non_local_map_type recv_map;
            for(int color : colors)
            {
                if(color >= 0)
                {
                    //recev the global neighbours as computed on color
                    recv_map[color] = mrComm.SendRecv(non_local_map[color], color, color );

                    for(auto& r_item : recv_map[color])
                    {
                        auto& recv_node = mr_model_part.Nodes()[r_item.first]; 
                        auto& neighbours = recv_node.GetValue(NEIGHBOUR_ELEMENTS);
                        for(auto& gp : r_item.second.GetContainer())
                            neighbours.push_back(gp);
                    }
                }
            } //after this loop is finished neighbours are ok for the owner nodes

            //fill back the recv_map with the updated information
            for(int color : colors)
            {
                if(color >= 0)
                {
                    for(auto& r_item : recv_map[color])
                    {
                        //r_item.first contains the id of the node
                        //r_item.second contains the list of neighbours
                        auto& recv_node = mr_model_part.Nodes()[r_item.first]; 
                        r_item.second = recv_node.GetValue(NEIGHBOUR_ELEMENTS);
                    }

                    //obtain the final list of neighbours for nodes owned on color
                    auto final_gp_map = mrComm.SendRecv(recv_map[color], color, color );

                    //update the local database
                    for(auto& r_item : final_gp_map)
                    {
                        auto& recv_node = mr_model_part.Nodes()[r_item.first]; 
                        recv_node.GetValue(NEIGHBOUR_ELEMENTS) = r_item.second;
                    }
                }
            }
        }

        auto constructor_functor =  Kratos::ComputeNeighbourListFunctor<
                            ModelPart::NodesContainerType, 
                            Variable< GlobalPointersVector< Element > >
                            > (rNodes, NEIGHBOUR_ELEMENTS);

        GlobalPointerCommunicator<Element> pointer_comm(mrComm, constructor_functor );
        auto id_proxy = pointer_comm.Apply(
                [](GlobalPointer<Element> const& gp){return gp->Id();}
        );

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mr_model_part.Nodes().size()); ++i)
        {
            auto it = mr_model_part.NodesBegin() + i;
            auto& neighbours = it->GetValue(NEIGHBOUR_ELEMENTS);
            neighbours.shrink_to_fit();
            std::sort(neighbours.ptr_begin(), neighbours.ptr_end(),
                [&id_proxy](GlobalPointer<Element> const& gp1, GlobalPointer<Element> const& gp2)
                {
                    return id_proxy.Get(gp1) < id_proxy.Get(gp2);
                }
            );
        }
    
    }

    void FindGlobalNodalElementalNeighboursProcess::ClearNeighbours()
    {
        NodesContainerType& rNodes = mr_model_part.Nodes();

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rNodes.size()); ++i)
        {
            auto in = rNodes.begin() + i;
            auto& rN = in->GetValue(NEIGHBOUR_ELEMENTS);
            rN = GlobalPointersVector< Element >();
        }
    }

    std::unordered_map<int, std::vector<int> > FindGlobalNodalElementalNeighboursProcess::GetNeighbourIds(
                ModelPart::NodesContainerType& rNodes
                )
    {
        std::unordered_map<int, std::vector<int> > output;

        auto constructor_functor =  Kratos::ComputeNeighbourListFunctor<
                            ModelPart::NodesContainerType, 
                            Variable< GlobalPointersVector< Element > >
                            > (rNodes, NEIGHBOUR_ELEMENTS);

        GlobalPointerCommunicator<Element> pointer_comm(mrComm, constructor_functor );

        auto result_proxy = pointer_comm.Apply(
                [](GlobalPointer<Element>& gp){return gp->Id();}
        );

        #pragma omp parallel for
        for(int k=0; k<static_cast<int>(rNodes.size()); ++k)
        {
            auto in = rNodes.begin() + k;
            auto& neighbours = in->GetValue(NEIGHBOUR_ELEMENTS);
            std::vector<int> tmp(neighbours.size());
            for(unsigned int i=0; i<neighbours.size(); ++i)
            {
                tmp[i] = result_proxy.Get(neighbours(i));
            }
            output[in->Id()] = tmp;  
        }


        return output;
    }



}  // namespace Kratos.


