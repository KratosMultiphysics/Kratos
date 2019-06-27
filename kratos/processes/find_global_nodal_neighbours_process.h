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


#if !defined(KRATOS_FIND_GLOBAL_NODAL_NEIGHBOURS_PROCESS_H_INCLUDED )
#define  KRATOS_FIND_GLOBAL_NODAL_NEIGHBOURS_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/global_pointer_variables.h"
#include "utilities/communication_coloring_utilities.h"
#include "utilities/pointer_communicator.h"
#include "utilities/global_pointer_utilities.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
typedef  ModelPart::NodesContainerType NodesContainerType;
typedef  ModelPart::ElementsContainerType ElementsContainerType;


///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class FindGlobalNodalNeighboursProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindGlobalNodalNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindGlobalNodalNeighboursProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /// avg_elems ------ expected number of neighbour elements per node.,
    /// avg_nodes ------ expected number of neighbour Nodes
    /// the better the guess for the quantities above the less memory occupied and the fastest the algorithm
    FindGlobalNodalNeighboursProcess(const DataCommunicator& rComm, ModelPart& model_part, unsigned int avg_nodes = 10)
        : mrComm(rComm),mr_model_part(model_part)
    {
        mavg_nodes = avg_nodes;
    }

    /// Destructor.
    ~FindGlobalNodalNeighboursProcess() override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        NodesContainerType& rNodes = mr_model_part.Nodes();

        //first of all the neighbour nodes and elements array are initialized to the guessed size
        //and empties the old entries
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            auto& rN = in->GetValue(NEIGHBOUR_NODES);
            rN = GlobalPointersVector< Node<3> >();
        }

        //adding the neighbouring nodes
        if(mrComm.IsDistributed() == false)
        {
            for(auto& relem : mr_model_part.Elements())
            {
                auto& rgeom = relem.GetGeometry();
                for(unsigned int i = 0; i < rgeom.size(); i++)
                {
                    for(unsigned int j = 0; j < rgeom.size(); j++)
                    {
                        if( j != i )
                        {
                            auto gp = GlobalPointer<Node<3>>(rgeom(j),0);
                            AddUniqueGlobalPointer< Node<3> >(rgeom[i].GetValue(NEIGHBOUR_NODES), gp);
                        }
                    }
                }
            }

            //TODO: sort according to nodal Id
            #pragma omp parallel for
            for(int i=0; i<static_cast<int>(mr_model_part.Nodes().size()); ++i)
            {
                auto it = mr_model_part.NodesBegin() + i;
                auto& neighbours = it->GetValue(NEIGHBOUR_NODES);
                neighbours.shrink_to_fit();
                std::sort(neighbours.ptr_begin(), neighbours.ptr_end(),
                          [](GlobalPointer<Node<3>> const& gp1, GlobalPointer<Node<3>> const& gp2)
                            {
                                return gp1->Id() < gp2->Id();
                            }
                         );

            }
        }
        else //mpi case!
        {
            int current_rank = mrComm.Rank();

            typedef std::unordered_map< int, std::vector<int> > map_of_sets;
            std::unordered_map<int, map_of_sets> neighbours_ids;

            for(auto& relem : mr_model_part.Elements())
            {
                const auto& rgeom = relem.GetGeometry();
                for(unsigned int i = 0; i < rgeom.size(); i++)
                {
                    const int i_owner_rank = rgeom[i].FastGetSolutionStepValue(PARTITION_INDEX);
                    auto& container = neighbours_ids[i_owner_rank][rgeom[i].Id()];
                    for(unsigned int j = 0; j < rgeom.size(); j++)
                        if( j != i  )
                            AddUnique(container,rgeom[j].Id());
                }
            }

            //here communicate non local data
            //compute communication plan
            std::vector<int> send_list;
            send_list.reserve( neighbours_ids.size() );
            for(auto& it : neighbours_ids)
                if(it.first != current_rank)
                    send_list.push_back( it.first );

            std::sort(send_list.begin(), send_list.end());
            auto colors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list, mrComm);

            //finalize computation of neighbour ids on owner nodes
            std::unordered_map<int, std::vector<int> > non_local_node_ids; //this will contain the id of the nodes that will need communicaiton
            for(int color : colors)
            {
                if(color >= 0)
                {
                    auto tmp = mrComm.SendRecv(neighbours_ids[color], color, color );
                    for(auto& item : tmp)
                    {
                        auto& ids = neighbours_ids[current_rank][item.first];
                        for(int neighbour_id : item.second)
                            AddUnique(ids,neighbour_id);
                        non_local_node_ids[color].push_back(item.first); //this are the nodes (ids) for which neihbours are needed
                    }
                }
            }

            for(auto& owner : neighbours_ids)
            {
                for(auto& item : owner.second)
                {
                    std::sort(item.second.begin(), item.second.end());
                    auto last = std::unique(item.second.begin(), item.second.end());
                    item.second.erase(last, item.second.end());
                }
            }

            //obtain all global pointers needed
            std::vector<int> all_ids;
            for(const auto& item : neighbours_ids[current_rank])
            {
                all_ids.push_back(item.first);
                for(int id : item.second)
                    all_ids.push_back(id);
            }
            std::sort(all_ids.begin(), all_ids.end());
            auto last = std::unique(all_ids.begin(), all_ids.end());
            all_ids.erase(last, all_ids.end());

            auto all_gps_map = GlobalPointerUtilities::RetrieveGlobalIndexedPointersMap(mr_model_part.Nodes(), all_ids, mrComm);

            //now construct the list of GlobalPointers - here neighbours are ok for locally owned nodes
            for(const auto& item : neighbours_ids[current_rank])
            {
                auto node_id = item.first;
                auto& r_node = mr_model_part.Nodes()[node_id];
                auto& neighbours = r_node.GetValue(NEIGHBOUR_NODES);
                neighbours.reserve(item.second.size());

                for(int id : item.second)
                {
                    auto found = all_gps_map.find(id);
                    KRATOS_DEBUG_ERROR_IF(found == all_gps_map.end()) << "id " << id << " not found in all_gps_map" << std::endl;

                    neighbours.push_back(found->second);
                }
                neighbours.shrink_to_fit();

            }

            //finalize computation by obtaining the neighbours for non-local nodes
            for(int color : colors)
            {
                if(color >= 0)
                {
                    std::unordered_map<int, GlobalPointersVector<Node<3>> > neighbours_to_send;

                    for(auto id : non_local_node_ids[color])
                    {
                        neighbours_to_send[id] = mr_model_part.Nodes()[id].GetValue(NEIGHBOUR_NODES);
                    }
                    auto received_neighbours = mrComm.SendRecv(neighbours_to_send, color, color );
                    for(auto& item : received_neighbours)
                    {
                        auto& r_node = mr_model_part.Nodes()[item.first];
                        r_node.SetValue(NEIGHBOUR_NODES, item.second);
                    }
                }
            }
        }



    }

    void ClearNeighbours()
    {
        NodesContainerType& rNodes = mr_model_part.Nodes();
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            auto& rN = in->GetValue(NEIGHBOUR_NODES);
            rN.erase(rN.begin(),rN.end() );
            rN.shrink_to_fit();
        }
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "FindGlobalNodalNeighboursProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindGlobalNodalNeighboursProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    std::unordered_map<int, std::vector<int> > GetNeighbourIds(
                ModelPart::NodesContainerType& rNodes
                )
    {
        std::unordered_map<int, std::vector<int> > output;

        GlobalPointersVector< Node<3> > gp_list;
        for(auto& node : rNodes)
            for(auto& gp : node.GetValue(NEIGHBOUR_NODES).GetContainer())
                gp_list.push_back(gp);
        gp_list.Unique();

        GlobalPointerCommunicator<Node<3>> pointer_comm(mrComm, gp_list);
        auto result_proxy = pointer_comm.Apply(
                [](GlobalPointer<Node<3>>& gp){return gp->Id();}
        );

        for(auto& node : rNodes)
        {
            auto& neighbours = node.GetValue(NEIGHBOUR_NODES);
            std::vector<int> tmp(neighbours.size());
            for(unsigned int i=0; i<neighbours.size(); ++i)
            {
                tmp[i] = result_proxy.Get(neighbours(i));
            }
            output[node.Id()] = tmp;
        }

        return output;
    }



    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    const DataCommunicator& mrComm;
    ModelPart& mr_model_part;
    unsigned int mavg_nodes;


    ///@}
    ///@name Private Operators
    ///@{

    //******************************************************************************************
    //******************************************************************************************
    void AddUnique(std::vector<int>& container, const int item)
    {
        bool found = false;
        for(const int& i : container)
        {
            if(i == item)
            {
                found = true;
                break;
            }
        }
        if(!found)
            container.push_back(item);
    }

    template< class TDataType > void  AddUniqueGlobalPointer
    (GlobalPointersVector< TDataType >& v, const GlobalPointer< TDataType >& candidate)
    {
        bool found = false;
        for(auto& gp : v.GetContainer())
        {
            if(&(*gp) == &(*candidate) && gp.GetRank() == candidate.GetRank())
            {
                found = true;
                break;
            }
        }
        if(!found)
            v.push_back(candidate);
    }

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    FindGlobalNodalNeighboursProcess& operator=(FindGlobalNodalNeighboursProcess const& rOther);

    /// Copy constructor.
    //FindGlobalNodalNeighboursProcess(FindGlobalNodalNeighboursProcess const& rOther);


    ///@}

}; // Class FindGlobalNodalNeighboursProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FindGlobalNodalNeighboursProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FindGlobalNodalNeighboursProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FIND_GLOBAL_NODAL_NEIGHBOURS_PROCESS_H_INCLUDED  defined


