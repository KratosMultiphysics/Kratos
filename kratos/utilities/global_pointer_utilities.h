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


#if !defined(KRATOS_GLOBAL_POINTER_UTILITIES_H_INCLUDED )
#define  KRATOS_GLOBAL_POINTER_UTILITIES_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/data_communicator.h"
#include "includes/global_pointer.h"
#include "containers/global_pointers_vector.h"
#include "containers/global_pointers_unordered_map.h"

#include "utilities/communication_coloring_utilities.h"

namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

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
class GlobalPointerUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GlobalPointerUtilities
    KRATOS_CLASS_POINTER_DEFINITION(GlobalPointerUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GlobalPointerUtilities()
    {}

    /// Destructor.
    virtual ~GlobalPointerUtilities() {}

    template< class TContainerType >
    static std::unordered_map< int, GlobalPointer<typename TContainerType::value_type> > RetrieveGlobalIndexedPointersMap(
        const TContainerType& container,
        const std::vector<int>& id_list,
        const DataCommunicator& rDataCommunicator)
    {
        typedef GlobalPointer<typename TContainerType::value_type> GPType;

        std::unordered_map< int, GPType > global_pointers_list;
        const int current_rank = rDataCommunicator.Rank();
        const int world_size = rDataCommunicator.Size();

        std::vector<int> remote_ids;

        if(rDataCommunicator.IsDistributed()) {
            // If the execution is distributed, for every entity id find if its in the container, and if it is, if it is local
            // to our partition.
            for(const int id : id_list ) {
                const auto it = container.find(id);

                if( it != container.end()) { 
                    if(ObjectIsLocal(*it, current_rank)) {
                        // Found locally
                        global_pointers_list.emplace(id,GPType(&*it, current_rank));
                    } else {
                        // Remote, but this is a lucky case since for those we know to which rank they  they belong
                        // TODO: optimize according to the comment just above
                        remote_ids.push_back(id);
                    }
                } else {
                    // Id not found and we have no clue of what node owns it
                    remote_ids.push_back(id);
                }
            }
        } else {
            // If the execution is not distributed, only check if the id is in the container.
            for(const int id : id_list ) {
                const auto it = container.find(id);
                if( it != container.end()) { 
                    // Found locally
                    global_pointers_list.emplace(id,GPType(&*it, current_rank));
                }
            }
        }

        // Gather everything onto master_rank processor
        int master_rank = 0;

        std::vector<int> all_remote_ids;
        std::vector< std::vector<int> > collected_remote_ids(world_size);
        std::unordered_map< int, GPType > all_non_local_gp_map;

        //STEP1 - here we send the id we need to the master_rank
        //NOTE: here we DO NOT use a collective since we need to keep distinguished the ids per rank
        for(int i=0; i<world_size; ++i)
        {
            if(i != master_rank)
            {
                if(current_rank == master_rank) //only master executes
                {
                    rDataCommunicator.Recv(collected_remote_ids[i],i);
                }
                else if(current_rank == i) //only processor i executes
                {
                    rDataCommunicator.Send(remote_ids,master_rank);
                }
            }
            else //no communication needed
            {
                if(current_rank == master_rank) //only master executes
                    collected_remote_ids[i] = remote_ids;
            }

            if(current_rank == master_rank)
            {
                for(const int id : collected_remote_ids[i])
                    all_remote_ids.push_back( id );
            }
        }

        // very useful for debugging. do not remove for now
        // if(current_rank == master_rank)
        // {
        //     std::cout << "collected ids " << std::endl;
        //     for(unsigned int rank=0; rank<collected_remote_ids.size(); ++rank)
        //     {
        //         std::cout << " r = " << rank << " - ";
        //         for(int id : collected_remote_ids[rank])
        //             std::cout << id << " " ;
        //         std::cout << std::endl;
        //     }

        //     std::cout << "all remote ids " << std::endl;
        //     for(int id : all_remote_ids)
        //         std::cout << id << " ";
        //     std::cout << std::endl;
        // }

        if(current_rank == master_rank)
        {
            std::sort(all_remote_ids.begin(), all_remote_ids.end());
            std::unique(all_remote_ids.begin(), all_remote_ids.end());
        }

        //communicate the size of all remote_ids and resize the vector accordingly
        int number_of_all_remote_ids = all_remote_ids.size();
        rDataCommunicator.Broadcast(number_of_all_remote_ids,master_rank);

        if(current_rank != master_rank)
            all_remote_ids.resize(number_of_all_remote_ids);

        //STEP2 - here we give to every processor the ids that are needed by someone
        rDataCommunicator.Broadcast(all_remote_ids,master_rank);

        //STEP3 - here we obtain the list of gps we own and we send it back to the master_rank
        //gather results on master_rank
        for(int i=0; i<world_size; ++i)
        {

            if(i != master_rank)
            {
                if(current_rank == master_rank)
                {
                    std::unordered_map< int, GPType > recv_gps;
                    rDataCommunicator.Recv(recv_gps, i);

                    for(auto& it : recv_gps)
                        all_non_local_gp_map.emplace(it.first, it.second);
                }
                else if(current_rank == i)
                {
                    auto non_local_gp_map = ComputeGpMap(container, all_remote_ids, rDataCommunicator);
                    rDataCommunicator.Send(non_local_gp_map,master_rank);
                }
            }
            else
            {
                auto recv_gps = ComputeGpMap(container, all_remote_ids, rDataCommunicator);

                for(auto& it : recv_gps)
                    all_non_local_gp_map.emplace(it.first, it.second);
            }
        }

        //STEP4 - here we obtain from the master_rank the list of gps we need
        //extract data and send to everyone
        for(int i=0; i<world_size; ++i)
        {
            if(i != master_rank)
            {
                if(current_rank == master_rank) //only master executes
                {
                    auto gp_list = ExtractById(all_non_local_gp_map,collected_remote_ids[i]);

                    //TODO: here we could use separately send and recv
                    rDataCommunicator.Send(gp_list,i);
                }
                else if(current_rank == i) //only processor i executes
                {
                    std::unordered_map< int, GPType > gp_list;
                    rDataCommunicator.Recv(gp_list, master_rank);

                    for(auto& it : gp_list)
                        global_pointers_list.emplace(it.first, it.second);
                }

            }
            else
            {
                auto gp_list = ExtractById(all_non_local_gp_map,collected_remote_ids[i]);

                for(auto& it : gp_list)
                    global_pointers_list.emplace(it.first, it.second);
            }
        }

        return global_pointers_list;

    }

    template< class TContainerType >
    static GlobalPointersVector< typename TContainerType::value_type > RetrieveGlobalIndexedPointers(
        const TContainerType& container,
        const std::vector<int>& id_list,
        const DataCommunicator& rDataCommunicator
    )
    {
        auto global_pointers_list = RetrieveGlobalIndexedPointersMap(container, id_list, rDataCommunicator);

        int current_rank = rDataCommunicator.Rank();

        //compute final array
        GlobalPointersVector< typename TContainerType::value_type > result;
        result.reserve(id_list.size());
        for(unsigned int i=0; i<id_list.size(); ++i)
        {
            auto it = global_pointers_list.find(id_list[i]);
            if(it != global_pointers_list.end())
                result.push_back( it->second );
            else
                KRATOS_ERROR << "id " << id_list[i] << " not found for processor " << current_rank << std::endl;
        }

        return result;

    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "GlobalPointerUtilities" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GlobalPointerUtilities";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

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

    static bool ObjectIsLocal(const Element& elem, const int CurrentRank)
    {
        return true; //if the iterator was found, then it is local!
    }

    static bool ObjectIsLocal(const Condition& cond, const int CurrentRank)
    {
        return true; //if the iterator was found, then it is local!
    }

    //particularizing to the case of nodes
    static bool ObjectIsLocal(const Node<3>& node, const int CurrentRank)
    {
        return node.FastGetSolutionStepValue(PARTITION_INDEX) == CurrentRank;
    }

    ///@}
    ///@name Private Operators
    ///@{


    template< class GPType >
    static std::unordered_map< int, GPType > ExtractById(
        std::unordered_map< int, GPType >& gp_list,
        const std::vector<int>& ids)
    {
        std::unordered_map< int, GPType > extracted_list;
        for(auto id : ids)
        {
            auto gp = gp_list[id];
            extracted_list[id] = gp;
        }
        return extracted_list;
    }

    template< class TContainerType >
    static std::unordered_map< int, GlobalPointer<typename TContainerType::value_type> > ComputeGpMap(
        const TContainerType& container,
        const std::vector<int>& ids,
        const DataCommunicator& rDataCommunicator)
    {
        const int current_rank = rDataCommunicator.Rank();
        std::unordered_map< int, GlobalPointer<typename TContainerType::value_type> > extracted_list;
        
        if(rDataCommunicator.IsDistributed()) {
            // If the execution is distributed, for every entity id find if its in the container, and if it is, if it is local
            // to our partition.
            for(auto id : ids) {
                const auto it = container.find(id);

                if( it != container.end()) { 
                    // Found locally
                    if(ObjectIsLocal(*it, current_rank)){
                        extracted_list.emplace(id, GlobalPointer<typename TContainerType::value_type>(&*it, current_rank));
                    }
                }
            }
        } else {
            // If the execution is not distributed, only check if the id is in the container.
            for(auto id : ids) {
                const auto it = container.find(id);

                if( it != container.end()) {
                    // Found locally                
                    extracted_list.emplace(id, GlobalPointer<typename TContainerType::value_type>(&*it, current_rank));
                }
            }
        }
        return extracted_list;
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
    GlobalPointerUtilities& operator=(GlobalPointerUtilities const& rOther) = delete;

    /// Copy constructor.
    GlobalPointerUtilities(GlobalPointerUtilities const& rOther) = delete;

    ///@}

}; // Class GlobalPointerUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  GlobalPointerUtilities& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GlobalPointerUtilities& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GLOBAL_POINTER_UTILITIES_H_INCLUDED  defined


