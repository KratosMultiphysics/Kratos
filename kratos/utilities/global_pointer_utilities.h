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
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/node.h"
#include "includes/geometrical_object.h"
#include "includes/variables.h"
#include "includes/data_communicator.h"
#include "includes/global_pointer.h"
#include "containers/global_pointers_vector.h"
#include "containers/global_pointers_unordered_map.h"
#include "utilities/communication_coloring_utilities.h"

namespace Kratos
{
///@addtogroup KratosCore
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

/**
 * @class GlobalPointerUtilities
 * @brief This class is used to manage global pointers. Overall, the GlobalPointerUtilities class provides a useful utility for retrieving global pointers to entities in a distributed Kratos simulation.
 * @details The class provides a method GetGlobalPointers that retrieves a map of global pointers corresponding to the given entity ids. The GetGlobalPointers method takes a reference to a container rContainer that holds the entities, a const reference to a vector of ids rIdList for the entities we need to retrieve, and a reference to a DataCommunicator object rDataCommunicator. The method returns an unordered map of global pointers corresponding to the given entity ids.
 * The method first checks if the execution is distributed, and if so, it finds the rank of the node that owns the entity if it is remote, and gathers everything onto the master_rank processor. It then creates a map of global pointers corresponding to the given entity ids, where the global pointers point to the entities in the container and are tagged with their global rank.
 * @ingroup KratosCore
 * @see GlobalPointer
 * @see GlobalPointersVector
 * @see GlobalPointersUnorderedMap
 * @author Riccardo Rossi
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

    /**
     * @brief Retrieves a map of global pointers corresponding to the given entity ids, where the global pointers point to the entities in the container and are tagged with their global rank.
     * @details If the execution is distributed, it finds the rank of the node that owns the entity if it is remote, and gathers everything onto the master_rank processor.
     * @param rContainer a reference to the container that holds the entities
     * @param rIdList a const reference to a vector of ids for the entities we need to retrieve
     * @param rDataCommunicator a reference to a DataCommunicator object
     * @return the unordered map of global pointers corresponding to the given entity ids in (id, global_pointer) pairs
     */
    template< class TContainerType >
    static std::unordered_map< int, GlobalPointer<typename TContainerType::value_type> > RetrieveGlobalIndexedPointersMap(
        const TContainerType& rContainer,
        const std::vector<int>& rIdList,
        const DataCommunicator& rDataCommunicator
        )
    {
        using GPType = GlobalPointer<typename TContainerType::value_type>;

        std::unordered_map< int, GPType > global_pointers_list;
        const int current_rank = rDataCommunicator.Rank();
        const int world_size = rDataCommunicator.Size();

        std::vector<int> remote_ids;

        if(rDataCommunicator.IsDistributed()) {
            // If the execution is distributed, for every entity id find if its in the container, and if it is, if it is local
            // to our partition.
            for(const int id : rIdList ) {
                const auto it = rContainer.find(id);

                if( it != rContainer.end()) {
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
            for(const int id : rIdList ) {
                const auto it = rContainer.find(id);
                if( it != rContainer.end()) {
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
        for(int i=0; i<world_size; ++i) {
            if(i != master_rank) {
                if(current_rank == master_rank) { //only master executes
                    rDataCommunicator.Recv(collected_remote_ids[i],i);
                } else if(current_rank == i) { //only processor i executes
                    rDataCommunicator.Send(remote_ids,master_rank);
                }
            } else { //no communication needed
                if(current_rank == master_rank) //only master executes
                    collected_remote_ids[i] = remote_ids;
            }

            if(current_rank == master_rank) {
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

        if(current_rank == master_rank) {
            std::sort(all_remote_ids.begin(), all_remote_ids.end());
            auto last = std::unique(all_remote_ids.begin(), all_remote_ids.end());
            all_remote_ids.erase(last, all_remote_ids.end());
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
        for(int i=0; i<world_size; ++i) {
            if(i != master_rank) {
                if(current_rank == master_rank) {
                    std::unordered_map< int, GPType > recv_gps;
                    rDataCommunicator.Recv(recv_gps, i);

                    for(auto& it : recv_gps)
                        all_non_local_gp_map.emplace(it.first, it.second);
                } else if(current_rank == i) {
                    auto non_local_gp_map = ComputeGpMap(rContainer, all_remote_ids, rDataCommunicator);
                    rDataCommunicator.Send(non_local_gp_map,master_rank);
                }
            } else {
                auto recv_gps = ComputeGpMap(rContainer, all_remote_ids, rDataCommunicator);

                for(auto& it : recv_gps)
                    all_non_local_gp_map.emplace(it.first, it.second);
            }
        }

        //STEP4 - here we obtain from the master_rank the list of gps we need
        //extract data and send to everyone
        for(int i=0; i<world_size; ++i) {
            if(i != master_rank) {
                if(current_rank == master_rank) { //only master executes
                    auto gp_list = ExtractById(all_non_local_gp_map,collected_remote_ids[i]);

                    //TODO: here we could use separately send and recv
                    rDataCommunicator.Send(gp_list,i);
                } else if(current_rank == i) { //only processor i executes
                    std::unordered_map< int, GPType > gp_list;
                    rDataCommunicator.Recv(gp_list, master_rank);

                    for(auto& it : gp_list)
                        global_pointers_list.emplace(it.first, it.second);
                }
            } else {
                auto gp_list = ExtractById(all_non_local_gp_map,collected_remote_ids[i]);

                for(auto& it : gp_list)
                    global_pointers_list.emplace(it.first, it.second);
            }
        }

        return global_pointers_list;
    }

    /**
     * @brief Retrieve global pointers for entities in container, given a data communicator. Only local entities are retrieved
     * @param rContainer container of entities from which to retrieve pointers
     * @param rDataCommunicator data communicator to use for retrieval
     * @return vector of global pointers to entities in container
     */
    template< class TContainerType >
    static GlobalPointersVector< typename TContainerType::value_type > LocalRetrieveGlobalPointers(
        const TContainerType& rContainer,
        const DataCommunicator& rDataCommunicator
        )
    {
        // Retrieve the ids
        std::vector<int> local_id_list;
        local_id_list.reserve(rContainer.size());
        for (const auto& r_entity : rContainer) {
            local_id_list.push_back(r_entity.Id());
        }

        // Retrieve the global pointers
        return RetrieveGlobalIndexedPointers(rContainer, local_id_list, rDataCommunicator);
    }

    /**
     * @brief Retrieve global pointers for entities in container, given a data communicator. All entities are retrieved
     * @param rContainer container of entities from which to retrieve pointers
     * @param rDataCommunicator data communicator to use for retrieval
     * @return vector of global pointers to entities in container
     */
    template< class TContainerType >
    static GlobalPointersVector< typename TContainerType::value_type > GlobalRetrieveGlobalPointers(
        const TContainerType& rContainer,
        const DataCommunicator& rDataCommunicator
        )
    {
        // Getting world size
        const int world_size = rDataCommunicator.Size();

        // Getting number of entities
        const std::size_t number_of_entities = rContainer.size();

        // Getting global number of points
        std::vector<int> number_of_entities_per_partition(world_size);
        std::vector<int> send_number_of_entities_per_partition(1, number_of_entities);
        rDataCommunicator.AllGather(send_number_of_entities_per_partition, number_of_entities_per_partition);

        // Retrieve the ids
        std::vector<int> global_id_list, local_id_list;
        local_id_list.reserve(number_of_entities);
        for (const auto& r_entity : rContainer) {
            local_id_list.push_back(r_entity.Id());
        }

        // Generate vectors with sizes for AllGatherv
        std::vector<int> recv_sizes(number_of_entities_per_partition);
        int message_size = 0;
        std::vector<int> recv_offsets(world_size, 0);
        for (int i_rank = 0; i_rank < world_size; i_rank++) {
            recv_offsets[i_rank] = message_size;
            message_size += recv_sizes[i_rank];
        }
        global_id_list.resize(message_size);

        // Invoque AllGatherv
        rDataCommunicator.AllGatherv(local_id_list, global_id_list, recv_sizes, recv_offsets);

        // Retrieve the global pointers
        return RetrieveGlobalIndexedPointers(rContainer, global_id_list, rDataCommunicator);
    }

    /**
     * @brief Retrieve global indexed pointers from container and data communicator
     * @param rContainer Container to retrieve pointers from
     * @param rIdList List of ids to retrieve pointers for
     * @param rDataCommunicator Communicator to retrieve data from
     * @return GlobalPointersVector containing pointers for all ids in rIdList
     * @throws KRATOS_ERROR if an id in rIdList is not found for the current processor rank
     */
    template< class TContainerType >
    static GlobalPointersVector< typename TContainerType::value_type > RetrieveGlobalIndexedPointers(
        const TContainerType& rContainer,
        const std::vector<int>& rIdList,
        const DataCommunicator& rDataCommunicator
    )
    {
        auto global_pointers_list = RetrieveGlobalIndexedPointersMap(rContainer, rIdList, rDataCommunicator);

        const int current_rank = rDataCommunicator.Rank();

        // Compute final array
        GlobalPointersVector< typename TContainerType::value_type > result;
        result.reserve(rIdList.size());
        for(unsigned int i=0; i<rIdList.size(); ++i) {
            auto it = global_pointers_list.find(rIdList[i]);
            if(it != global_pointers_list.end())
                result.push_back( it->second );
            else
                KRATOS_ERROR << "The id " << rIdList[i] << " was not found for processor " << current_rank << std::endl;
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

    /**
     * @brief Determines if an object meets a certain geometrical object and is located locally.
     * @param rGeometricalObject the geometrical object that the object must meet
     * @param CurrentRank the current rank of the object
     * @return true if the object is local, false otherwise
     */
    static bool ObjectIsLocal(const GeometricalObject& rGeometricalObject, const int CurrentRank)
    {
        return true; //if the iterator was found, then it is local!
    }

    /**
     * @brief Determines if an object is local based on its partition index.
     * @param rNode The Node object to check.
     * @param CurrentRank The current partition index to compare with.
     * @return True if the object is local, false otherwise.
     */
    static bool ObjectIsLocal(const Node& rNode, const int CurrentRank)
    {
        return rNode.FastGetSolutionStepValue(PARTITION_INDEX) == CurrentRank;
    }

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief Extracts values from an unordered map by their ids and returns them in a new map.
     * @param rGPList the unordered map to extract values from
     * @param rIds a vector of ids to extract values for
     * @return a new unordered map containing the extracted values
     */
    template< class GPType >
    static std::unordered_map< int, GPType > ExtractById(
        std::unordered_map< int, GPType >& rGPList,
        const std::vector<int>& rIds)
    {
        std::unordered_map< int, GPType > extracted_list;
        for(auto id : rIds){
            auto gp = rGPList[id];
            extracted_list[id] = gp;
        }
        return extracted_list;
    }

    /**
     * @brief Compute a map of global pointers to entities with ids contained in rIds.
     * @param rContainer A reference to a constant TContainerType object containing the entities.
     * @param rIds A vector of integers representing the ids of the entities to extract.
     * @param rDataCommunicator A reference to a constant DataCommunicator object used for distributed execution.
     * @return An unordered map of integers and GlobalPointer objects to entities with ids contained in rIds.
     */
    template< class TContainerType >
    static std::unordered_map< int, GlobalPointer<typename TContainerType::value_type> > ComputeGpMap(
        const TContainerType& rContainer,
        const std::vector<int>& rIds,
        const DataCommunicator& rDataCommunicator)
    {
        const int current_rank = rDataCommunicator.Rank();
        std::unordered_map< int, GlobalPointer<typename TContainerType::value_type> > extracted_list;

        if(rDataCommunicator.IsDistributed()) {
            // If the execution is distributed, for every entity id find if its in the container, and if it is, if it is local
            // to our partition.
            for(auto id : rIds) {
                const auto it = rContainer.find(id);

                if( it != rContainer.end()) {
                    // Found locally
                    if(ObjectIsLocal(*it, current_rank)){
                        extracted_list.emplace(id, GlobalPointer<typename TContainerType::value_type>(&*it, current_rank));
                    }
                }
            }
        } else {
            // If the execution is not distributed, only check if the id is in the container.
            for(auto id : rIds) {
                const auto it = rContainer.find(id);

                if( it != rContainer.end()) {
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