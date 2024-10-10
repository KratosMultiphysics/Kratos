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
//                   Michael Andre, https://github.com/msandre
//                   Jordi Cotela Dalmau
//                   Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/variable.h"
#include "containers/array_1d.h"
#include "includes/model_part.h"

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
 * @class GatherModelPartUtility
 * @ingroup KratosMPI
 * @brief This function gathers a model part from a given rank to the master rank
 * @details The objective of this class is to gather a model part from a given rank to the master rank
 * @author Riccardo Rossi
 */
namespace Kratos
{
class KRATOS_API(KRATOS_MPI_CORE) GatherModelPartUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Nodes container definition
    using NodesContainerType = ModelPart::NodesContainerType;

    /// Elements container definition
    using ElementsContainerType = ModelPart::ElementsContainerType;

    /// Conditions container definition
    using ConditionsContainerType = ModelPart::ConditionsContainerType;

    /// Pointer definition of GatherModelPartUtility
    KRATOS_CLASS_POINTER_DEFINITION(GatherModelPartUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @details This function is designed to obtain data from "origin_model_part.GetMesh(mesh_id)", copy it to a new model part and make rank "gather_rank" to have a copy of it. Transferred nodes will be treated as ghost on the gather_rank
     * @param GatherRank MPI rank to which the model part is gathered
     * @param rOriginModelPart Model part on which the origin mesh is contained
     * @param rDestinationModelPart Model part to which we gather the data
     */
    GatherModelPartUtility(
        const int GatherRank,
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart
        ) : GatherModelPartUtility(GatherRank, rOriginModelPart, 0, rDestinationModelPart)
    {
        // Nothing to do here
    }

    /**
     * @brief Default constructor.
     * @details This function is designed to obtain data from "origin_model_part.GetMesh(mesh_id)", copy it to a new model part and make rank "gather_rank" to have a copy of it. Transferred nodes will be treated as ghost on the gather_rank
     * @param GatherRank MPI rank to which the model part is gathered
     * @param rOriginModelPart Model part on which the origin mesh is contained
     * @param MeshId Id of the mesh which contains the data
     * @param rDestinationModelPart Model part to which we gather the data
     */
    GatherModelPartUtility(
        const int GatherRank,
        ModelPart& rOriginModelPart,
        const int MeshId,
        ModelPart& rDestinationModelPart
        );

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function gathers the model part on the master rank
     */
    void GatherOnMaster();

    /**
     * @brief The function gathers the model part from the master rank to the other ranks
     * @param rVariable The variable to be gathered
     */
    template <class TDataType>
    void GatherOnMaster(const Variable<TDataType>& rVariable);

    /**
     * @brief The function scatters the model part from the master rank to the other ranks
     * @param rVariable The variable to be scattered
     */
    template <class TDataType>
    void ScatterFromMaster(const Variable<TDataType>& rVariable);

    /**
     * @brief Function to gather entities from other partitions
     * @details This function is intended to gather entities from other partitions. The map indicates the partitions to gather and the vector the entities to gather from each partition. In the current serial case it does nothing. The map represents the following. The key (int) the rank where to bring the entities, and the vector of indexes, the ids of the entities to bring
     * @param rModelPart Model part to gather entities from other partitions
     * @param rNodesToBring Nodes to gather from other partitions. It must be an ordered map to ensure the order of the ranks. This is a synchronous MPI implementation, and order matters until asynchronous is implemented.
     * @param rElementsToBring Elements to gather from other partitions. It must be an ordered map to ensure the order of the ranks. This is a synchronous MPI implementation, and order matters until asynchronous is implemented.
     * @param rConditionsToBring Conditions to gather from other partitions. It must be an ordered map to ensure the order of the ranks. This is a synchronous MPI implementation, and order matters until asynchronous is implemented.
     * @param CallExecuteAfterBringingEntities Call Execute after gathering entities
     * @param EchoLevel The verbosity level
     */
    static void GatherEntitiesFromOtherPartitions(
        ModelPart& rModelPart,
        const std::map<int, std::vector<std::size_t>>& rNodesToBring,
        const std::map<int, std::vector<std::size_t>>& rElementsToBring,
        const std::map<int, std::vector<std::size_t>>& rConditionsToBring,
        const bool CallExecuteAfterBringingEntities = true,
        const int EchoLevel = 0
        );

    /**
     * @brief Function to gather nodes from other partitions
     * @details This function is intended to gather nodes from other partitions. The map indicates the partitions to gather and the vector the nodes to gather from each partition. In the current serial case it does nothing.
     * @param rModelPart Model part to gather nodes from other partitions
     * @param rNodesToBring Nodes to gather from other partitions. It must be an ordered map to ensure the order of the ranks. This is a synchronous MPI implementation, and order matters until asynchronous is implemented. The map represents the following. The key (int) the rank where to bring the nodes, and the vector of indexes, the ids of the nodes to bring
     * @param CallExecuteAfterBringingEntities Call Execute after gathering nodes
     * @param EchoLevel The verbosity level
     */
    static void GatherNodesFromOtherPartitions(
        ModelPart& rModelPart,
        const std::map<int, std::vector<std::size_t>>& rNodesToBring,
        const bool CallExecuteAfterBringingEntities = true,
        const int EchoLevel = 0
        );

    /**
     * @brief Function to gather elements from other partitions
     * @details This function is intended to gather elements from other partitions. The map indicates the partitions to gather and the vector the elements to gather from each partition. In the current serial case it does nothing.
     * @param rModelPart Model part to gather elements from other partitions
     * @param rElementsToBring Elements to gather from other partitions. It must be an ordered map to ensure the order of the ranks. This is a synchronous MPI implementation, and order matters until asynchronous is implemented. The map represents the following. The key (int) the rank where to bring the elements, and the vector of indexes, the ids of the elements to bring
     * @param CallExecuteAfterBringingEntities Call Execute after gathering elements
     * @param EchoLevel The verbosity level
     */
    static void GatherElementsFromOtherPartitions(
        ModelPart& rModelPart,
        const std::map<int, std::vector<std::size_t>>& rElementsToBring,
        const bool CallExecuteAfterBringingEntities = true,
        const int EchoLevel = 0
        );

    /**
     * @brief Function to gather conditions from other partitions
     * @details This function is intended to gather conditions from other partitions. The map indicates the partitions to gather and the vector the conditions to gather from each partition. In the current serial case it does nothing.
     * @param rModelPart Model part to gather entities from other partitions
     * @param rConditionsToBring Conditions to gather from other partitions. It must be an ordered map to ensure the order of the ranks. This is a synchronous MPI implementation, and order matters until asynchronous is implemented. The map represents the following. The key (int) the rank where to bring the conditions, and the vector of indexes, the ids of the conditions to bring
     * @param CallExecuteAfterBringingEntities Call Execute after gathering conditions
     * @param EchoLevel The verbosity level
     */
    static void GatherConditionsFromOtherPartitions(
        ModelPart& rModelPart,
        const std::map<int, std::vector<std::size_t>>& rConditionsToBring,
        const bool CallExecuteAfterBringingEntities = true,
        const int EchoLevel = 0
        );

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
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart; /// The model part to be gathered
    int mGatherRank;        /// The rank to which the model part is gathered
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Function to gather entities from other partitions
     * @details This function is intended to gather entities from other partitions. The map indicates the partitions to gather and the vector the entities to gather from each partition. In the current serial case it does nothing.
     * @param rModelPart Model part to gather entities from other partitions
     * @param rEntitiesToBring Entities to gather from other partitions. It must be an ordered map to ensure the order of the ranks. This is a synchronous MPI implementation, and order matters until asynchronous is implemented.
     * @param EchoLevel The verbosity level
     */
    template <class TObjectType>
    static void GatherEntityFromOtherPartitions(
        ModelPart& rModelPart,
        const std::map<int, std::vector<std::size_t>>& rEntitiesToBring,
        const int EchoLevel = 0
        );

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

    GatherModelPartUtility& operator=(GatherModelPartUtility const& rOther) = delete;

    /// Copy constructor.

    GatherModelPartUtility(GatherModelPartUtility const& rOther) = delete;

    ///@}

}; // Class GatherModelPartUtility

///@}
///@name Type Definitions
///@{

extern template void GatherModelPartUtility::GatherOnMaster(const Variable<double>&);
extern template void GatherModelPartUtility::GatherOnMaster(const Variable<array_1d<double, 3>>&);
extern template void GatherModelPartUtility::ScatterFromMaster(const Variable<double>&);
extern template void GatherModelPartUtility::ScatterFromMaster(const Variable<array_1d<double, 3>>&);

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  GatherModelPartUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const GatherModelPartUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.
