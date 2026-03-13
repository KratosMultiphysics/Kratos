//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) OptAppModelPartUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Static operations
    ///@{

    /**
     * @brief Get Newly Created or Existing Model Parts List With Common Reference Entities Between Reference List And Examined List
     *
     * This method returns newly created or gets an existing sub-model part in each model part in rReferenceModelParts with entities
     * common between that model part and all the model parts in rExaminedModelPartsList with following
     * entities.
     *      1. If AreNodesConsidered is made to true, then common nodes are found and added to the resulting
     *         sub-model part.
     *      2. If AreConditionsConsidered is made to true, then common conditions are found and added to the resulting
     *         sub-model part. All the nodes belonging to added conditions are also added to the sub-model part.
     *      3. If AreElementsConsidered is made to true, then common elements are found and added to the resulting
     *         sub-model part. All the nodes belonging to added elements are also added to the sub-model part.
     *      4. If AreNeighboursConsidered is made to true, then first common nodes are found. Then neighbour conditions (if AreConditionsConsidered is
     *         set to true) and elements (if AreElementsConsidered is set to true) for those common nodes are found and
     *         added to the sub-model part. All nodes belonging to added conditions and elements are also added to the sub-model part.
     *
     * The created sub-model part is only populated with the entities from the refrence model part. Hence, to find common
     * conditions and elements, each entity from each model part in rExaminedModelPartsList is searched within each of the
     * model parts in rReferenceModelParts.
     *
     * This method creates the sub-model parts with unique name generated for all the input arguments. If it finds an existing
     * sub-model part with the same name, then it is retrieved. Hence this method only creates these model parts with expensive operations
     * only one. If it is required to force create, then the sub-model parts should be removed. All the sub-model parts created
     * by this method have names starting with "<OPTIMIZATION_APP_AUTO>" so it is easier to search and remove them if required.
     *
     * This method keeps the communicator types, process info, properties and tables as same as the rReferenceModelParts.
     * The local mesh, interface mesh and ghost meshes are properly updated. Hence this method is compatible with OpenMP
     * and MPI.
     *
     * This method does not create or destroy nodes, conditions and elements.
     *
     * @param rExaminedModelPartsList       List of input model parts where the entities are searched between.
     * @param rReferenceModelParts          List of reference model parts where common entities are added to output model parts.
     * @param AreNodesConsidered            If true, common nodes from reference model parts are added to output model parts.
     * @param AreConditionsConsidered       If true, common conditions from reference model parts are added to output model parts.
     * @param AreElementsConsidered         If true, common elements from reference model parts are added to output model parts.
     * @param AreNeighboursConsidered       If true, neighbour conditions and elements from reference model parts for common nodes are added to output model parts.
     * @param EchoLevel                     Echo level for printing info.
     * @return std::vector<ModelPart*>      List of output model parts.
     */
    static std::vector<ModelPart*> GetModelPartsWithCommonReferenceEntities(
        const std::vector<ModelPart*>& rExaminedModelPartsList,
        const std::vector<ModelPart*>& rReferenceModelParts,
        const bool AreNodesConsidered,
        const bool AreConditionsConsidered,
        const bool AreElementsConsidered,
        const bool AreNeighboursConsidered,
        const IndexType EchoLevel = 0);

    /**
     * @brief Removes automatically generated model parts
     *
     * This method removes all the automatically generated sub model parts in the given list of rModelParts
     *
     * @param rModelParts       The list of model parts to look for automatically generated sub model parts.
     */
    static void RemoveModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(
        const std::vector<ModelPart*> rModelParts);

    static void LogModelPartStatus(
        ModelPart& rModelPart,
        const std::string& rStatus);

    static std::vector<std::string> GetModelPartStatusLog(ModelPart& rModelPart);

    static bool CheckModelPartStatus(
        const ModelPart& rModelPart,
        const std::string& rStatus);

    /**
     * @brief Generate model part elements using the geometries of the conditions given in the origin container.
     *
     * @param rOriginContainer          Origin container with conditions.
     * @param rDestinationModelPart     Destination model part.
     * @param rReferenceElement         Reference element to be used for new element creation.
     */
    static void GenerateModelPart(
        ModelPart::ConditionsContainerType& rOriginContainer,
        ModelPart& rDestinationModelPart,
        const Element& rReferenceElement);

    ///@}
};

///@}
} // namespace Kratos