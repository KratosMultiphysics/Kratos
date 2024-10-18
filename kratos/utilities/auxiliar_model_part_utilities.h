//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#pragma once

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class AuxiliarModelPartUtilities
 * @ingroup KratosCore
 * @brief This utility includes auxiliar methods not included in the model part to avoid increase more than necessary the API
 * @todo Typo, Auxiliar is not English, it is Auxiliary, please replace it.
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) AuxiliarModelPartUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// The index type definition
    using IndexType = std::size_t;

    /// Counted pointer of AuxiliarModelPartUtilities
    KRATOS_CLASS_POINTER_DEFINITION( AuxiliarModelPartUtilities );

    using DataLocation = Globals::DataLocation;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * The default constructor
     */
    AuxiliarModelPartUtilities(ModelPart& rModelPart):
        mrModelPart(rModelPart)
    {
    }

    virtual ~AuxiliarModelPartUtilities()= default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method adds the given element and the belonging nodes
     * @param pNewElement The new element added
     */
    void AddElementWithNodes(Element::Pointer pNewElement);

    /**
     * @brief Inserts a list of elements and the belonging nodes to a submodelpart provided their Id. Does nothing if applied to the top model part
     * @param rElementIds The ids of the elements
     */
    void AddElementsWithNodes(const std::vector<IndexType>& rElementIds);

    /**
     * @brief Inserts a list of pointers to elements and the belonging nodes
     * @param ItElementsBegin The begin iterator
     * @param ItElementsEnd The end iterator
     * @tparam TIteratorType The class of iterator considered
     */
    template<class TIteratorType >
    void AddElementsWithNodes(
        TIteratorType ItElementsBegin,
        TIteratorType ItElementsEnd
        )
    {
        KRATOS_TRY

        // Using auxiliay method
        ModelPart* p_root_model_part = &mrModelPart.GetRootModelPart();
        std::vector<IndexType> list_of_nodes;
        ModelPart::ElementsContainerType new_elements_to_add ;
        AuxiliaryAddEntitiesWithNodes<ModelPart::ElementsContainerType, TIteratorType>(p_root_model_part->Elements(), new_elements_to_add , list_of_nodes, ItElementsBegin, ItElementsEnd);

        // Add to all of the leaves
        ModelPart* p_current_part = &mrModelPart;
        while(p_current_part->IsSubModelPart()) {
            for(auto it_elem = new_elements_to_add.begin(); it_elem!=new_elements_to_add.end(); ++it_elem) {
                p_current_part->Elements().push_back( *(it_elem.base()) );
            }
            p_current_part->AddNodes(list_of_nodes);

            p_current_part->Elements().Unique();
            p_current_part = &(p_current_part->GetParentModelPart());
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This method adds the given condition and the belonging nodes
     * @param pNewCondition The new condition added
     */
    void AddConditionWithNodes(Condition::Pointer pNewCondition);

    /**
     * @brief Inserts a list of conditions and the belonging nodes to a submodelpart provided their Id. Does nothing if applied to the top model part
     * @param rConditionIds The ids of the conditions
     */
    void AddConditionsWithNodes(const std::vector<IndexType>& rConditionIds);

    /**
     * @brief Inserts a list of pointers to conditions and the belonging nodes
     * @param ItConditionsBegin The begin iterator
     * @param ItConditionsEnd The end iterator
     * @tparam TIteratorType The class of iterator considered
     */
    template<class TIteratorType >
    void AddConditionsWithNodes(
        TIteratorType ItConditionsBegin,
        TIteratorType ItConditionsEnd
        )
    {
        KRATOS_TRY

        // Using auxiliay method
        ModelPart* p_root_model_part = &mrModelPart.GetRootModelPart();
        std::vector<IndexType> list_of_nodes;
        ModelPart::ConditionsContainerType new_conditions_to_add ;
        AuxiliaryAddEntitiesWithNodes<ModelPart::ConditionsContainerType, TIteratorType>(p_root_model_part->Conditions(), new_conditions_to_add, list_of_nodes, ItConditionsBegin, ItConditionsEnd);

        // Add to all of the leaves
        ModelPart* p_current_part = &mrModelPart;
        while(p_current_part->IsSubModelPart()) {
            for(auto it_cond = new_conditions_to_add.begin(); it_cond!=new_conditions_to_add.end(); ++it_cond) {
                p_current_part->Conditions().push_back( *(it_cond.base()) );
            }
            p_current_part->AddNodes(list_of_nodes);

            p_current_part->Conditions().Unique();
            p_current_part = &(p_current_part->GetParentModelPart());
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This method copies the structure of submodelparts
     * @param rModelPartToCopyFromIt The model part to copy from it
     * @param rModelPartToCopyIntoIt The model part where to copy the structure of the submodelparts
     */
    static void CopySubModelPartStructure(const ModelPart& rModelPartToCopyFromIt, ModelPart& rModelPartToCopyIntoIt);

    /**
     * @brief This method ensured that the properties of elements and conditions are on the model part (it does recursively in all model parts)
     * @param RemovePreviousProperties If we clear previous properties and ensure only the properties existing in the elements and conditions (true by default)
     */
    void RecursiveEnsureModelPartOwnsProperties(const bool RemovePreviousProperties = true);

    /**
     * @brief This method ensured that the properties of elements and conditions are on the model part
     * @param RemovePreviousProperties If we clear previous properties and ensure only the properties existing in the elements and conditions (true by default)
     */
    void EnsureModelPartOwnsProperties(const bool RemovePreviousProperties = true);

    /**
     * @brief Remove the element with given Id from mesh with ThisIndex in this modelpart and all its subs.
     * @details This method removes belonging entities too.
     * This means that if the element has nodes defining a condition, and the nodes defining that condition are removed the condition is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param ElementId The id of the element to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongings(IndexType ElementId, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove given element from mesh with ThisIndex in this modelpart and all its subs.
     * @details This method removes belonging entities too.
     * This means that if the element has nodes defining a condition, and the nodes defining that condition are removed the condition is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param rThisElement The reference of the element
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongings(Element& rThisElement, const Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove given element from mesh with ThisIndex in this modelpart and all its subs.
     * @details This method removes belonging entities too.
     * This means that if the element has nodes defining a condition, and the nodes defining that condition are removed the condition is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param pThisElement The pointer to the element to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongings(Element::Pointer pThisElement, const Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove the element with given Id from mesh with ThisIndex in parents, itself and children.
     * @details This method removes belonging entities too.
     * This means that if the element has nodes defining a condition, and the nodes defining that condition are removed the condition is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param ElementId The id of the element to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongingsFromAllLevels(IndexType ElementId, const Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove given element from mesh with ThisIndex in parents, itself and children.
     * @details This method removes belonging entities too.
     * This means that if the element has nodes defining a condition, and the nodes defining that condition are removed the condition is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param rThisElement The reference of the element
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongingsFromAllLevels(Element& rThisElement, const Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove given element from mesh with ThisIndex in parents, itself and children.
     * @details This method removes belonging entities too.
     * This means that if the element has nodes defining a condition, and the nodes defining that condition are removed the condition is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param pThisElement The pointer to the element to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongingsFromAllLevels(Element::Pointer pThisElement, const Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief  It erases all elements identified by "IdentifierFlag" by removing the pointer.
     * @details This method removes belonging entities too.
     * This means that if the element has nodes defining a condition, and the nodes defining that condition are removed the condition is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * Pointers are erased from this level downwards nodes will be automatically destructured when no pointer is left to them
     * @param IdentifierFlag The flag that identifies the entities to remove
     */
    void RemoveElementsAndBelongings(Flags IdentifierFlag = TO_ERASE);

    /**
     * @brief It erases all elements identified by "IdentifierFlag" by removing the pointer.
     * @details This method removes belonging entities too.
     * This means that if the element has nodes defining a condition, and the nodes defining that condition are removed the condition is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * Pointers are erased from this level downwards nodes will be automatically destructured when no pointer is left to them
     * @param IdentifierFlag The flag that identifies the entities to remove
     */
    void RemoveElementsAndBelongingsFromAllLevels(const Flags IdentifierFlag = TO_ERASE);

    /**
     * @brief Remove the condition with given Id from mesh with ThisIndex in this modelpart and all its subs.
     * @details This method removes belonging entities too.
     * This means that if the condition has nodes defining an element, and the nodes defining that element are removed the element is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param ConditionId The ID of the condition to remove
     * @param IdentifierFlag The flag that identifies the entities to removeentities too
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveConditionAndBelongings(IndexType ConditionId, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove given condition from mesh with ThisIndex in this modelpart and all its subs.
     * @details This method removes belonging entities too.
     * This means that if the condition has nodes defining an element, and the nodes defining that element are removed the element is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param rThisCondition The reference to the condition to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveConditionAndBelongings(Condition& ThisCondition, const Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove given condition from mesh with ThisIndex in this modelpart and all its subs. This method removes belonging entities too
     * @details The method check that when removing the nodes doesn't affect to other entities
     * @param pThisCondition The pointer to the condition to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveConditionAndBelongings(Condition::Pointer pThisCondition, const Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove the condition with given Id from mesh with ThisIndex in parents, itself and children.
     * @details This method removes belonging entities too.
     * This means that if the condition has nodes defining an element, and the nodes defining that element are removed the element is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param ConditionId The ID of the condition to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveConditionAndBelongingsFromAllLevels(IndexType ConditionId, const Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove given condition from mesh with ThisIndex in parents, itself and children.
     * @details This method removes belonging entities too.
     * This means that if the condition has nodes defining an element, and the nodes defining that element are removed the element is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param rThisCondition The reference to the condition to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveConditionAndBelongingsFromAllLevels(Condition& rThisCondition, const Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove given condition from mesh with ThisIndex in parents, itself and children.
     * @details This method removes belonging entities too.
     * This means that if the condition has nodes defining an element, and the nodes defining that element are removed the element is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param pThisCondition The pointer to the condition to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveConditionAndBelongingsFromAllLevels(Condition::Pointer pThisCondition, const Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief It erases all conditions identified by "IdentifierFlag" by removing the pointer.
     * @details This method removes belonging entities too.
     * This means that if the condition has nodes defining an element, and the nodes defining that element are removed the element is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * Pointers are erased from this level downwards nodes will be automatically destructured when no pointer is left to them
     * @param IdentifierFlag The flag that identifies the entities to remove
     */
    void RemoveConditionsAndBelongings(Flags IdentifierFlag = TO_ERASE);

    /**
     * @brief It erases all conditions identified by "IdentifierFlag" by removing the pointer.
     * @details This method removes belonging entities too.
     * This means that if the condition has nodes defining an element, and the nodes defining that element are removed the element is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * Pointers are erased from this level downwards nodes will be automatically destructured when no pointer is left to them
     * @param IdentifierFlag The flag that identifies the entities to remove
     */
    void RemoveConditionsAndBelongingsFromAllLevels(const Flags IdentifierFlag = TO_ERASE);

    /**
     * @brief This method removed nodes from submodelparts not contained neither in the elements or conditions
     */
    void RemoveOrphanNodesFromSubModelParts();

    /// To Export a Scalar data (Double/int/...)
    template<class TContainerType>
    void GetScalarData(
        const Variable<typename TContainerType::value_type>& rVariable,
        const DataLocation DataLoc,
        TContainerType& data) const
    {
        KRATOS_TRY

        switch (DataLoc)
        {
        case (DataLocation::NodeHistorical):{
            data.resize(mrModelPart.NumberOfNodes());

            auto inodebegin = mrModelPart.NodesBegin();

            IndexPartition<IndexType>(mrModelPart.NumberOfNodes()).for_each([&](IndexType Index){
                auto inode = inodebegin + Index;

                data[Index] = inode->FastGetSolutionStepValue(rVariable);
            });

            break;
        }
        case (DataLocation::NodeNonHistorical):{
            data.resize(mrModelPart.NumberOfNodes());

            GetScalarDataFromContainer(mrModelPart.Nodes(), rVariable, data);
            break;
        }
        case (DataLocation::Element):{
            data.resize(mrModelPart.NumberOfElements());

            GetScalarDataFromContainer(mrModelPart.Elements(), rVariable, data);
            break;
        }
        case (DataLocation::Condition):{
            data.resize(mrModelPart.NumberOfConditions());

            GetScalarDataFromContainer(mrModelPart.Conditions(), rVariable, data);
            break;
        }
        case (DataLocation::ModelPart):{
            data.resize(1);
            data[0] = mrModelPart[rVariable];
            break;
        }
        case (DataLocation::ProcessInfo):{
            data.resize(1);
            data[0] = mrModelPart.GetProcessInfo()[rVariable];
            break;
        }
        default:{
            KRATOS_ERROR << "unknown Datalocation" << std::endl;
            break;
        }
        }

        KRATOS_CATCH("")
    }

    /// To Export a Vector data (std::vector/array/..)
    template<class TContainerType, class TVarType>
    void GetVectorData(
        const Variable<TVarType>& rVariable,
        const DataLocation DataLoc,
        TContainerType& data) const
    {
        KRATOS_TRY

        switch (DataLoc)
        {
        case (DataLocation::NodeHistorical):{
            unsigned int TSize = mrModelPart.NumberOfNodes() > 0 ? mrModelPart.NodesBegin()->FastGetSolutionStepValue(rVariable).size() : 0;

            TSize = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(TSize);
            data.resize(mrModelPart.NumberOfNodes()*TSize);

            auto inodebegin = mrModelPart.NodesBegin();

            IndexPartition<IndexType>(mrModelPart.NumberOfNodes()).for_each([&](IndexType Index){
                auto inode = inodebegin + Index;

                const auto& r_val = inode->FastGetSolutionStepValue(rVariable);
                for(std::size_t dim = 0 ; dim < TSize ; dim++){
                    data[(Index*TSize) + dim] = r_val[dim];
                }
            });

            break;
        }
        case (DataLocation::NodeNonHistorical):{
            unsigned int TSize = mrModelPart.NumberOfNodes() > 0 ? mrModelPart.NodesBegin()->GetValue(rVariable).size() : 0;

            TSize = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(TSize);

            data.resize(mrModelPart.NumberOfNodes()*TSize);

            GetVectorDataFromContainer(mrModelPart.Nodes(), TSize, rVariable, data);
            break;
        }
        case (DataLocation::Element):{
            unsigned int TSize = mrModelPart.NumberOfElements() > 0 ? mrModelPart.ElementsBegin()->GetValue(rVariable).size() : 0;

            TSize = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(TSize);

            data.resize(mrModelPart.NumberOfElements()*TSize);

            GetVectorDataFromContainer(mrModelPart.Elements(), TSize, rVariable, data);
            break;
        }
        case (DataLocation::Condition):{
            unsigned int TSize = mrModelPart.NumberOfConditions() > 0 ? mrModelPart.ConditionsBegin()->GetValue(rVariable).size() : 0;

            TSize = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(TSize);

            data.resize(mrModelPart.NumberOfConditions()*TSize);

            GetVectorDataFromContainer(mrModelPart.Conditions(), TSize, rVariable, data);
            break;
        }
        case (DataLocation::ModelPart):{
            std::size_t TSize = mrModelPart[rVariable].size();
            data.resize(TSize);

            IndexType counter = 0;
            auto& r_val = mrModelPart[rVariable];
            for(std::size_t dim = 0 ; dim < TSize ; dim++){
                    data[counter++] = r_val[dim];
                }
            break;
        }
        case (DataLocation::ProcessInfo):{
            const std::size_t TSize = mrModelPart.GetProcessInfo()[rVariable].size();
            data.resize(TSize);

            IndexType counter = 0;
            auto& r_val = mrModelPart.GetProcessInfo()[rVariable];
            for(std::size_t dim = 0 ; dim < TSize ; dim++){
                    data[counter++] = r_val[dim];
                }
            break;
        }
        default:{
            KRATOS_ERROR << "unknown Datalocation" << std::endl;
            break;
        }
        }

        KRATOS_CATCH("")
    }

    /// To Import a Scalar data (Double/int/...)
    template<class TContainerType>
    void SetScalarData(
        const Variable<typename TContainerType::value_type>& rVariable,
        const DataLocation DataLoc,
        const TContainerType& rData)
    {
        KRATOS_TRY

        switch (DataLoc)
        {
        case (DataLocation::NodeHistorical):{
            auto inodebegin = mrModelPart.NodesBegin();
            IndexPartition<IndexType>(mrModelPart.NumberOfNodes()).for_each([&](IndexType Index){
                auto inode = inodebegin + Index;

                auto& r_val = inode->FastGetSolutionStepValue(rVariable);
                r_val = rData[Index];
            });

            break;
        }
        case (DataLocation::NodeNonHistorical):{
            SetScalarDataFromContainer(mrModelPart.Nodes(), rVariable, rData);
            break;
        }
        case (DataLocation::Element):{
            SetScalarDataFromContainer(mrModelPart.Elements(), rVariable, rData);
            break;
        }
        case (DataLocation::Condition):{
            SetScalarDataFromContainer(mrModelPart.Conditions(), rVariable, rData);
            break;
        }
        case (DataLocation::ModelPart):{
            mrModelPart[rVariable]= rData[0];
            break;
        }
        case (DataLocation::ProcessInfo):{
            mrModelPart.GetProcessInfo()[rVariable] = rData[0] ;
            break;
        }
        default:{
            KRATOS_ERROR << "unknown Datalocation" << std::endl;
            break;
        }
        }

        KRATOS_CATCH("")
    }

    /// To Import a Vector data (std::vector/array/..)
    template<class TContainerType, class TVarType>
    void SetVectorData(
        const Variable<TVarType>& rVariable,
        const DataLocation DataLoc,
        const TContainerType& rData)
    {
        KRATOS_TRY

        switch (DataLoc)
        {
        case (DataLocation::NodeHistorical):{
            unsigned int size = mrModelPart.NumberOfNodes() > 0 ? mrModelPart.NodesBegin()->FastGetSolutionStepValue(rVariable).size() : 0;

            size = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(size);

            auto inodebegin = mrModelPart.NodesBegin();
            IndexPartition<IndexType>(mrModelPart.NumberOfNodes()).for_each([&](IndexType Index){
                auto inode = inodebegin + Index;
                auto& r_val = inode->FastGetSolutionStepValue(rVariable);

                KRATOS_DEBUG_ERROR_IF(r_val.size() != size) << "mismatch in size!" << std::endl;

                for(std::size_t dim = 0 ; dim < size ; dim++){
                    r_val[dim] = rData[(Index*size) + dim];
                }
            });

            break;
        }
        case (DataLocation::NodeNonHistorical):{
            unsigned int size = mrModelPart.NumberOfNodes() > 0 ? mrModelPart.NodesBegin()->GetValue(rVariable).size() : 0;

            size = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(size);

            SetVectorDataFromContainer(mrModelPart.Nodes(), size, rVariable, rData);
            break;
        }
        case (DataLocation::Element):{
            unsigned int size = mrModelPart.NumberOfElements() > 0 ? mrModelPart.ElementsBegin()->GetValue(rVariable).size() : 0;

            size = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(size);

            SetVectorDataFromContainer(mrModelPart.Elements(), size, rVariable, rData);
            break;
        }
        case (DataLocation::Condition):{
            unsigned int size = mrModelPart.NumberOfConditions() > 0 ? mrModelPart.ConditionsBegin()->GetValue(rVariable).size() : 0;

            size = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(size);

            SetVectorDataFromContainer(mrModelPart.Conditions(), size, rVariable, rData);
            break;
        }
        case (DataLocation::ModelPart):{
            const std::size_t size = mrModelPart[rVariable].size();

            IndexType counter = 0;
            auto& r_val = mrModelPart[rVariable];
                for(std::size_t dim = 0 ; dim < size ; dim++){
                    r_val[dim] = rData[counter++];
                }
            break;
            }
        case (DataLocation::ProcessInfo):{
            const std::size_t size = mrModelPart.GetProcessInfo()[rVariable].size();

            IndexType counter = 0;
            auto& r_val = mrModelPart.GetProcessInfo()[rVariable];
            for(std::size_t dim = 0 ; dim < size ; dim++){
                    r_val[dim] = rData[counter++];
                }
            break;
        }
        default:{
            KRATOS_ERROR << "unknown Datalocation" << std::endl;
            break;
        }

        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This method deep copies a whole model part
     * @details When a pointer to Model is provided the provided Model will be considered for the copy, otherwise the Model of the current ModelPart will be considered. The last is the default behaviour.
     * This is deep copy, meaning that every entity is deep copied, so created from scratch. The only thing that would be equal will be the Model if not custom Model is provided
     * @param rNewModelPartName The name of the new model part
     * @param pModel The pointer to the Model that will host the new ModelPart, if nullptr, the current Model will be used.
     * @return The deep copied model part
     */
    ModelPart& DeepCopyModelPart(
        const std::string& rNewModelPartName,
        Model* pModel = nullptr
        );

    /**
     * @brief This method deep copies a entities
     * @details Only works with Element and Condition due to the lack of consistency of the entities Clone methods
     * @param rModelPart The model part to copy the entities
     * @param rEntities The entities to be copied
     * @param rReferenceEntities The entities to be copied
     * @param rGeometryPointerDatabase The database of geometries
     * @tparam TClassContainer rEntities type
     * @tparam TReferenceClassContainer rReferenceEntities type
     */
    template<class TClassContainer, class TReferenceClassContainer>
    void DeepCopyEntities(
        ModelPart& rModelPart,
        TClassContainer& rEntities,
        TReferenceClassContainer& rReferenceEntities,
        std::unordered_map<Geometry<Node>::Pointer,Geometry<Node>::Pointer>& rGeometryPointerDatabase
        )
    {
        KRATOS_TRY

        auto& r_properties= rModelPart.rProperties();
        rEntities.SetMaxBufferSize(rReferenceEntities.GetMaxBufferSize());
        rEntities.SetSortedPartSize(rReferenceEntities.GetSortedPartSize());
        const auto& r_reference_entities_container = rReferenceEntities.GetContainer();
        auto& r_entities_container = rEntities.GetContainer();
        const IndexType number_entities = r_reference_entities_container.size();
        r_entities_container.resize(number_entities);
        const auto it_ent_begin = r_reference_entities_container.begin();
        IndexPartition<std::size_t>(number_entities).for_each([&it_ent_begin,&r_entities_container,&rGeometryPointerDatabase,&r_properties](std::size_t i) {
            auto it_ent = it_ent_begin + i;
            auto& p_old_ent = (*it_ent);
            auto p_new_ent = p_old_ent->Create(p_old_ent->Id(), rGeometryPointerDatabase[p_old_ent->pGetGeometry()], r_properties(p_old_ent->pGetProperties()->Id()));
            p_new_ent->SetData(p_old_ent->GetData());
            p_new_ent->Set(Flags(*p_old_ent));
            r_entities_container[i] = p_new_ent;
        });

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Input and output
    ///@{
  
    /**
    * @brief Retrieve the IDs of neighboring elements for each element.
    * @details This function retrieves the IDs of neighboring elements for each element in the model part.
    * The IDs are stored in an unordered map where the key is the ID of the element and the value
    * is a vector containing the IDs of its neighboring elements.
    * @return An unordered map containing the IDs of neighboring elements for each element.
    */
    std::unordered_map<IndexType, std::vector<IndexType>> RetrieveElementsNeighbourElementsIds();

    /**
    * @brief Retrieve the IDs of neighboring conditions for each condition.
    * @details This function retrieves the IDs of neighboring conditions for each condition in the model part.
    * The IDs are stored in an unordered map where the key is the ID of the condition and the value
    * is a vector containing the IDs of its neighboring conditions.
    * @return An unordered map containing the IDs of neighboring conditions for each condition.
    */
    std::unordered_map<IndexType, std::vector<IndexType>> RetrieveConditionsNeighbourConditionsIds();

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "AuxiliarModelPartUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info() << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << Info() << std::endl;
    }

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    template<typename TDataType, class TContainerType, class TDataContainerType>
    void GetScalarDataFromContainer(
        const TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        TDataContainerType& data) const
    {
        KRATOS_TRY

        DataSizeCheck(rContainer.size(), data.size());

        IndexPartition<std::size_t>(rContainer.size()).for_each([&](std::size_t index){
            const auto& r_entity = *(rContainer.begin() + index);
            data[index] = r_entity.GetValue(rVariable);
        });

        KRATOS_CATCH("")
    }

    template<typename TDataType, class TContainerType, class TDataContainerType>
    void GetVectorDataFromContainer(
        const TContainerType& rContainer,
        const std::size_t VectorSize,
        const Variable<TDataType>& rVariable,
        TDataContainerType& data) const
    {
        KRATOS_TRY

        DataSizeCheck(rContainer.size()*VectorSize, data.size());

        IndexPartition<std::size_t>(rContainer.size()).for_each([&](std::size_t index){
            const auto& r_entity = *(rContainer.begin() + index);
            const auto& r_val = r_entity.GetValue(rVariable);
            for(std::size_t dim = 0 ; dim < VectorSize ; dim++){
                data[(VectorSize*index) + dim] = r_val[dim];
            }
        });

        KRATOS_CATCH("")
    }

    template<typename TDataType, class TContainerType, class TDataContainerType>
    void SetScalarDataFromContainer(
        TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const TDataContainerType& rData) const
    {
        KRATOS_TRY

        DataSizeCheck(rContainer.size(), rData.size());

        IndexPartition<std::size_t>(rContainer.size()).for_each([&](std::size_t index){
            auto& r_entity = *(rContainer.begin() + index);
            r_entity.SetValue(rVariable, rData[index]);
        });

        KRATOS_CATCH("")
    }

    template<typename TDataType, class TContainerType, class TDataContainerType>
    void SetVectorDataFromContainer(
        TContainerType& rContainer,
        const std::size_t VectorSize,
        const Variable<TDataType>& rVariable,
        const TDataContainerType& rData) const
    {
        KRATOS_TRY

        DataSizeCheck(rContainer.size()*VectorSize, rData.size());

        IndexPartition<std::size_t>(rContainer.size()).for_each([&](std::size_t index){
            auto& r_entity = *(rContainer.begin() + index);
            TDataType aux;
            KRATOS_DEBUG_ERROR_IF(aux.size() != VectorSize) << "mismatch in size!" << std::endl;
            for(std::size_t dim = 0 ; dim < VectorSize ; dim++){
                aux[dim] = rData[(VectorSize*index) + dim];
            }
            r_entity.SetValue(rVariable, aux);
        });

        KRATOS_CATCH("")
    }

    void DataSizeCheck(
        const std::size_t ContainerSize,
        const std::size_t DataSize) const
    {
        KRATOS_ERROR_IF(ContainerSize != DataSize) << "Mismatch in size! Container size: " << ContainerSize << " | Data size: " << DataSize << std::endl;
    }

    /**
     * @brief Inserts a list of entities and the belonging nodes to a submodelpart provided their Id. Does nothing if applied to the top model part
	 * @param rEntitiesContainer The entities to be added
     * @param rEntitiesIds The ids of the entities
     */
    template<class TEntitiesContainer>
    void AuxiliaryAddEntitiesWithNodes(
        TEntitiesContainer& rEntitiesContainer,
        const std::vector<IndexType>& rEntitiesIds
        )
    {
        KRATOS_TRY
        
        // Obtain from the root model part the corresponding list of nodes
        const auto it_ent_end = rEntitiesContainer.end();
        std::unordered_set<IndexType> set_of_node_ids;
        for(IndexType i=0; i<rEntitiesIds.size(); ++i) {
          auto it_ent = rEntitiesContainer.find(rEntitiesIds[i]);
          if(it_ent!=it_ent_end) {
            const auto& r_geom = it_ent->GetGeometry();
            for (IndexType j = 0; j < r_geom.size(); ++j) {
              set_of_node_ids.insert(r_geom[j].Id());
            }
          } else {
            KRATOS_ERROR << "The entity with Id " << rEntitiesIds[i] << " does not exist in the root model part";
          }
        }

        // Adding nodes
        std::vector<IndexType> list_of_nodes;
        list_of_nodes.insert(list_of_nodes.end(), set_of_node_ids.begin(), set_of_node_ids.end());
        mrModelPart.AddNodes(list_of_nodes);

        // Add to all of the leaves
        ModelPart* p_current_part = &mrModelPart;
        while(p_current_part->IsSubModelPart()) {
          p_current_part->AddNodes(list_of_nodes);
          p_current_part = &(p_current_part->GetParentModelPart());
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Inserts a list of pointers to elements and the belonging nodes
     * @param rEntitiesContainer The entities to be added
     * @param ItElementsBegin The begin iterator
     * @param ItElementsEnd The end iterator
     * @tparam TEntitiesContainer The class of entities considered
     * @tparam TIteratorType The class of iterator considered
     */
    template<class TEntitiesContainer, class TIteratorType>
    void AuxiliaryAddEntitiesWithNodes(
        TEntitiesContainer& rEntitiesContainer,
        TEntitiesContainer& rAux,
        std::vector<IndexType>& rListOfNodes,
        TIteratorType ItEntitiesBegin,
        TIteratorType ItEntitiesEnd
        )
    {
        KRATOS_TRY
        
        TEntitiesContainer aux_root;
        std::unordered_set<IndexType> set_of_nodes;
        
        const auto it_ent_end = rEntitiesContainer.end();
        for(TIteratorType it_ent = ItEntitiesBegin; it_ent!=ItEntitiesEnd; ++it_ent) {
            auto it_ent_found = rEntitiesContainer.find(it_ent->Id());
            if(it_ent_found == it_ent_end) { // Entity does not exist in the top model part
                aux_root.push_back( *(it_ent.base()) );
                rAux.push_back( *(it_ent.base()) );
                const auto& r_geom = it_ent->GetGeometry();
                for (IndexType i = 0; i < r_geom.size(); ++i) {
                    set_of_nodes.insert(r_geom[i].Id());
                }
            } else { // If it_ent does exist verify it_ent is the same entity
                if(&(*it_ent_found) != &(*it_ent)) { //check if the pointee coincides
                    KRATOS_ERROR << "Attempting to add a new entity wit_enth Id :" << it_ent_found->Id() << ", unfortunately a (different) entity wit_enth the same Id already exists" << std::endl;
                } else {
                    rAux.push_back( *(it_ent.base()) );
                    const auto& r_geom = it_ent->GetGeometry();
                    for (IndexType i = 0; i < r_geom.size(); ++i) {
                        set_of_nodes.insert(r_geom[i].Id());
                    }
                }
            }
        }

        // Adding nodes
        rListOfNodes.insert(rListOfNodes.end(), set_of_nodes.begin(), set_of_nodes.end());
        mrModelPart.AddNodes(rListOfNodes);

        for(auto it_ent = aux_root.begin(); it_ent!=aux_root.end(); ++it_ent) {
            rEntitiesContainer.push_back( *(it_ent.base()) );
        }
        rEntitiesContainer.Unique();

        KRATOS_CATCH("")
    }
    
    /**
     * @brief This method copies the submodelpart structure from the original model part to the new one.
     * @details This method is called recursively
     * @param rOriginalModelPart The original model part
     * @param rNewModelPart The new model part
     */
    void DeepCopySubModelPart(
        const ModelPart& rOldModelPart,
        ModelPart& rNewModelPart
        );

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{
    ///@}
};// class AuxiliarModelPartUtilities

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
