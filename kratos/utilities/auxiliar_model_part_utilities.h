//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_AUXILIAR_MODEL_PART_UTILITIES)
#define KRATOS_AUXILIAR_MODEL_PART_UTILITIES

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "includes/serializer.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The index type definition
    typedef std::size_t IndexType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

/**
 * @class AuxiliarModelPartUtilities
 * @ingroup KratosCore
 * @brief This uility includes auxiliar methods not included in the model part to avoid increase more than necessary the API
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) AuxiliarModelPartUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of AuxiliarModelPartUtilities
    KRATOS_CLASS_POINTER_DEFINITION( AuxiliarModelPartUtilities );

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
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method adds the given element and the belonging nodes
     * @param pNewElement The new element added
     * @param ThisIndex The mesh index
     */
    void AddElementWithNodes(
        Element::Pointer pNewElement,
        IndexType ThisIndex = 0
        );

    /**
     * @brief Inserts a list of elements and the belonging nodes to a submodelpart provided their Id. Does nothing if applied to the top model part
     * @param rElementIds The ids of the elements
     * @param ThisIndex The mesh index
     */
    void AddElementsWithNodes(
        std::vector<IndexType> const& rElementIds,
        IndexType ThisIndex = 0
        );

    /**
     * @brief Inserts a list of pointers to elements and the belonging nodes
     * @param ThisIndex The mesh index
     */
    template<class TIteratorType >
    void AddElementsWit_elemhNodes(
        TIteratorType ItElementsBegin,
        TIteratorType ItElementsEnd,
        IndexType ThisIndex = 0
        )
    {
        KRATOS_TRY
        ModelPart::ElementsContainerType aux;
        ModelPart::ElementsContainerType aux_root;
        ModelPart* p_root_model_part = &mrModelPart.GetRootModelPart();
        std::unordered_set<IndexType> set_of_nodes;

        for(TIteratorType it_elem = ItElementsBegin; it_elem!=ItElementsEnd; ++it_elem) {
            auto it_elem_found = p_root_model_part->Elements().find(it_elem->Id());
            if(it_elem_found == p_root_model_part->ElementsEnd()) { // Element does not exist in the top model part
                aux_root.push_back( *(it_elem.base()) );
                aux.push_back( *(it_elem.base()) );
                const auto& r_geom = it_elem->GetGeometry();
                for (IndexType i = 0; i < r_geom.size(); ++i) {
                    set_of_nodes.insert(r_geom[i].Id());
                }
            } else { // If it_elem does exist verify it_elem is the same element
                if(&(*it_elem_found) != &(*it_elem)) { //check if the pointee coincides
                    KRATOS_ERROR << "Attempting to add a new element wit_elemh Id :" << it_elem_found->Id() << ", unfortunately a (different) element wit_elemh the same Id already exists" << std::endl;
                } else {
                    aux.push_back( *(it_elem.base()) );
                    const auto& r_geom = it_elem->GetGeometry();
                    for (IndexType i = 0; i < r_geom.size(); ++i) {
                        set_of_nodes.insert(r_geom[i].Id());
                    }
                }
            }
        }

        // Adding nodes
        std::vector<IndexType> list_of_nodes;
        list_of_nodes.insert(list_of_nodes.end(), set_of_nodes.begin(), set_of_nodes.end());
        mrModelPart.AddNodes(list_of_nodes);

        for(auto it_elem = aux_root.begin(); it_elem!=aux_root.end(); ++it_elem) {
            p_root_model_part->Elements().push_back( *(it_elem.base()) );
        }
        p_root_model_part->Elements().Unique();

        // Add to all of the leaves
        ModelPart* p_current_part = &mrModelPart;
        while(p_current_part->IsSubModelPart()) {
            for(auto it_elem = aux.begin(); it_elem!=aux.end(); ++it_elem) {
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
     * @param ThisIndex The mesh index
     */
    void AddConditionWithNodes(
        Condition::Pointer pNewCondition,
        IndexType ThisIndex = 0
        );

    /**
     * @brief Inserts a list of conditions and the belonging nodes to a submodelpart provided their Id. Does nothing if applied to the top model part
     * @param rConditionIds The ids of the conditions
     * @param ThisIndex The mesh index
     */
    void AddConditionsWithNodes(
        std::vector<IndexType> const& rConditionIds,
        IndexType ThisIndex = 0
        );

    /**
     * @brief Inserts a list of pointers to conditions and the belonging nodes
     * @param ThisIndex The mesh index
     */
    template<class TIteratorType >
    void AddConditionsWithNodes(
        TIteratorType ItConditionsBegin,
        TIteratorType ItConditionsEnd,
        IndexType ThisIndex = 0
        )
    {
        KRATOS_TRY
        ModelPart::ConditionsContainerType aux;
        ModelPart::ConditionsContainerType aux_root;
        ModelPart* p_root_model_part = &mrModelPart.GetRootModelPart();
        std::unordered_set<IndexType> set_of_nodes;

        for(TIteratorType it_cond = ItConditionsBegin; it_cond!=ItConditionsEnd; ++it_cond) {
            auto it_cond_found = p_root_model_part->Conditions().find(it_cond->Id());
            if(it_cond_found == p_root_model_part->ConditionsEnd()) { // Condition does not exist in the top model part
                aux_root.push_back( *(it_cond.base()) );
                aux.push_back( *(it_cond.base()) );
                const auto& r_geom = it_cond->GetGeometry();
                for (IndexType i = 0; i < r_geom.size(); ++i) {
                    set_of_nodes.insert(r_geom[i].Id());
                }
            } else { // If it_cond does exist verify it_cond is the same condition
                if(&(*it_cond_found) != &(*it_cond)) { //check if the pointee coincides
                    KRATOS_ERROR << "Attempting to add a new condition wit_condh Id :" << it_cond_found->Id() << ", unfortunately a (different) condition wit_condh the same Id already exists" << std::endl;
                } else {
                    aux.push_back( *(it_cond.base()) );
                    const auto& r_geom = it_cond->GetGeometry();
                    for (IndexType i = 0; i < r_geom.size(); ++i) {
                        set_of_nodes.insert(r_geom[i].Id());
                    }
                }
            }
        }

        // Adding nodes
        std::vector<IndexType> list_of_nodes;
        list_of_nodes.insert(list_of_nodes.end(), set_of_nodes.begin(), set_of_nodes.end());
        mrModelPart.AddNodes(list_of_nodes);

        for(auto it_cond = aux_root.begin(); it_cond!=aux_root.end(); ++it_cond) {
            p_root_model_part->Conditions().push_back( *(it_cond.base()) );
        }
        p_root_model_part->Conditions().Unique();

        // Add to all of the leaves
        ModelPart* p_current_part = &mrModelPart;
        while(p_current_part->IsSubModelPart()) {
            for(auto it_cond = aux.begin(); it_cond!=aux.end(); ++it_cond) {
                p_current_part->Conditions().push_back( *(it_cond.base()) );
            }
            p_current_part->AddNodes(list_of_nodes);

            p_current_part->Conditions().Unique();
            p_current_part = &(p_current_part->GetParentModelPart());
        }

        KRATOS_CATCH("")
    }

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
    void RemoveElementAndBelongings(Element& rThisElement, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove given element from mesh with ThisIndex in this modelpart and all its subs.
     * @details This method removes belonging entities too.
     * This means that if the element has nodes defining a condition, and the nodes defining that condition are removed the condition is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param pThisElement The pointer to the element to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongings(Element::Pointer pThisElement, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove the element with given Id from mesh with ThisIndex in parents, itself and children.
     * @details This method removes belonging entities too.
     * This means that if the element has nodes defining a condition, and the nodes defining that condition are removed the condition is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param ElementId The id of the element to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongingsFromAllLevels(IndexType ElementId, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove given element from mesh with ThisIndex in parents, itself and children.
     * @details This method removes belonging entities too.
     * This means that if the element has nodes defining a condition, and the nodes defining that condition are removed the condition is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param rThisElement The reference of the element
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongingsFromAllLevels(Element& rThisElement, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove given element from mesh with ThisIndex in parents, itself and children.
     * @details This method removes belonging entities too.
     * This means that if the element has nodes defining a condition, and the nodes defining that condition are removed the condition is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param pThisElement The pointer to the element to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongingsFromAllLevels(Element::Pointer pThisElement, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

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
    void RemoveElementsAndBelongingsFromAllLevels(Flags IdentifierFlag = TO_ERASE);

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
    void RemoveConditionAndBelongings(Condition& ThisCondition, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove given condition from mesh with ThisIndex in this modelpart and all its subs. This method removes belonging entities too
     * @details The method check that when removing the nodes doesn't affect to other entities
     * @param pThisCondition The pointer to the condition to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveConditionAndBelongings(Condition::Pointer pThisCondition, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove the condition with given Id from mesh with ThisIndex in parents, itself and children.
     * @details This method removes belonging entities too.
     * This means that if the condition has nodes defining an element, and the nodes defining that element are removed the element is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param ConditionId The ID of the condition to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveConditionAndBelongingsFromAllLevels(IndexType ConditionId, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove given condition from mesh with ThisIndex in parents, itself and children.
     * @details This method removes belonging entities too.
     * This means that if the condition has nodes defining an element, and the nodes defining that element are removed the element is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param rThisCondition The reference to the condition to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveConditionAndBelongingsFromAllLevels(Condition& rThisCondition, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /**
     * @brief Remove given condition from mesh with ThisIndex in parents, itself and children.
     * @details This method removes belonging entities too.
     * This means that if the condition has nodes defining an element, and the nodes defining that element are removed the element is removed too
     * Basically the method checks that when removing the nodes doesn't affect to other entities
     * @param pThisCondition The pointer to the condition to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveConditionAndBelongingsFromAllLevels(Condition::Pointer pThisCondition, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

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
    void RemoveConditionsAndBelongingsFromAllLevels(Flags IdentifierFlag = TO_ERASE);

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

    ModelPart& mrModelPart;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
    }

    void load(Serializer& rSerializer)
    {
    }

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
#endif /* KRATOS_AUXILIAR_MODEL_PART_UTILITIES defined */
