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
     * @brief Remove the element with given Id from mesh with ThisIndex in this modelpart and all its subs. This method removes belonging entities too
     * @details The method check that when removing the nodes doesn't affect to other entities 
     * @param ElementId The id of the element to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongings(IndexType ElementId, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /** 
     * @brief Remove given element from mesh with ThisIndex in this modelpart and all its subs. This method removes belonging entities too
     * @details The method check that when removing the nodes doesn't affect to other entities 
     * @param rThisElement The reference of the element
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongings(Element& rThisElement, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /** 
     * @brief Remove given element from mesh with ThisIndex in this modelpart and all its subs. This method removes belonging entities too
     * @details The method check that when removing the nodes doesn't affect to other entities 
     * @param pThisElement The pointer to the element to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongings(Element::Pointer pThisElement, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /** 
     * @brief Remove the element with given Id from mesh with ThisIndex in parents, itself and children. This method removes belonging entities too
     * @details The method check that when removing the nodes doesn't affect to other entities 
     * @param ElementId The id of the element to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongingsFromAllLevels(IndexType ElementId, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /** 
     * @brief Remove given element from mesh with ThisIndex in parents, itself and children. This method removes belonging entities too
     * @details The method check that when removing the nodes doesn't affect to other entities 
     * @param rThisElement The reference of the element
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongingsFromAllLevels(Element& rThisElement, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /** 
     * @brief Remove given element from mesh with ThisIndex in parents, itself and children. This method removes belonging entities too
     * @details The method check that when removing the nodes doesn't affect to other entities 
     * @param pThisElement The pointer to the element to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveElementAndBelongingsFromAllLevels(Element::Pointer pThisElement, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /** 
    * @brief  It erases all elements identified by "IdentifierFlag" by removing the pointer. This method removes belonging entities too
    * @details Pointers are erased from this level downwards nodes will be automatically destructured when no pointer is left to them
    * The method check that when removing the nodes doesn't affect to other entities
    * @param IdentifierFlag The flag that identifies the entities to remove
    */
    void RemoveElementsAndBelongings(Flags IdentifierFlag = TO_ERASE);

    /** 
     * @brief It erases all elements identified by "IdentifierFlag" by removing the pointer. This method removes belonging entities too
     * @details Pointers will be erase from all levels nodes will be automatically destructured when no pointer is left to them
     * The method check that when removing the nodes doesn't affect to other entities
     * @param IdentifierFlag The flag that identifies the entities to remove
     */
    void RemoveElementsAndBelongingsFromAllLevels(Flags IdentifierFlag = TO_ERASE);

    /** 
     * @brief Remove the condition with given Id from mesh with ThisIndex in this modelpart and all its subs. This method removes belonging
     * @details The method check that when removing the nodes doesn't affect to other entities 
     * @param ConditionId The ID of the condition to remove
     * @param IdentifierFlag The flag that identifies the entities to removeentities too
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveConditionAndBelongings(IndexType ConditionId, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /** 
     * @brief Remove given condition from mesh with ThisIndex in this modelpart and all its subs. This method removes belonging entities too
     * @details The method check that when removing the nodes doesn't affect to other entities 
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
     * @brief Remove the condition with given Id from mesh with ThisIndex in parents, itself and children. This method removes belonging entities too
     * @details The method check that when removing the nodes doesn't affect to other entities 
     * @param ConditionId The ID of the condition to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveConditionAndBelongingsFromAllLevels(IndexType ConditionId, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /** 
     * @brief Remove given condition from mesh with ThisIndex in parents, itself and children. This method removes belonging entities too
     * @details The method check that when removing the nodes doesn't affect to other entities 
     * @param rThisCondition The reference to the condition to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveConditionAndBelongingsFromAllLevels(Condition& rThisCondition, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /** 
     * @brief Remove given condition from mesh with ThisIndex in parents, itself and children. This method removes belonging entities too
     * @details The method check that when removing the nodes doesn't affect to other entities 
     * @param pThisCondition The pointer to the condition to remove
     * @param IdentifierFlag The flag that identifies the entities to remove
     * @param ThisIndex The index of the mesh where remove the entity
     */
    void RemoveConditionAndBelongingsFromAllLevels(Condition::Pointer pThisCondition, Flags IdentifierFlag = TO_ERASE, IndexType ThisIndex = 0);

    /** 
     * @brief It erases all conditions identified by "IdentifierFlag" by removing the pointer. This method removes belonging entities too
     * @details Pointers are erased from this level downwards nodes will be automatically destructured when no pointer is left to them
     * The method check that when removing the nodes doesn't affect to other entities 
     * @param IdentifierFlag The flag that identifies the entities to remove
     */
    void RemoveConditionsAndBelongings(Flags IdentifierFlag = TO_ERASE);

    /** 
     * @brief It erases all conditions identified by "IdentifierFlag" by removing the pointer. This method removes belonging entities too
     * @details Pointers will be erase from all levels nodes will be automatically destructured when no pointer is left to them
     * The method check that when removing the nodes doesn't affect to other entities 
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
