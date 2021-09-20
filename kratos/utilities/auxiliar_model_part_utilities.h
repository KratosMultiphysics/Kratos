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
#include "utilities/parallel_utilities.h"

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

enum class DataLocation {
    NodeHistorical,
    NodeNonHistorical,
    Element,
    Condition,
    ModelPart,
    ProcessInfo
};

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


    /// To Export a Scalar data (Double/int/...)
    template<typename TDataType>
    void GetScalarData(
        const Variable<TDataType>& rVariable,
        const DataLocation DataLoc,
        std::vector<TDataType>& data) const
    {
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
    }

    /// To Export a Vector data (std::vector/array/..)
    template<class TDataType>
    void GetVectorData(
        const Variable<TDataType>& rVariable,
        const DataLocation DataLoc,
        std::vector<double>& data) const
    {
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
    }

    /// To Import a Scalar data (Double/int/...)
    template<typename TDataType>
    void SetScalarData(
        const Variable<TDataType>& rVariable,
        const DataLocation DataLoc,
        const std::vector<TDataType>& rData)
    {
        switch (DataLoc)
        {
        case (DataLocation::NodeHistorical):{
            ImportDataSizeCheck(mrModelPart.NumberOfNodes(), rData.size());

            auto inodebegin = mrModelPart.NodesBegin();
            IndexPartition<IndexType>(mrModelPart.NumberOfNodes()).for_each([&](IndexType Index){
                auto inode = inodebegin + Index;

                auto& r_val = inode->FastGetSolutionStepValue(rVariable);
                r_val = rData[Index];
            });

            break;
        }
        case (DataLocation::NodeNonHistorical):{
            ImportDataSizeCheck(mrModelPart.NumberOfNodes(), rData.size());

            SetScalarDataFromContainer(mrModelPart.Nodes(), rVariable, rData);
            break;
        }
        case (DataLocation::Element):{
            ImportDataSizeCheck(mrModelPart.NumberOfElements(), rData.size());

            SetScalarDataFromContainer(mrModelPart.Elements(), rVariable, rData);
            break;
        }
        case (DataLocation::Condition):{
            ImportDataSizeCheck(mrModelPart.NumberOfConditions(), rData.size());

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

    }

    /// To Import a Vector data (std::vector/array/..)
    template<class TDataType>
    void SetVectorData(
        const Variable<TDataType>& rVariable,
        const DataLocation DataLoc,
        const std::vector<double>& rData)
    {
        switch (DataLoc)
        {
        case (DataLocation::NodeHistorical):{
            unsigned int size = mrModelPart.NumberOfNodes() > 0 ? mrModelPart.NodesBegin()->FastGetSolutionStepValue(rVariable).size() : 0;

            size = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(size);
            ImportDataSizeCheckVector(mrModelPart.NumberOfNodes()*size , rData.size());

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

            ImportDataSizeCheckVector(mrModelPart.NumberOfNodes()*size , rData.size());

            SetVectorDataFromContainer(mrModelPart.Nodes(), size, rVariable, rData);
            break;
        }
        case (DataLocation::Element):{
            unsigned int size = mrModelPart.NumberOfElements() > 0 ? mrModelPart.ElementsBegin()->GetValue(rVariable).size() : 0;

            size = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(size);

            ImportDataSizeCheckVector(mrModelPart.NumberOfElements()*size , rData.size());

            SetVectorDataFromContainer(mrModelPart.Elements(), size, rVariable, rData);
            break;
        }
        case (DataLocation::Condition):{
            unsigned int size = mrModelPart.NumberOfConditions() > 0 ? mrModelPart.ConditionsBegin()->GetValue(rVariable).size() : 0;

            size = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(size);

            ImportDataSizeCheckVector(mrModelPart.NumberOfConditions()*size , rData.size());

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

    }

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

    template<typename TDataType, class TContainerType>
    void GetScalarDataFromContainer(TContainerType& rContainer, const Variable<TDataType>& rVariable, std::vector<TDataType>& data) const
    {
        IndexPartition<std::size_t>(rContainer.size()).for_each([&](std::size_t index){
            auto& r_entity = *(rContainer.begin() + index);
            data[index] = r_entity.GetValue(rVariable);
        });
    }

    template<typename TDataType, class TContainerType>
    void GetVectorDataFromContainer(TContainerType& rContainer, const std::size_t TSize, const Variable<TDataType>& rVariable, std::vector<double>& data) const
    {
        IndexPartition<std::size_t>(rContainer.size()).for_each([&](std::size_t index){
            const auto& r_entity = *(rContainer.begin() + index);
            const auto& r_val = r_entity.GetValue(rVariable);
            for(std::size_t dim = 0 ; dim < TSize ; dim++){
                data[(TSize*index) + dim] = r_val[dim];
            }
        });
    }

    template<typename TDataType, class TContainerType>
    void SetScalarDataFromContainer(TContainerType& rContainer, const Variable<TDataType>& rVariable, const std::vector<TDataType>& rData)
    {
        IndexPartition<std::size_t>(rContainer.size()).for_each([&](std::size_t index){
            auto& r_entity = *(rContainer.begin() + index);
            r_entity.SetValue(rVariable,rData[index]);
        });
    }

    template<typename TDataType, class TContainerType>
    void SetVectorDataFromContainer(TContainerType& rContainer, const std::size_t size, const Variable<TDataType>& rVariable, const std::vector<double>& rData)
    {
        IndexPartition<std::size_t>(rContainer.size()).for_each([&](std::size_t index){
            auto& r_entity = *(rContainer.begin() + index);
            TDataType aux;
            KRATOS_DEBUG_ERROR_IF(aux.size() != size) << "mismatch in size!" << std::endl;
            for(std::size_t dim = 0 ; dim < size ; dim++){
                aux[dim] = rData[(size*index) + dim];
            }
            r_entity.SetValue(rVariable, aux);
        });
    }

    // Only for SetScalarData()
    void ImportDataSizeCheck(std::size_t rContainerSize, std::size_t rSize){
        KRATOS_ERROR_IF(rContainerSize != rSize) << "mismatch in size! Expected size: " << rContainerSize << std::endl;
    }

    // Only for SetVectorData()
    void ImportDataSizeCheckVector(std::size_t rContainerSize, std::size_t rSize){
        KRATOS_ERROR_IF(rContainerSize != rSize) << "mismatch in size! Expected size: " << rContainerSize << std::endl;
    }


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
