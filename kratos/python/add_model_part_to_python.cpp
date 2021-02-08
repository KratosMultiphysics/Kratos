//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//


// System includes

// External includes
#include <pybind11/stl.h>


// Project includes
#include "includes/define_python.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/kratos_components.h"
#include "includes/process_info.h"
#include "utilities/quaternion.h"
#include "python/add_model_part_to_python.h"
#include "python/containers_interface.h"

namespace Kratos
{

namespace Python
{

template<class TDataType>
void AddNodalSolutionStepVariable(ModelPart& rModelPart, Variable<TDataType> const& rThisVariable)
{
    rModelPart.AddNodalSolutionStepVariable(rThisVariable);
}

template<class TDataType>
bool HasNodalSolutionStepVariable(ModelPart& rModelPart, Variable<TDataType> const& rThisVariable)
{
    return rModelPart.HasNodalSolutionStepVariable(rThisVariable);
}

void SetModelPartName(ModelPart& rModelPart, std::string const& NewName)
{
    rModelPart.Name() = NewName;
}

const std::string GetModelPartName(ModelPart const& rModelPart)
{
    return rModelPart.Name();
}

ProcessInfo& GetProcessInfo(ModelPart& rModelPart)
{
    return rModelPart.GetProcessInfo();
}

void SetProcessInfo(ModelPart& rModelPart, ProcessInfo& NewProcessInfo)
{
    rModelPart.SetProcessInfo(NewProcessInfo);
}

ModelPart::MeshType::Pointer ModelPartGetMesh(ModelPart& rModelPart)
{
    return rModelPart.pGetMesh();
}

ModelPart::MeshType::Pointer ModelPartGetMesh2(ModelPart& rModelPart, ModelPart::IndexType MeshIndex)
{
    ModelPart::IndexType number_of_meshes = rModelPart.NumberOfMeshes();
    // adding necessary meshes to the model part.
    ModelPart::MeshType empty_mesh;
    for(ModelPart::IndexType i = number_of_meshes ; i < MeshIndex + 1 ; i++)
        rModelPart.GetMeshes().push_back(Kratos::make_shared<ModelPart::MeshType>(empty_mesh.Clone()));

    return rModelPart.pGetMesh(MeshIndex);
}

Node < 3 > ::Pointer ModelPartCreateNewNode(ModelPart& rModelPart, int Id, double x, double y, double z)
{
    return rModelPart.CreateNewNode(Id, x, y, z);
}

Element::Pointer ModelPartCreateNewElement(ModelPart& rModelPart, const std::string ElementName, ModelPart::IndexType Id, std::vector< ModelPart::IndexType >& NodeIdList, ModelPart::PropertiesType::Pointer pProperties)
{
    if (!KratosComponents<Element>::Has(ElementName)) {
        std::stringstream msg;
        KratosComponents<Element> instance; // creating an instance for using "PrintData"
        instance.PrintData(msg);

        KRATOS_ERROR << "The Element \"" << ElementName << "\" is not registered!\nMaybe you need to import the application where it is defined?\nThe following Elements are registered:\n" << msg.str() << std::endl;
    }

    Geometry< Node < 3 > >::PointsArrayType pElementNodeList;

    for(unsigned int i = 0; i < NodeIdList.size(); i++) {
        pElementNodeList.push_back(rModelPart.pGetNode(NodeIdList[i]));
    }

    return rModelPart.CreateNewElement(ElementName, Id, pElementNodeList, pProperties);
}

Condition::Pointer ModelPartCreateNewCondition(ModelPart& rModelPart, const std::string ConditionName, ModelPart::IndexType Id, std::vector< ModelPart::IndexType >& NodeIdList, ModelPart::PropertiesType::Pointer pProperties)
{
    if (!KratosComponents<Condition>::Has(ConditionName)) {
        std::stringstream msg;
        KratosComponents<Condition> instance; // creating an instance for using "PrintData"
        instance.PrintData(msg);

        KRATOS_ERROR << "The Condition \"" << ConditionName << "\" is not registered!\nMaybe you need to import the application where it is defined?\nThe following Conditions are registered:\n" << msg.str() << std::endl;
    }

    Geometry< Node < 3 > >::PointsArrayType pConditionNodeList;

    for(unsigned int i = 0; i <NodeIdList.size(); i++) {
        pConditionNodeList.push_back(rModelPart.pGetNode(NodeIdList[i]));
    }

    return rModelPart.CreateNewCondition(ConditionName, Id, pConditionNodeList, pProperties);
}


// Nodes

ModelPart::SizeType ModelPartNumberOfNodes1(ModelPart& rModelPart)
{
    return rModelPart.NumberOfNodes();
}

ModelPart::NodesContainerType::Pointer ModelPartGetNodes1(ModelPart& rModelPart)
{
    return rModelPart.pNodes();
}

ModelPart::NodesContainerType::Pointer ModelPartGetNodes2(ModelPart& rModelPart, ModelPart::IndexType ThisIndex)
{
    return rModelPart.pNodes(ThisIndex);
}

void ModelPartSetNodes1(ModelPart& rModelPart, ModelPart::NodesContainerType::Pointer pOtherNodes)
{
    rModelPart.SetNodes(pOtherNodes);
}

void ModelPartSetNodes2(ModelPart& rModelPart, ModelPart::NodesContainerType::Pointer pOtherNodes, ModelPart::IndexType ThisIndex)
{
    rModelPart.SetNodes(pOtherNodes, ThisIndex);
}

ModelPart::NodeType::Pointer ModelPartGetNode1(ModelPart& rModelPart, ModelPart::IndexType NodeId)
{
    return rModelPart.pGetNode(NodeId);
}

ModelPart::NodeType::Pointer ModelPartGetNode2(ModelPart& rModelPart, ModelPart::IndexType NodeId, ModelPart::IndexType ThisIndex)
{
    return rModelPart.pGetNode(NodeId, ThisIndex);
}

void ModelPartRemoveNode1(ModelPart& rModelPart, ModelPart::IndexType NodeId)
{
    rModelPart.RemoveNode(NodeId);
}

void ModelPartRemoveNode2(ModelPart& rModelPart, ModelPart::IndexType NodeId, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemoveNode(NodeId, ThisIndex);
}

void ModelPartRemoveNode3(ModelPart& rModelPart, ModelPart::NodeType::Pointer pThisNode)
{
    rModelPart.RemoveNode(pThisNode);
}

void ModelPartRemoveNode4(ModelPart& rModelPart, ModelPart::NodeType::Pointer pThisNode, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemoveNode(pThisNode, ThisIndex);
}

void ModelPartRemoveNodeFromAllLevels1(ModelPart& rModelPart, ModelPart::IndexType NodeId)
{
    rModelPart.RemoveNodeFromAllLevels(NodeId);
}

void ModelPartRemoveNodeFromAllLevels2(ModelPart& rModelPart, ModelPart::IndexType NodeId, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemoveNodeFromAllLevels(NodeId, ThisIndex);
}

void ModelPartRemoveNodeFromAllLevels3(ModelPart& rModelPart, ModelPart::NodeType::Pointer pThisNode)
{
    rModelPart.RemoveNodeFromAllLevels(pThisNode);
}

void ModelPartRemoveNodeFromAllLevels4(ModelPart& rModelPart, ModelPart::NodeType::Pointer pThisNode, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemoveNodeFromAllLevels(pThisNode, ThisIndex);
}

void ModelPartRemoveNodesFromAllLevels(ModelPart& rModelPart, Flags identifier_flag)
{
    rModelPart.RemoveNodesFromAllLevels(identifier_flag);
}


// Properties

ModelPart::SizeType  ModelPartNumberOfProperties1(ModelPart& rModelPart)
{
    return rModelPart.NumberOfProperties();
}

void ModelPartAddProperties1(ModelPart& rModelPart, Properties::Pointer pNewProperties)
{
    rModelPart.AddProperties(pNewProperties);
}

void ModelPartAddProperties2(ModelPart& rModelPart, Properties::Pointer pNewProperties, ModelPart::IndexType ThisIndex)
{
    rModelPart.AddProperties(pNewProperties, ThisIndex);
}

bool ModelPartHasProperties1(const ModelPart& rModelPart, const unsigned int PropertiesId, const unsigned int MeshId)
{
    return rModelPart.HasProperties(PropertiesId, MeshId);
}

bool ModelPartHasProperties2(const ModelPart& rModelPart, const unsigned int PropertiesId)
{
    return rModelPart.HasProperties(PropertiesId, 0);
}

bool ModelPartHasSubProperties1(const ModelPart& rModelPart, const std::string& rAdress, const unsigned int MeshId)
{
    return rModelPart.HasProperties(rAdress, MeshId);
}

bool ModelPartHasSubProperties2(const ModelPart& rModelPart, const std::string& rAdress)
{
    return rModelPart.HasProperties(rAdress, 0);
}

bool ModelPartRecursivelyHasProperties1(const ModelPart& rModelPart, const unsigned int PropertiesId, const unsigned int MeshId)
{
    return rModelPart.RecursivelyHasProperties(PropertiesId, MeshId);
}

bool ModelPartRecursivelyHasProperties2(const ModelPart& rModelPart, const unsigned int PropertiesId)
{
    return rModelPart.RecursivelyHasProperties(PropertiesId, 0);
}

Properties::Pointer ModelPartCreateNewProperties1(ModelPart& rModelPart, unsigned int PropertiesId, unsigned int MeshId)
{
    return rModelPart.CreateNewProperties(PropertiesId, MeshId);
}

Properties::Pointer ModelPartCreateNewProperties2(ModelPart& rModelPart, unsigned int PropertiesId)
{
    return rModelPart.CreateNewProperties(PropertiesId, 0);
}

Properties::Pointer ModelPartGetPropertiesDirect1(ModelPart& rModelPart, unsigned int PropertiesId, unsigned int MeshId)
{
    return rModelPart.pGetProperties(PropertiesId, MeshId);
}

Properties::Pointer ModelPartGetPropertiesDirect2(ModelPart& rModelPart, unsigned int PropertiesId)
{
    return rModelPart.pGetProperties(PropertiesId);
}

Properties::Pointer ModelPartGetSubProperties1(ModelPart& rModelPart, const std::string& rAdress, unsigned int MeshId)
{
    return rModelPart.pGetProperties(rAdress, MeshId);
}

Properties::Pointer ModelPartSubProperties2(ModelPart& rModelPart, const std::string& rAdress)
{
    return rModelPart.pGetProperties(rAdress);
}

ModelPart::PropertiesContainerType::Pointer ModelPartGetProperties1(ModelPart& rModelPart)
{
    return rModelPart.pProperties();
}

ModelPart::PropertiesContainerType::Pointer ModelPartGetProperties2(ModelPart& rModelPart, ModelPart::IndexType ThisIndex)
{
    return rModelPart.pProperties(ThisIndex);
}

void ModelPartSetProperties1(ModelPart& rModelPart, ModelPart::PropertiesContainerType::Pointer pOtherProperties)
{
    rModelPart.SetProperties(pOtherProperties);
}

void ModelPartSetProperties2(ModelPart& rModelPart, ModelPart::PropertiesContainerType::Pointer pOtherProperties, ModelPart::IndexType ThisIndex)
{
    rModelPart.SetProperties(pOtherProperties, ThisIndex);
}

void ModelPartRemoveProperties1(ModelPart& rModelPart, ModelPart::IndexType PropertiesId)
{
    rModelPart.RemoveProperties(PropertiesId);
}

void ModelPartRemoveProperties2(ModelPart& rModelPart, ModelPart::IndexType PropertiesId, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemoveProperties(PropertiesId, ThisIndex);
}

void ModelPartRemoveProperties3(ModelPart& rModelPart, ModelPart::PropertiesType::Pointer pThisProperties)
{
    rModelPart.RemoveProperties(pThisProperties);
}

void ModelPartRemoveProperties4(ModelPart& rModelPart, ModelPart::PropertiesType::Pointer pThisProperties, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemoveProperties(pThisProperties, ThisIndex);
}


void ModelPartRemovePropertiesFromAllLevels1(ModelPart& rModelPart, ModelPart::IndexType PropertiesId)
{
    rModelPart.RemovePropertiesFromAllLevels(PropertiesId);
}

void ModelPartRemovePropertiesFromAllLevels2(ModelPart& rModelPart, ModelPart::IndexType PropertiesId, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemovePropertiesFromAllLevels(PropertiesId, ThisIndex);
}

void ModelPartRemovePropertiesFromAllLevels3(ModelPart& rModelPart, ModelPart::PropertiesType::Pointer pThisProperties)
{
    rModelPart.RemovePropertiesFromAllLevels(pThisProperties);
}

void ModelPartRemovePropertiesFromAllLevels4(ModelPart& rModelPart, ModelPart::PropertiesType::Pointer pThisProperties, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemovePropertiesFromAllLevels(pThisProperties, ThisIndex);
}

// Elements

ModelPart::SizeType ModelPartNumberOfElements1(ModelPart& rModelPart)
{
    return rModelPart.NumberOfElements();
}

ModelPart::ElementsContainerType::Pointer ModelPartGetElements1(ModelPart& rModelPart)
{
    return rModelPart.pElements();
}

ModelPart::ElementsContainerType::Pointer ModelPartGetElements2(ModelPart& rModelPart, ModelPart::IndexType ThisIndex)
{
    return rModelPart.pElements(ThisIndex);
}

void ModelPartSetElements1(ModelPart& rModelPart, ModelPart::ElementsContainerType::Pointer pOtherElements)
{
    rModelPart.SetElements(pOtherElements);
}

void ModelPartSetElements2(ModelPart& rModelPart, ModelPart::ElementsContainerType::Pointer pOtherElements, ModelPart::IndexType ThisIndex)
{
    rModelPart.SetElements(pOtherElements, ThisIndex);
}

ModelPart::ElementType::Pointer ModelPartGetElement1(ModelPart& rModelPart, ModelPart::IndexType ElementId)
{
    return rModelPart.pGetElement(ElementId);
}

ModelPart::ElementType::Pointer ModelPartGetElement2(ModelPart& rModelPart, ModelPart::IndexType ElementId, ModelPart::IndexType ThisIndex)
{
    return rModelPart.pGetElement(ElementId, ThisIndex);
}

void ModelPartRemoveElement1(ModelPart& rModelPart, ModelPart::IndexType ElementId)
{
    rModelPart.RemoveElement(ElementId);
}

void ModelPartRemoveElement2(ModelPart& rModelPart, ModelPart::IndexType ElementId, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemoveElement(ElementId, ThisIndex);
}

void ModelPartRemoveElement3(ModelPart& rModelPart, ModelPart::ElementType::Pointer pThisElement)
{
    rModelPart.RemoveElement(pThisElement);
}

void ModelPartRemoveElement4(ModelPart& rModelPart, ModelPart::ElementType::Pointer pThisElement, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemoveElement(pThisElement, ThisIndex);
}

void ModelPartRemoveElementFromAllLevels1(ModelPart& rModelPart, ModelPart::IndexType ElementId)
{
    rModelPart.RemoveElementFromAllLevels(ElementId);
}

void ModelPartRemoveElementFromAllLevels2(ModelPart& rModelPart, ModelPart::IndexType ElementId, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemoveElementFromAllLevels(ElementId, ThisIndex);
}

void ModelPartRemoveElementFromAllLevels3(ModelPart& rModelPart, ModelPart::ElementType::Pointer pThisElement)
{
    rModelPart.RemoveElementFromAllLevels(pThisElement);
}

void ModelPartRemoveElementFromAllLevels4(ModelPart& rModelPart, ModelPart::ElementType::Pointer pThisElement, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemoveElementFromAllLevels(pThisElement, ThisIndex);
}

void ModelPartRemoveElementsFromAllLevels(ModelPart& rModelPart, Flags identifier_flag)
{
    rModelPart.RemoveElementsFromAllLevels(identifier_flag);
}
// Conditions

ModelPart::SizeType ModelPartNumberOfConditions1(ModelPart& rModelPart)
{
    return rModelPart.NumberOfConditions();
}

ModelPart::ConditionsContainerType::Pointer ModelPartGetConditions1(ModelPart& rModelPart)
{
    return rModelPart.pConditions();
}

ModelPart::ConditionsContainerType::Pointer ModelPartGetConditions2(ModelPart& rModelPart, ModelPart::IndexType ThisIndex)
{
    return rModelPart.pConditions(ThisIndex);
}

void ModelPartSetConditions1(ModelPart& rModelPart, ModelPart::ConditionsContainerType::Pointer pOtherConditions)
{
    rModelPart.SetConditions(pOtherConditions);
}

void ModelPartSetConditions2(ModelPart& rModelPart, ModelPart::ConditionsContainerType::Pointer pOtherConditions, ModelPart::IndexType ThisIndex)
{
    rModelPart.SetConditions(pOtherConditions, ThisIndex);
}

void ModelPartAddCondition1(ModelPart& rModelPart, Condition::Pointer newCondition)
{
    rModelPart.AddCondition( newCondition );
}

void ModelPartAddCondition2(ModelPart& rModelPart, Condition::Pointer newCondition, unsigned int ThisIndex)
{
    rModelPart.AddCondition( newCondition, ThisIndex );
}

ModelPart::ConditionType::Pointer ModelPartGetCondition1(ModelPart& rModelPart, ModelPart::IndexType ConditionId)
{
    return rModelPart.pGetCondition(ConditionId);
}

ModelPart::ConditionType::Pointer ModelPartGetCondition2(ModelPart& rModelPart, ModelPart::IndexType ConditionId, ModelPart::IndexType ThisIndex)
{
    return rModelPart.pGetCondition(ConditionId, ThisIndex);
}

void ModelPartRemoveCondition1(ModelPart& rModelPart, ModelPart::IndexType ConditionId)
{
    rModelPart.RemoveCondition(ConditionId);
}

void ModelPartRemoveCondition2(ModelPart& rModelPart, ModelPart::IndexType ConditionId, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemoveCondition(ConditionId, ThisIndex);
}

void ModelPartRemoveCondition3(ModelPart& rModelPart, ModelPart::ConditionType::Pointer pThisCondition)
{
    rModelPart.RemoveCondition(pThisCondition);
}

void ModelPartRemoveCondition4(ModelPart& rModelPart, ModelPart::ConditionType::Pointer pThisCondition, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemoveCondition(pThisCondition, ThisIndex);
}

void ModelPartRemoveConditionFromAllLevels1(ModelPart& rModelPart, ModelPart::IndexType ConditionId)
{
    rModelPart.RemoveConditionFromAllLevels(ConditionId);
}

void ModelPartRemoveConditionFromAllLevels2(ModelPart& rModelPart, ModelPart::IndexType ConditionId, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemoveConditionFromAllLevels(ConditionId, ThisIndex);
}

void ModelPartRemoveConditionFromAllLevels3(ModelPart& rModelPart, ModelPart::ConditionType::Pointer pThisCondition)
{
    rModelPart.RemoveConditionFromAllLevels(pThisCondition);
}

void ModelPartRemoveConditionFromAllLevels4(ModelPart& rModelPart, ModelPart::ConditionType::Pointer pThisCondition, ModelPart::IndexType ThisIndex)
{
    rModelPart.RemoveConditionFromAllLevels(pThisCondition, ThisIndex);
}

void ModelPartRemoveConditionsFromAllLevels(ModelPart& rModelPart, Flags identifier_flag)
{
    rModelPart.RemoveConditionsFromAllLevels(identifier_flag);
}


// Geometries

void ModelPartAddGeometry1(ModelPart& rModelPart, ModelPart::GeometryType::Pointer pNewGeometry)
{
    rModelPart.AddGeometry(pNewGeometry);
}

ModelPart::GeometryType::Pointer ModelPartGetGeometry1(ModelPart& rModelPart, ModelPart::IndexType GeometryId)
{
    return rModelPart.pGetGeometry(GeometryId);
}

ModelPart::GeometryType::Pointer ModelPartGetGeometry2(ModelPart& rModelPart, const std::string& GeometryName)
{
    return rModelPart.pGetGeometry(GeometryName);
}

bool ModelPartHasGeometry1(ModelPart& rModelPart, ModelPart::IndexType GeometryId)
{
    return rModelPart.HasGeometry(GeometryId);
}

bool ModelPartHasGeometry2(ModelPart& rModelPart, const std::string& GeometryName)
{
    return rModelPart.HasGeometry(GeometryName);
}

void ModelPartRemoveGeometry1(ModelPart& rModelPart, ModelPart::IndexType GeometryId)
{
    rModelPart.RemoveGeometry(GeometryId);
}

void ModelPartRemoveGeometry2(ModelPart& rModelPart, const std::string& GeometryName)
{
    rModelPart.RemoveGeometry(GeometryName);
}

void ModelPartRemoveGeometryFromAllLevels1(ModelPart& rModelPart, ModelPart::IndexType GeometryId)
{
    rModelPart.RemoveGeometryFromAllLevels(GeometryId);
}

void ModelPartRemoveGeometryFromAllLevels2(ModelPart& rModelPart, const std::string& GeometryName)
{
    rModelPart.RemoveGeometryFromAllLevels(GeometryName);
}

// Master slave constraints
/* // Try with perfect forwarding
template <typename ... Args>
ModelPart::MasterSlaveConstraintType::Pointer CreateNewMasterSlaveConstraint(ModelPart& rModelPart, Args&& ... args)
{
    return rModelPart.CreateNewMasterSlaveConstraint(std::forward<Args>(args) ...);
}*/

ModelPart::MasterSlaveConstraintType::Pointer CreateNewMasterSlaveConstraint1(ModelPart& rModelPart,
                                                                              std::string ConstraintName,
                                                                              ModelPart::IndexType Id,
                                                                              ModelPart::DofsVectorType& rMasterDofsVector,
                                                                              ModelPart::DofsVectorType& rSlaveDofsVector,
                                                                              ModelPart::MatrixType RelationMatrix,
                                                                              ModelPart::VectorType ConstantVector)
{
    if (!KratosComponents<MasterSlaveConstraint>::Has(ConstraintName)) {
        std::stringstream msg;
        KratosComponents<MasterSlaveConstraint> instance; // creating an instance for using "PrintData"
        instance.PrintData(msg);

        KRATOS_ERROR << "The Constraint \"" << ConstraintName << "\" is not registered!\nMaybe you need to import the application where it is defined?\nThe following Constraints are registered:\n" << msg.str() << std::endl;
    }

    return rModelPart.CreateNewMasterSlaveConstraint(ConstraintName, Id, rMasterDofsVector, rSlaveDofsVector, RelationMatrix, ConstantVector);
}

// Master slave constraints

ModelPart::MasterSlaveConstraintType::Pointer CreateNewMasterSlaveConstraint2(ModelPart& rModelPart,
                                                                              std::string ConstraintName,
                                                                              ModelPart::IndexType Id,
                                                                              ModelPart::NodeType& rMasterNode,
                                                                              ModelPart::DoubleVariableType& rMasterVariable,
                                                                              ModelPart::NodeType& rSlaveNode,
                                                                              ModelPart::DoubleVariableType& rSlaveVariable,
                                                                              double Weight,
                                                                              double Constant)
{
    if (!KratosComponents<MasterSlaveConstraint>::Has(ConstraintName)) {
        std::stringstream msg;
        KratosComponents<MasterSlaveConstraint> instance; // creating an instance for using "PrintData"
        instance.PrintData(msg);

        KRATOS_ERROR << "The Constraint \"" << ConstraintName << "\" is not registered!\nMaybe you need to import the application where it is defined?\nThe following Constraints are registered:\n" << msg.str() << std::endl;
    }

    return rModelPart.CreateNewMasterSlaveConstraint(ConstraintName, Id, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant);
}

void ModelPartAddMasterSlaveConstraint(ModelPart& rModelPart, ModelPart::MasterSlaveConstraintType::Pointer pMasterSlaveConstraint)
{
    rModelPart.AddMasterSlaveConstraint(pMasterSlaveConstraint);
}

void AddMasterSlaveConstraintsByIds(ModelPart& rModelPart, std::vector< ModelPart::IndexType >& ConstraintIds )
{
    rModelPart.AddMasterSlaveConstraints(ConstraintIds);
}


const ModelPart::MasterSlaveConstraintContainerType& ModelPartGetMasterSlaveConstraints1(ModelPart& rModelPart)
{
    return rModelPart.MasterSlaveConstraints();
}

ModelPart::SizeType ModelPartNumberOfMasterSlaveConstraints1(ModelPart& rModelPart)
{
    return rModelPart.NumberOfMasterSlaveConstraints();
}

ModelPart::MasterSlaveConstraintType::Pointer ModelPartGetMasterSlaveConstraint1(ModelPart& rModelPart, ModelPart::IndexType MasterSlaveConstraintId)
{
    return rModelPart.pGetMasterSlaveConstraint(MasterSlaveConstraintId);
}

void ModelPartRemoveMasterSlaveConstraint1(ModelPart& rModelPart, ModelPart::IndexType MasterSlaveConstraintId)
{
    rModelPart.RemoveMasterSlaveConstraint(MasterSlaveConstraintId);
}

void ModelPartRemoveMasterSlaveConstraint2(ModelPart& rModelPart, ModelPart::MasterSlaveConstraintType& rOtherMasterSlaveConstraint)
{
    rModelPart.RemoveMasterSlaveConstraint(rOtherMasterSlaveConstraint);
}

void ModelPartRemoveMasterSlaveConstraintFromAllLevels1(ModelPart& rModelPart, ModelPart::IndexType MasterSlaveConstraintId)
{
    rModelPart.RemoveMasterSlaveConstraintFromAllLevels(MasterSlaveConstraintId);
}

void ModelPartRemoveMasterSlaveConstraintFromAllLevels2(ModelPart& rModelPart, ModelPart::MasterSlaveConstraintType& rMasterSlaveConstraint)
{
    rModelPart.RemoveMasterSlaveConstraintFromAllLevels(rMasterSlaveConstraint);
}

Communicator&  ModelPartGetCommunicator(ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator();
}

void RemoveSubModelPart1(ModelPart& rModelPart, std::string const& ThisSubModelPartName)
{
    rModelPart.RemoveSubModelPart(ThisSubModelPartName);
}

void RemoveSubModelPart2(ModelPart& rModelPart, ModelPart& ThisSubModelPart)
{
    rModelPart.RemoveSubModelPart(ThisSubModelPart);
}

void AddNodesByIds(ModelPart& rModelPart, std::vector< ModelPart::IndexType >& NodesIds )
{
    rModelPart.AddNodes(NodesIds);
}

void AddConditionsByIds(ModelPart& rModelPart,std::vector< ModelPart::IndexType >& ConditionsIds )
{
    rModelPart.AddConditions(ConditionsIds);
}

void AddElementsByIds(ModelPart& rModelPart, std::vector< ModelPart::IndexType >& ElementsIds )
{
    rModelPart.AddElements(ElementsIds);
}

const ModelPart::SubModelPartIterator GetSubModelPartBegin(ModelPart& rModelPart)
{
    return rModelPart.SubModelPartsBegin();
}

const ModelPart::SubModelPartIterator GetSubModelPartEnd(ModelPart& rModelPart)
{
    return rModelPart.SubModelPartsEnd();
}

/** Retrieve the variable names of the entities in the given container.
 * 
 * Retrieve the variable names of the entities in `rContainer`. If the
 * `doFullSearch` is enabled, it will iterate and check all the entities 
 * in the container. If not enabled it will be assumed that first entity of
 * the container is representative of the list of variables in every intenty
 */
template<class TContainerType>
const std::unordered_set<std::string> GetNonHistoricalVariablesNames(ModelPart& rModelPart, TContainerType& rContainer, bool doFullSearch=false) {

    std::unordered_set<std::string> variable_names;

    if(doFullSearch) {
        if(rContainer.size() == 0) {
            KRATOS_WARNING("DEBUG") << "Checking and empty container" << std::endl;
        } else {
            for(auto & variable: rContainer.begin()->Data()) {
                variable_names.insert(variable.first->Name());
            }
        }
    } else {
        if(rContainer.size() == 0) {
            KRATOS_WARNING("DEBUG") << "Checking and empty container" << std::endl;
        }
        for(auto & entity : rContainer) {
            for(auto & variable: entity.Data()) {
                variable_names.insert(variable.first->Name());
            }
        }
    }

    return variable_names;
}

void AddModelPartToPython(pybind11::module& m)
{

    ModelPart::IndexType(ModelPart::*pointer_to_clone_time_step_1)(void) = &ModelPart::CloneTimeStep;
    ModelPart::IndexType(ModelPart::*pointer_to_clone_time_step_2)(double) = &ModelPart::CloneTimeStep;
    ProcessInfo::Pointer(ModelPart::*pointer_to_get_process_info)(void) = &ModelPart::pGetProcessInfo;
    void (ModelPart::*pointer_to_set_process_info)(ProcessInfo::Pointer) = &ModelPart::SetProcessInfo;
    // ModelPart::MeshType::Pointer (ModelPart::*pointer_to_get_mesh)() = &ModelPart::pGetMesh;
    //      std::string& (ModelPart::*pointer_to_name)(void) = &ModelPart::Name;


    namespace py = pybind11;

    py::class_<typename ModelPart::SubModelPartsContainerType >(m, "SubModelPartsContainerType")
        .def("__iter__", [](typename ModelPart::SubModelPartsContainerType& self){ return py::make_iterator(self.begin(), self.end());},  py::keep_alive<0,1>())
        ;

    MapInterface<ModelPart::GeometriesMapType>().CreateInterface(m,"GeometriesMapType");
    PointerVectorSetPythonInterface<ModelPart::MasterSlaveConstraintContainerType>().CreateInterface(m,"MasterSlaveConstraintsArray");

    py::class_<ModelPart, Kratos::shared_ptr<ModelPart>, DataValueContainer, Flags >(m,"ModelPart")
        .def_property("Name", GetModelPartName, SetModelPartName)
        .def("FullName", &ModelPart::FullName)
        //  .def_property("ProcessInfo", GetProcessInfo, SetProcessInfo)
        .def_property("ProcessInfo", pointer_to_get_process_info, pointer_to_set_process_info)
        .def("CreateSolutionStep", &ModelPart::CreateSolutionStep)
        .def("CloneSolutionStep", &ModelPart::CloneSolutionStep)
        .def("CreateTimeStep", &ModelPart::CreateTimeStep)
        .def("ReduceTimeStep", &ModelPart::ReduceTimeStep)
        .def("CloneTimeStep", pointer_to_clone_time_step_1)
        .def("CloneTimeStep", pointer_to_clone_time_step_2)
        //       .def("CopySolutionStepData",&ModelPart::CopySolutionStepData)
        .def("NumberOfNodes", &ModelPart::NumberOfNodes)
        .def("NumberOfNodes", ModelPartNumberOfNodes1)
        .def("SetBufferSize", &ModelPart::SetBufferSize)
        .def("GetBufferSize", &ModelPart::GetBufferSize)
        .def("NumberOfElements", ModelPartNumberOfElements1)
        .def("NumberOfElements", &ModelPart::NumberOfElements)
        .def("NumberOfConditions", ModelPartNumberOfConditions1)
        .def("NumberOfConditions", &ModelPart::NumberOfConditions)
        .def("NumberOfGeometries", &ModelPart::NumberOfGeometries)
        .def("NumberOfMasterSlaveConstraints", ModelPartNumberOfMasterSlaveConstraints1)
        .def("NumberOfMasterSlaveConstraints", &ModelPart::NumberOfMasterSlaveConstraints)
        .def("NumberOfMeshes", &ModelPart::NumberOfMeshes)
        .def("NumberOfProperties", &ModelPart::NumberOfProperties)
        .def("NumberOfProperties", ModelPartNumberOfProperties1)
        .def("GetMesh", ModelPartGetMesh)
        .def("GetMesh", ModelPartGetMesh2)
        .def_property("Nodes", ModelPartGetNodes1, ModelPartSetNodes1)
        .def("GetNode", ModelPartGetNode1)
        .def("GetNode", ModelPartGetNode2)
        .def("GetNodes", ModelPartGetNodes1)
        .def("SetNodes", ModelPartSetNodes1)
        .def("GetNodes", ModelPartGetNodes2)
        .def("SetNodes", ModelPartSetNodes2)
        .def("RemoveNode", ModelPartRemoveNode1)
        .def("RemoveNode", ModelPartRemoveNode2)
        .def("RemoveNode", ModelPartRemoveNode3)
        .def("RemoveNode", ModelPartRemoveNode4)
        .def("RemoveNodes", &ModelPart::RemoveNodes)
        .def("RemoveNodeFromAllLevels", ModelPartRemoveNodeFromAllLevels1)
        .def("RemoveNodeFromAllLevels", ModelPartRemoveNodeFromAllLevels2)
        .def("RemoveNodeFromAllLevels", ModelPartRemoveNodeFromAllLevels3)
        .def("RemoveNodeFromAllLevels", ModelPartRemoveNodeFromAllLevels4)
        .def("RemoveNodesFromAllLevels", ModelPartRemoveNodesFromAllLevels)
        .def("NodesArray", &ModelPart::NodesArray, py::return_value_policy::reference_internal)
        .def("NumberOfTables", &ModelPart::NumberOfTables)
        .def("AddTable", &ModelPart::AddTable)
        .def("GetTable", &ModelPart::pGetTable)
        .def("HasProperties", ModelPartHasProperties1)
        .def("HasProperties", ModelPartHasProperties2)
        .def("HasProperties", ModelPartHasSubProperties1)
        .def("HasProperties", ModelPartHasSubProperties2)
        .def("RecursivelyHasProperties", ModelPartRecursivelyHasProperties1)
        .def("RecursivelyHasProperties", ModelPartRecursivelyHasProperties2)
        .def("CreateNewProperties", ModelPartCreateNewProperties1)
        .def("CreateNewProperties", ModelPartCreateNewProperties2)
        .def("GetProperties", ModelPartGetPropertiesDirect1)
        .def("GetProperties", ModelPartGetSubProperties1)
        .def_property("Properties", ModelPartGetProperties1, ModelPartSetProperties1)
        .def("AddProperties", ModelPartAddProperties1)
        .def("AddProperties", ModelPartAddProperties2)
        .def("GetProperties", ModelPartGetProperties1)
        .def("SetProperties", ModelPartSetProperties1)
        .def("GetProperties", ModelPartGetProperties2)
        .def("SetProperties", ModelPartSetProperties2)
        .def("RemoveProperties", ModelPartRemoveProperties1)
        .def("RemoveProperties", ModelPartRemoveProperties2)
        .def("RemoveProperties", ModelPartRemoveProperties3)
        .def("RemoveProperties", ModelPartRemoveProperties4)
        .def("RemovePropertiesFromAllLevels", ModelPartRemovePropertiesFromAllLevels1)
        .def("RemovePropertiesFromAllLevels", ModelPartRemovePropertiesFromAllLevels2)
        .def("RemovePropertiesFromAllLevels", ModelPartRemovePropertiesFromAllLevels3)
        .def("RemovePropertiesFromAllLevels", ModelPartRemovePropertiesFromAllLevels4)
        .def("PropertiesArray", &ModelPart::PropertiesArray, py::return_value_policy::reference_internal)
        .def_property("Elements", ModelPartGetElements1, ModelPartSetElements1)
        .def("GetElement", ModelPartGetElement1)
        .def("GetElement", ModelPartGetElement2)
        .def("GetElements", ModelPartGetElements1)
        .def("SetElements", ModelPartSetElements1)
        .def("GetElements", ModelPartGetElements2)
        .def("SetElements", ModelPartSetElements2)
        .def("RemoveElement", ModelPartRemoveElement1)
        .def("RemoveElement", ModelPartRemoveElement2)
        .def("RemoveElement", ModelPartRemoveElement3)
        .def("RemoveElement", ModelPartRemoveElement4)
        .def("RemoveElements", &ModelPart::RemoveElements)
        .def("RemoveElementFromAllLevels", ModelPartRemoveElementFromAllLevels1)
        .def("RemoveElementFromAllLevels", ModelPartRemoveElementFromAllLevels2)
        .def("RemoveElementFromAllLevels", ModelPartRemoveElementFromAllLevels3)
        .def("RemoveElementFromAllLevels", ModelPartRemoveElementFromAllLevels4)
        .def("RemoveElementsFromAllLevels", ModelPartRemoveElementsFromAllLevels)
        .def("ElementsArray", &ModelPart::ElementsArray, py::return_value_policy::reference_internal)
        .def_property("Conditions", ModelPartGetConditions1, ModelPartSetConditions1)
        .def("GetCondition", ModelPartGetCondition1)
        .def("GetCondition", ModelPartGetCondition2)
        .def("GetConditions", ModelPartGetConditions1)
        .def("SetConditions", ModelPartSetConditions1)
        .def("GetConditions", ModelPartGetConditions2)
        .def("SetConditions", ModelPartSetConditions2)
        .def("AddCondition", ModelPartAddCondition1)
        .def("AddCondition", ModelPartAddCondition2)
        .def("RemoveCondition", ModelPartRemoveCondition1)
        .def("RemoveCondition", ModelPartRemoveCondition2)
        .def("RemoveCondition", ModelPartRemoveCondition3)
        .def("RemoveCondition", ModelPartRemoveCondition4)
        .def("RemoveConditions", &ModelPart::RemoveConditions)
        .def("RemoveConditionFromAllLevels", ModelPartRemoveConditionFromAllLevels1)
        .def("RemoveConditionFromAllLevels", ModelPartRemoveConditionFromAllLevels2)
        .def("RemoveConditionFromAllLevels", ModelPartRemoveConditionFromAllLevels3)
        .def("RemoveConditionFromAllLevels", ModelPartRemoveConditionFromAllLevels4)
        .def("RemoveConditionsFromAllLevels", ModelPartRemoveConditionsFromAllLevels)
        .def("AddGeometry", ModelPartAddGeometry1)
        .def("GetGeometry", ModelPartGetGeometry1)
        .def("GetGeometry", ModelPartGetGeometry2)
        .def("HasGeometry", ModelPartHasGeometry1)
        .def("HasGeometry", ModelPartHasGeometry2)
        .def("RemoveGeometry", ModelPartRemoveGeometry1)
        .def("RemoveGeometry", ModelPartRemoveGeometry2)
        .def("RemoveGeometryFromAllLevels", ModelPartRemoveGeometryFromAllLevels1)
        .def("RemoveGeometryFromAllLevels", ModelPartRemoveGeometryFromAllLevels2)
        .def_property("Geometries", [](ModelPart& self) { return self.Geometries(); },
            [](ModelPart& self, ModelPart::GeometriesMapType& geometries) {
                KRATOS_ERROR << "Setting geometries is not allowed! Trying to set value of ModelPart::Geometries."; })
        .def("CreateSubModelPart", &ModelPart::CreateSubModelPart, py::return_value_policy::reference_internal)
        .def("NumberOfSubModelParts", &ModelPart::NumberOfSubModelParts)
        .def("GetSubModelPart", &ModelPart::GetSubModelPart, py::return_value_policy::reference_internal)
        .def("RemoveSubModelPart", RemoveSubModelPart1)
        .def("RemoveSubModelPart", RemoveSubModelPart2)
        .def("HasSubModelPart", &ModelPart::HasSubModelPart)
        .def("ConditionsArray", &ModelPart::ConditionsArray, py::return_value_policy::reference_internal)
        .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<bool>)
        .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<int>)
        .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<double>)
        .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<array_1d<double, 3 > >)
        .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<Vector>)
        .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<Matrix>)
        .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<Quaternion<double> >)
        .def("HasNodalSolutionStepVariable", HasNodalSolutionStepVariable<bool>)
        .def("HasNodalSolutionStepVariable", HasNodalSolutionStepVariable<int>)
        .def("HasNodalSolutionStepVariable", HasNodalSolutionStepVariable<double>)
        .def("HasNodalSolutionStepVariable", HasNodalSolutionStepVariable<array_1d<double, 3 > >)
        .def("HasNodalSolutionStepVariable", HasNodalSolutionStepVariable<Vector>)
        .def("HasNodalSolutionStepVariable", HasNodalSolutionStepVariable<Matrix>)
        .def("HasNodalSolutionStepVariable", HasNodalSolutionStepVariable<Quaternion<double> >)
        .def("GetNodalSolutionStepDataSize", &ModelPart::GetNodalSolutionStepDataSize)
        .def("GetNodalSolutionStepTotalDataSize", &ModelPart::GetNodalSolutionStepTotalDataSize)
        .def("OverwriteSolutionStepData", &ModelPart::OverwriteSolutionStepData)
        .def("CreateNewNode", ModelPartCreateNewNode)
        .def("CreateNewElement", ModelPartCreateNewElement)
        .def("CreateNewElement", [](ModelPart& rModelPart, const std::string& ElementName, ModelPart::IndexType Id,
            ModelPart::GeometryType::Pointer pGeometry, ModelPart::PropertiesType::Pointer pProperties)
            {return rModelPart.CreateNewElement(ElementName, Id, pGeometry, pProperties);})
        .def("CreateNewCondition", ModelPartCreateNewCondition)
        .def("CreateNewCondition", [](ModelPart& rModelPart, const std::string& ConditionName, ModelPart::IndexType Id,
            ModelPart::GeometryType::Pointer pGeometry, ModelPart::PropertiesType::Pointer pProperties)
            {return rModelPart.CreateNewCondition(ConditionName, Id, pGeometry, pProperties);})
        .def("GetCommunicator", ModelPartGetCommunicator, py::return_value_policy::reference_internal)
        .def("Check", &ModelPart::Check)
        .def("IsSubModelPart", &ModelPart::IsSubModelPart)
        .def("IsDistributed", &ModelPart::IsDistributed)
        .def("AddNode", &ModelPart::AddNode)
        .def("AddNodes",AddNodesByIds)
        .def("AddCondition", &ModelPart::AddCondition)
        .def("AddConditions",AddConditionsByIds)
        .def("AddElement", &ModelPart::AddElement)
        .def("AddElements",AddElementsByIds)
        .def("GetParentModelPart", [](ModelPart& self) -> ModelPart& {return self.GetParentModelPart();}, py::return_value_policy::reference_internal)
        .def("GetRootModelPart",   [](ModelPart& self) -> ModelPart& {return self.GetRootModelPart();},   py::return_value_policy::reference_internal)
        .def("GetModel", &ModelPart::GetModel, py::return_value_policy::reference_internal)
        .def_property("SubModelParts",  [](ModelPart& self){ return self.SubModelParts(); },
                                        [](ModelPart& self, ModelPart::SubModelPartsContainerType& subs){ KRATOS_ERROR << "setting submodelparts is not allowed"; })
        .def_property_readonly("MasterSlaveConstraints", ModelPartGetMasterSlaveConstraints1)
        .def("GetHistoricalVariablesNames", [](ModelPart& rModelPart) -> std::unordered_set<std::string> {
            std::unordered_set<std::string> variable_names;
            for(auto & variable: rModelPart.GetNodalSolutionStepVariablesList()) {
                variable_names.insert(variable.Name());
            }
            return variable_names;
        })
        .def("GetNonHistoricalVariablesNames", [](ModelPart& rModelPart, ModelPart::NodesContainerType& rContainer, bool doFullSearch=false) -> std::unordered_set<std::string> {
            return GetNonHistoricalVariablesNames(rModelPart, rContainer, doFullSearch);
        })
        .def("GetNonHistoricalVariablesNames", [](ModelPart& rModelPart, ModelPart::ElementsContainerType& rContainer, bool doFullSearch=false) -> std::unordered_set<std::string> {
            return GetNonHistoricalVariablesNames(rModelPart, rContainer, doFullSearch);
        })
        .def("GetNonHistoricalVariablesNames", [](ModelPart& rModelPart, ModelPart::ConditionsContainerType& rContainer, bool doFullSearch=false) -> std::unordered_set<std::string> {
            return GetNonHistoricalVariablesNames(rModelPart, rContainer, doFullSearch);
        })
        .def("GetMasterSlaveConstraint", ModelPartGetMasterSlaveConstraint1)
        .def("GetMasterSlaveConstraints", ModelPartGetMasterSlaveConstraints1)
        .def("RemoveMasterSlaveConstraint", ModelPartRemoveMasterSlaveConstraint1)
        .def("RemoveMasterSlaveConstraint", ModelPartRemoveMasterSlaveConstraint2)
        .def("RemoveMasterSlaveConstraintFromAllLevels", ModelPartRemoveMasterSlaveConstraintFromAllLevels1)
        .def("RemoveMasterSlaveConstraintFromAllLevels", ModelPartRemoveMasterSlaveConstraintFromAllLevels2)
        .def("RemoveMasterSlaveConstraints", &ModelPart::RemoveMasterSlaveConstraints)
        .def("RemoveMasterSlaveConstraintsFromAllLevels", &ModelPart::RemoveMasterSlaveConstraintsFromAllLevels)
        .def("AddMasterSlaveConstraint", ModelPartAddMasterSlaveConstraint)
        .def("AddMasterSlaveConstraints", AddMasterSlaveConstraintsByIds)
        .def("CreateNewMasterSlaveConstraint",CreateNewMasterSlaveConstraint1, py::return_value_policy::reference_internal)
        .def("CreateNewMasterSlaveConstraint",CreateNewMasterSlaveConstraint2, py::return_value_policy::reference_internal)
        .def("__str__", PrintObject<ModelPart>)
        ;
}

} // namespace Python.

} // Namespace Kratos
