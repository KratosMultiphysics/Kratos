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
#include "includes/model_part.h"
#include "python/add_model_part_to_python.h"
#include "includes/process_info.h"
#include "utilities/quaternion.h"
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
    Geometry< Node < 3 > >::PointsArrayType pElementNodeList;

    for(unsigned int i = 0; i < NodeIdList.size(); i++) {
        pElementNodeList.push_back(rModelPart.pGetNode(NodeIdList[i]));
    }

    return rModelPart.CreateNewElement(ElementName, Id, pElementNodeList, pProperties);
}

Condition::Pointer ModelPartCreateNewCondition(ModelPart& rModelPart, const std::string ConditionName, ModelPart::IndexType Id, std::vector< ModelPart::IndexType >& NodeIdList, ModelPart::PropertiesType::Pointer pProperties)
{
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

Properties::Pointer ModelPartGetPropertiesById(ModelPart& rModelPart, unsigned int property_id, unsigned int mesh_id)
{
    return rModelPart.pGetProperties(property_id, mesh_id);
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
    return rModelPart.CreateNewMasterSlaveConstraint(ConstraintName, Id, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant);
}

// Master slave constraints
ModelPart::MasterSlaveConstraintType::Pointer CreateNewMasterSlaveConstraint3(ModelPart& rModelPart,
                                                                              std::string ConstraintName,
                                                                              ModelPart::IndexType Id,
                                                                              ModelPart::NodeType& rMasterNode,
                                                                              ModelPart::VariableComponentType& rMasterVariable,
                                                                              ModelPart::NodeType& rSlaveNode,
                                                                              ModelPart::VariableComponentType& rSlaveVariable,
                                                                              double Weight,
                                                                              double Constant)
{
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



// Communicator

ModelPart::MeshType& CommunicatorGetLocalMesh(Communicator& rCommunicator)
{
    return rCommunicator.LocalMesh();
}


ModelPart::MeshType& CommunicatorGetLocalMeshWithIndex(Communicator& rCommunicator, Communicator::IndexType Index)
{
    return rCommunicator.LocalMesh(Index);
}

ModelPart::MeshType& CommunicatorGetGhostMesh(Communicator& rCommunicator)
{
    return rCommunicator.GhostMesh();
}

ModelPart::MeshType& CommunicatorGetGhostMeshWithIndex(Communicator& rCommunicator, Communicator::IndexType Index)
{
    return rCommunicator.GhostMesh(Index);
}

ModelPart::MeshType& CommunicatorGetInterfaceMesh(Communicator& rCommunicator)
{
    return rCommunicator.InterfaceMesh();
}

ModelPart::MeshType& CommunicatorGetInterfaceMeshWithIndex(Communicator& rCommunicator, Communicator::IndexType Index)
{
    return rCommunicator.InterfaceMesh(Index);
}

Communicator::NeighbourIndicesContainerType const&  NeighbourIndicesConst(Communicator& rCommunicator)
{
    return rCommunicator.NeighbourIndices();
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

template<class TDataType>
bool CommunicatorAssembleCurrentData(Communicator& rCommunicator, Variable<TDataType> const& ThisVariable)
{
    return rCommunicator.AssembleCurrentData(ThisVariable);
}

template<class TDataType>
bool CommunicatorAssembleNonHistoricalData(Communicator& rCommunicator, Variable<TDataType> const& ThisVariable)
{
    return rCommunicator.AssembleNonHistoricalData(ThisVariable);
}

template<class TDataType>
TDataType CommunicatorSumAll(Communicator& rCommunicator, const TDataType& rValue)
{
    TDataType Value = rValue;
    rCommunicator.SumAll(Value);
    return Value;
}

template<class TDataType>
TDataType CommunicatorMinAll(Communicator& rCommunicator, const TDataType& rValue)
{
    TDataType Value = rValue;
    rCommunicator.MinAll(Value);
    return Value;
}

template<class TDataType>
TDataType CommunicatorMaxAll(Communicator& rCommunicator, const TDataType& rValue)
{
    TDataType Value = rValue;
    rCommunicator.MaxAll(Value);
    return Value;
}

template<class TDataType>
TDataType CommunicatorScanSum(Communicator& rCommunicator, const TDataType rSendPartial, TDataType rReceiveAccumulated)
{
    TDataType SendPartial = rSendPartial;
    TDataType ReceiveAccumulated = rReceiveAccumulated;
    rCommunicator.ScanSum(SendPartial, ReceiveAccumulated);
    return ReceiveAccumulated;
}




void AddModelPartToPython(pybind11::module& m)
{

    ModelPart::IndexType(ModelPart::*pointer_to_clone_time_step_1)(void) = &ModelPart::CloneTimeStep;
    ModelPart::IndexType(ModelPart::*pointer_to_clone_time_step_2)(double) = &ModelPart::CloneTimeStep;
    ProcessInfo::Pointer(ModelPart::*pointer_to_get_process_info)(void) = &ModelPart::pGetProcessInfo;
    void (ModelPart::*pointer_to_set_process_info)(ProcessInfo::Pointer) = &ModelPart::SetProcessInfo;
    // ModelPart::MeshType::Pointer (ModelPart::*pointer_to_get_mesh)() = &ModelPart::pGetMesh;
    //      std::string& (ModelPart::*pointer_to_name)(void) = &ModelPart::Name;


    using namespace pybind11;

    class_<Communicator > (m,"Communicator")
        .def(init<>())
        .def("MyPID", &Communicator::MyPID)
        .def("Barrier", &Communicator::Barrier)
        .def("TotalProcesses", &Communicator::TotalProcesses)
        .def("GetNumberOfColors", &Communicator::GetNumberOfColors)
        .def("NeighbourIndices", NeighbourIndicesConst, return_value_policy::reference_internal)
        .def("SynchronizeNodalSolutionStepsData", &Communicator::SynchronizeNodalSolutionStepsData)
        .def("SynchronizeDofs", &Communicator::SynchronizeDofs)
        .def("LocalMesh", CommunicatorGetLocalMesh, return_value_policy::reference_internal )
        .def("LocalMesh", CommunicatorGetLocalMeshWithIndex, return_value_policy::reference_internal )
        .def("GhostMesh", CommunicatorGetGhostMesh, return_value_policy::reference_internal )
        .def("GhostMesh", CommunicatorGetGhostMeshWithIndex, return_value_policy::reference_internal )
        .def("InterfaceMesh", CommunicatorGetInterfaceMesh, return_value_policy::reference_internal )
        .def("InterfaceMesh", CommunicatorGetInterfaceMeshWithIndex, return_value_policy::reference_internal )
        .def("SumAll", CommunicatorSumAll<int> )
        .def("SumAll", CommunicatorSumAll<double> )
        .def("SumAll", CommunicatorSumAll<array_1d<double,3> > )
        .def("MinAll", CommunicatorMinAll<int> )
        .def("MinAll", CommunicatorMinAll<double> )
        .def("MaxAll", CommunicatorMaxAll<int> )
        .def("MaxAll", CommunicatorMaxAll<double> )
        .def("ScanSum", CommunicatorScanSum<int> )
        .def("ScanSum", CommunicatorScanSum<double> )
        .def("AssembleCurrentData", CommunicatorAssembleCurrentData<int> )
        .def("AssembleCurrentData", CommunicatorAssembleCurrentData<double> )
        .def("AssembleCurrentData", CommunicatorAssembleCurrentData<array_1d<double,3> > )
        .def("AssembleCurrentData", CommunicatorAssembleCurrentData<Vector> )
        .def("AssembleCurrentData", CommunicatorAssembleCurrentData<Matrix> )
        .def("AssembleNonHistoricalData", CommunicatorAssembleNonHistoricalData<int> )
        .def("AssembleNonHistoricalData", CommunicatorAssembleNonHistoricalData<double> )
        .def("AssembleNonHistoricalData", CommunicatorAssembleNonHistoricalData<array_1d<double,3> > )
        .def("AssembleNonHistoricalData", CommunicatorAssembleNonHistoricalData<Vector> )
        .def("AssembleNonHistoricalData", CommunicatorAssembleNonHistoricalData<Matrix> )
        ;

        class_<typename ModelPart::SubModelPartsContainerType >(m, "SubModelPartsContainerType")
        .def("__iter__", [](typename ModelPart::SubModelPartsContainerType& self){ return make_iterator(self.begin(), self.end());},  keep_alive<0,1>())
        ;

    PointerVectorSetPythonInterface<ModelPart::MasterSlaveConstraintContainerType>().CreateInterface(m,"MasterSlaveConstraintsArray");

    class_<ModelPart, Kratos::shared_ptr<ModelPart>, DataValueContainer, Flags >(m,"ModelPart")
        .def(init<std::string const&>())
        .def(init<>())
        .def_property("Name", GetModelPartName, SetModelPartName)
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
        .def("RemoveNodeFromAllLevels", ModelPartRemoveNodeFromAllLevels1)
        .def("RemoveNodeFromAllLevels", ModelPartRemoveNodeFromAllLevels2)
        .def("RemoveNodeFromAllLevels", ModelPartRemoveNodeFromAllLevels3)
        .def("RemoveNodeFromAllLevels", ModelPartRemoveNodeFromAllLevels4)
        .def("RemoveNodesFromAllLevels", ModelPartRemoveNodesFromAllLevels)
        .def("NodesArray", &ModelPart::NodesArray, return_value_policy::reference_internal)
        .def("NumberOfTables", &ModelPart::NumberOfTables)
        .def("AddTable", &ModelPart::AddTable)
        .def("GetTable", &ModelPart::pGetTable)
        .def("GetProperties", ModelPartGetPropertiesById) //new method where one asks for one specific property on one given mesh
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
        .def("PropertiesArray", &ModelPart::PropertiesArray, return_value_policy::reference_internal)
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
        .def("ElementsArray", &ModelPart::ElementsArray, return_value_policy::reference_internal)
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
        .def("CreateSubModelPart", &ModelPart::CreateSubModelPart, return_value_policy::reference_internal)
        .def("NumberOfSubModelParts", &ModelPart::NumberOfSubModelParts)
        .def("GetSubModelPart", &ModelPart::GetSubModelPart, return_value_policy::reference_internal)
        .def("RemoveSubModelPart", RemoveSubModelPart1)
        .def("RemoveSubModelPart", RemoveSubModelPart2)
        .def("HasSubModelPart", &ModelPart::HasSubModelPart)
        .def("ConditionsArray", &ModelPart::ConditionsArray, return_value_policy::reference_internal)
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
        .def("CreateNewCondition", ModelPartCreateNewCondition)
        .def("GetCommunicator", ModelPartGetCommunicator, return_value_policy::reference_internal)
        .def("Check", &ModelPart::Check)
        .def("IsSubModelPart", &ModelPart::IsSubModelPart)
        .def("AddNode", &ModelPart::AddNode)
        .def("AddNodes",AddNodesByIds)
        .def("AddCondition", &ModelPart::AddCondition)
        .def("AddConditions",AddConditionsByIds)
        .def("AddElement", &ModelPart::AddElement)
        .def("AddElements",AddElementsByIds)
        .def("GetParentModelPart", &ModelPart::GetParentModelPart, return_value_policy::reference_internal)
        .def("GetRootModelPart", &ModelPart::GetRootModelPart, return_value_policy::reference_internal)
        .def_property("SubModelParts",  [](ModelPart& self){ return self.SubModelParts(); },
                                        [](ModelPart& self, ModelPart::SubModelPartsContainerType& subs){ KRATOS_ERROR << "setting submodelparts is not allowed"; })

        .def_property_readonly("MasterSlaveConstraints", ModelPartGetMasterSlaveConstraints1)
        .def("GetMasterSlaveConstraint", ModelPartGetMasterSlaveConstraint1)
        .def("GetMasterSlaveConstraints", ModelPartGetMasterSlaveConstraints1)
        .def("RemoveMasterSlaveConstraint", ModelPartRemoveMasterSlaveConstraint1)
        .def("RemoveMasterSlaveConstraint", ModelPartRemoveMasterSlaveConstraint2)
        .def("RemoveMasterSlaveConstraintFromAllLevels", ModelPartRemoveMasterSlaveConstraintFromAllLevels1)
        .def("RemoveMasterSlaveConstraintFromAllLevels", ModelPartRemoveMasterSlaveConstraintFromAllLevels2)
        .def("AddMasterSlaveConstraint", ModelPartAddMasterSlaveConstraint)
        .def("AddMasterSlaveConstraints", AddMasterSlaveConstraintsByIds)
        .def("CreateNewMasterSlaveConstraint",CreateNewMasterSlaveConstraint1, return_value_policy::reference_internal)
        .def("CreateNewMasterSlaveConstraint",CreateNewMasterSlaveConstraint2, return_value_policy::reference_internal)
        .def("CreateNewMasterSlaveConstraint",CreateNewMasterSlaveConstraint3, return_value_policy::reference_internal)
        KRATOS_DEF_PYTHON_STR(ModelPart)
        ;
        //KRATOS_DEF_PYTHON_STR(ModelPart, py_class);
}

} // namespace Python.

} // Namespace Kratos
