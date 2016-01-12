/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//
//   Project Name:        Kratos
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-06-20 17:38:25 $
//   Revision:            $Revision: 1.14 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "python/add_model_part_to_python.h"
#include "includes/process_info.h"

namespace Kratos
{

namespace Python
{

template<class TDataType>
void AddNodalSolutionStepVariable(ModelPart& rModelPart, Variable<TDataType> const& rThisVariable)
{
    rModelPart.AddNodalSolutionStepVariable(rThisVariable);
}

void SetModelPartName(ModelPart& rModelPart, std::string const& NewName)
{
    rModelPart.Name() = NewName;
}

std::string GetModelPartName(ModelPart const& rModelPart)
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
        rModelPart.GetMeshes().push_back(empty_mesh.Clone());

    return rModelPart.pGetMesh(MeshIndex);
}

Node < 3 > ::Pointer ModelPartCreateNewNode(ModelPart& rModelPart, int Id, double x, double y, double z)
{
    return rModelPart.CreateNewNode(Id, x, y, z);
}

Element::Pointer ModelPartCreateNewElement(ModelPart& rModelPart, std::string ElementName, ModelPart::IndexType Id, boost::python::list& NodeIdList, ModelPart::PropertiesType::Pointer pProperties)
{
    Geometry< Node < 3 > >::PointsArrayType pElementNodeList;
    
    for(unsigned int i = 0; i < len(NodeIdList); i++) {
        pElementNodeList.push_back(rModelPart.pGetNode(boost::python::extract<int>(NodeIdList[i])));
    }

    return rModelPart.CreateNewElement(ElementName, Id, pElementNodeList, pProperties);
}

Condition::Pointer ModelPartCreateNewCondition(ModelPart& rModelPart, std::string ConditionName, ModelPart::IndexType Id, boost::python::list& NodeIdList, ModelPart::PropertiesType::Pointer pProperties)
{  
    Geometry< Node < 3 > >::PointsArrayType pConditionNodeList;
    
    for(unsigned int i = 0; i < len(NodeIdList); i++) {
        pConditionNodeList.push_back(rModelPart.pGetNode(boost::python::extract<int>(NodeIdList[i])));
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



void AddModelPartToPython()
{

    ModelPart::IndexType(ModelPart::*pointer_to_clone_time_step_1)(void) = &ModelPart::CloneTimeStep;
    ModelPart::IndexType(ModelPart::*pointer_to_clone_time_step_2)(double) = &ModelPart::CloneTimeStep;
    ProcessInfo::Pointer(ModelPart::*pointer_to_get_process_info)(void) = &ModelPart::pGetProcessInfo;
    void (ModelPart::*pointer_to_set_process_info)(ProcessInfo::Pointer) = &ModelPart::SetProcessInfo;
    // ModelPart::MeshType::Pointer (ModelPart::*pointer_to_get_mesh)() = &ModelPart::pGetMesh;
    //	  std::string& (ModelPart::*pointer_to_name)(void) = &ModelPart::Name;

    using namespace boost::python;

    class_<Communicator > ("Communicator")
    .def(init<>())
    .def("MyPID", &Communicator::MyPID)
    .def("TotalProcesses", &Communicator::TotalProcesses)
    .def("GetNumberOfColors", &Communicator::GetNumberOfColors)
    .def("NeighbourIndices", NeighbourIndicesConst, return_internal_reference<>())
    .def("SynchronizeNodalSolutionStepsData", &Communicator::SynchronizeNodalSolutionStepsData)
    .def("SynchronizeDofs", &Communicator::SynchronizeDofs)
    .def("LocalMesh", CommunicatorGetLocalMesh, return_internal_reference<>() )
    .def("LocalMesh", CommunicatorGetLocalMeshWithIndex, return_internal_reference<>() )
    .def("GhostMesh", CommunicatorGetGhostMesh, return_internal_reference<>() )
    .def("GhostMesh", CommunicatorGetGhostMeshWithIndex, return_internal_reference<>() )
    .def("InterfaceMesh", CommunicatorGetInterfaceMesh, return_internal_reference<>() )
    .def("InterfaceMesh", CommunicatorGetInterfaceMeshWithIndex, return_internal_reference<>() )
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



	class_<ModelPart, bases<DataValueContainer, Flags> >("ModelPart")
		.def(init<std::string const&>())
		.def(init<>())
		.add_property("Name", GetModelPartName, SetModelPartName)
		//  .add_property("ProcessInfo", GetProcessInfo, SetProcessInfo)
		.add_property("ProcessInfo", pointer_to_get_process_info, pointer_to_set_process_info)
		.def("CreateSolutionStep", &ModelPart::CreateSolutionStep)
		.def("CloneSolutionStep", &ModelPart::CloneSolutionStep)
		.def("CreateTimeStep", &ModelPart::CreateTimeStep)
		.def("ReduceTimeStep", &ModelPart::ReduceTimeStep)
		.def("CloneTimeStep", pointer_to_clone_time_step_1)
		.def("CloneTimeStep", pointer_to_clone_time_step_2)
		// 	  .def("CopySolutionStepData",&ModelPart::CopySolutionStepData)
		.def("NumberOfNodes", &ModelPart::NumberOfNodes)
		.def("NumberOfNodes", ModelPartNumberOfNodes1)
		.def("SetBufferSize", &ModelPart::SetBufferSize)
		.def("GetBufferSize", &ModelPart::GetBufferSize)
		.def("NumberOfElements", ModelPartNumberOfElements1)
		.def("NumberOfElements", &ModelPart::NumberOfElements)
		.def("NumberOfConditions", ModelPartNumberOfConditions1)
		.def("NumberOfConditions", &ModelPart::NumberOfConditions)
		.def("NumberOfMeshes", &ModelPart::NumberOfMeshes)
		.def("NumberOfProperties", &ModelPart::NumberOfProperties)
		.def("NumberOfProperties", ModelPartNumberOfProperties1)
		.def("GetMesh", ModelPartGetMesh)
		.def("GetMesh", ModelPartGetMesh2)
		.add_property("Nodes", ModelPartGetNodes1, ModelPartSetNodes1)
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
		.def("NodesArray", &ModelPart::NodesArray, return_internal_reference<>())
		.def("NumberOfTables", &ModelPart::NumberOfTables)
		.def("AddTable", &ModelPart::AddTable)
		.def("GetTable", &ModelPart::pGetTable)
		.add_property("Properties", ModelPartGetProperties1, ModelPartSetProperties1)
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
		.def("PropertiesArray", &ModelPart::PropertiesArray, return_internal_reference<>())
		.add_property("Elements", ModelPartGetElements1, ModelPartSetElements1)
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
		.def("RemoveElementFromAllLevels", ModelPartRemoveElementFromAllLevels1)
		.def("RemoveElementFromAllLevels", ModelPartRemoveElementFromAllLevels2)
		.def("RemoveElementFromAllLevels", ModelPartRemoveElementFromAllLevels3)
		.def("RemoveElementFromAllLevels", ModelPartRemoveElementFromAllLevels4)
		.def("ElementsArray", &ModelPart::ElementsArray, return_internal_reference<>())
		.add_property("Conditions", ModelPartGetConditions1, ModelPartSetConditions1)
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
		.def("RemoveConditionFromAllLevels", ModelPartRemoveConditionFromAllLevels1)
		.def("RemoveConditionFromAllLevels", ModelPartRemoveConditionFromAllLevels2)
		.def("RemoveConditionFromAllLevels", ModelPartRemoveConditionFromAllLevels3)
		.def("RemoveConditionFromAllLevels", ModelPartRemoveConditionFromAllLevels4)
		.def("CreateSubModelPart", &ModelPart::CreateSubModelPart, return_internal_reference<>())
		.def("NumberOfSubModelParts", &ModelPart::NumberOfSubModelParts)
		.def("GetSubModelPart", &ModelPart::GetSubModelPart, return_internal_reference<>())
		.def("RemoveSubModelPart", RemoveSubModelPart1)
		.def("RemoveSubModelPart", RemoveSubModelPart2)
		.def("HasSubModelPart", &ModelPart::HasSubModelPart)
		.def("ConditionsArray", &ModelPart::ConditionsArray, return_internal_reference<>())
		.def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<bool>)
		.def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<int>)
		.def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<double>)
		.def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<array_1d<double, 3 > >)
		.def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<Vector>)
		.def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<Matrix>)
		.def("GetNodalSolutionStepDataSize", &ModelPart::GetNodalSolutionStepDataSize)
		.def("GetNodalSolutionStepTotalDataSize", &ModelPart::GetNodalSolutionStepTotalDataSize)
		.def("OverwriteSolutionStepData", &ModelPart::OverwriteSolutionStepData)
		.def("CreateNewNode", ModelPartCreateNewNode)
		.def("CreateNewElement", ModelPartCreateNewElement)
		.def("CreateNewCondition", ModelPartCreateNewCondition)
		.def("GetCommunicator", ModelPartGetCommunicator, return_internal_reference<>())
		.def("Check", &ModelPart::Check)
		.def("IsSubModelPart", &ModelPart::IsSubModelPart)
		//.def("",&ModelPart::)
		.def(self_ns::str(self))
		;
}

} // namespace Python.

} // Namespace Kratos

