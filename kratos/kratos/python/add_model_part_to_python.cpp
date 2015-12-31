// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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

boost::python::list GetMeshesFromModelPart( ModelPart& rModelPart )
{
    boost::python::list  meshes_list;
    
    ModelPart::MeshesContainerType& meshes = rModelPart.GetMeshes();
    unsigned int number_of_meshes = meshes.size();
    
    for(unsigned int i=0; i<number_of_meshes; i++)
        meshes_list.append( rModelPart.pGetMesh(i) );
    
    return( meshes_list );
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

Element::Pointer ModelPartCreateNewElement(ModelPart& rModelPart, std::string ElementName, ModelPart::IndexType Id, boost::python::list& NodeList, ModelPart::PropertiesType::Pointer pProperties)
{
    Geometry< Node < 3 > >::PointsArrayType pElementNodeList;
    
    for(unsigned int i = 0; i < len(NodeList); i++) {
        pElementNodeList.push_back(rModelPart.pGetNode(boost::python::extract<int>(NodeList[i])));
    }
    
    return rModelPart.CreateNewElement(ElementName, Id, pElementNodeList, pProperties);
}

Condition::Pointer ModelPartCreateNewCondition(ModelPart& rModelPart, std::string ConditionName, ModelPart::IndexType Id, boost::python::list& NodeList, ModelPart::PropertiesType::Pointer pProperties)
{  
    Geometry< Node < 3 > >::PointsArrayType pConditionNodeList;
    
    for(unsigned int i = 0; i < len(NodeList); i++) {
        pConditionNodeList.push_back(rModelPart.pGetNode(boost::python::extract<int>(NodeList[i])));
    }
    
    return rModelPart.CreateNewCondition(ConditionName, Id, pConditionNodeList, pProperties);
}


// Nodes

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

// Properties

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

// Elements

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

// Conditions

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



    class_<ModelPart, bases<DataValueContainer, Flags> > ("ModelPart")
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
    .def("SetBufferSize", &ModelPart::SetBufferSize)
    .def("GetBufferSize", &ModelPart::GetBufferSize)
    .def("NumberOfElements", &ModelPart::NumberOfElements)
    .def("NumberOfConditions", &ModelPart::NumberOfConditions)
    .def("NumberOfMeshes", &ModelPart::NumberOfMeshes)
    .def("NumberOfProperties", &ModelPart::NumberOfProperties)
    .def("GetMesh", ModelPartGetMesh)
    .def("GetMesh", ModelPartGetMesh2)
    .def("GetListOfMeshes", GetMeshesFromModelPart)    
    .add_property("Nodes", ModelPartGetNodes1, ModelPartSetNodes1)
    .def("GetNodes", ModelPartGetNodes1)
    .def("SetNodes", ModelPartSetNodes1)
    .def("GetNodes", ModelPartGetNodes2)
    .def("SetNodes", ModelPartSetNodes2)
    .def("NodesArray", &ModelPart::NodesArray, return_internal_reference<>())
    .def("AddTable", &ModelPart::AddTable)
    .def("GetTable", &ModelPart::pGetTable)
    .add_property("Properties", ModelPartGetProperties1, ModelPartSetProperties1)
    .def("GetProperties", ModelPartGetProperties1)
    .def("SetProperties", ModelPartSetProperties1)
    .def("GetProperties", ModelPartGetProperties2)
    .def("SetProperties", ModelPartSetProperties2)
    .def("PropertiesArray", &ModelPart::PropertiesArray, return_internal_reference<>())
    .add_property("Elements", ModelPartGetElements1, ModelPartSetElements1)
    .def("GetElements", ModelPartGetElements1)
    .def("SetElements", ModelPartSetElements1)
    .def("GetElements", ModelPartGetElements2)
    .def("SetElements", ModelPartSetElements2)
    .def("ElementsArray", &ModelPart::ElementsArray, return_internal_reference<>())
    .add_property("Conditions", ModelPartGetConditions1, ModelPartSetConditions1)
    .def("GetConditions", ModelPartGetConditions1)
    .def("SetConditions", ModelPartSetConditions1)
    .def("GetConditions", ModelPartGetConditions2)
    .def("SetConditions", ModelPartSetConditions2)
    .def("AddCondition", ModelPartAddCondition1)
    .def("AddCondition", ModelPartAddCondition2)
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
    //.def("",&ModelPart::)
    .def(self_ns::str(self))
    ;
}

} // namespace Python.

} // Namespace Kratos

