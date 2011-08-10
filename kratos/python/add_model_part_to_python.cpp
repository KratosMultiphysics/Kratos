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

        Node < 3 > ::Pointer ModelPartCreateNewNode(ModelPart& rModelPart, int Id, double x, double y, double z)
        {
            return rModelPart.CreateNewNode(Id, x, y, z);
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
                    .def("AssembleCurrentData", CommunicatorAssembleCurrentData<double> )
                    .def("AssembleCurrentData", CommunicatorAssembleCurrentData<array_1d<double,3> > )
                    .def("AssembleCurrentData", CommunicatorAssembleCurrentData<Vector> )
                    .def("AssembleCurrentData", CommunicatorAssembleCurrentData<Matrix> )
                   ;



            class_<ModelPart > ("ModelPart")
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
                    .def("GetMesh", ModelPartGetMesh)
                    .add_property("Nodes", ModelPartGetNodes1, ModelPartSetNodes1)
                    .def("GetNodes", ModelPartGetNodes1)
                    .def("SetNodes", ModelPartSetNodes1)
                    .def("GetNodes", ModelPartGetNodes2)
                    .def("SetNodes", ModelPartSetNodes2)
                    .def("NodesArray", &ModelPart::NodesArray, return_internal_reference<>())
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
                    .def("ConditionsArray", &ModelPart::ConditionsArray, return_internal_reference<>())
                    .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<bool>)
                    .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<int>)
                    .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<double>)
                    .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<array_1d<double, 3 > >)
                    .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<Vector>)
                    .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<Matrix>)
                    .def("OverwriteSolutionStepData", &ModelPart::OverwriteSolutionStepData)
                    .def("CreateNewNode", ModelPartCreateNewNode)
                    .def("GetCommunicator", ModelPartGetCommunicator, return_internal_reference<>())
                    .def("Check", &ModelPart::Check)
                    //.def("",&ModelPart::)
                    .def(self_ns::str(self))
                    ;
        }

    } // namespace Python.

} // Namespace Kratos

