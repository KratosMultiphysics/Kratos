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
    ProcessInfo& GetProcessInfo(ModelPart& rModelPart )
    {	  return rModelPart.GetProcessInfo();	 }
    void SetProcessInfo(ModelPart& rModelPart, ProcessInfo& NewProcessInfo)
    {	  rModelPart.SetProcessInfo(NewProcessInfo);	  }
    
    ModelPart::MeshType::Pointer ModelPartGetMesh(ModelPart& rModelPart)
    {
	return rModelPart.pGetMesh();
    }
    
    
    Node<3>::Pointer ModelPartCreateNewNode(ModelPart& rModelPart, int Id, double x, double y, double z)
    {
	return rModelPart.CreateNewNode(Id,x,y,z);
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

  void  AddModelPartToPython()
  {

	  ModelPart::IndexType (ModelPart::*pointer_to_clone_time_step_1)(void) = &ModelPart::CloneTimeStep;
	  ModelPart::IndexType (ModelPart::*pointer_to_clone_time_step_2)(double) = &ModelPart::CloneTimeStep;
	  ProcessInfo& (ModelPart::*pointer_to_get_process_info)(void) = &ModelPart::GetProcessInfo;
	    // ModelPart::MeshType::Pointer (ModelPart::*pointer_to_get_mesh)() = &ModelPart::pGetMesh;
//	  std::string& (ModelPart::*pointer_to_name)(void) = &ModelPart::Name;

  using namespace boost::python;

  class_<ModelPart>("ModelPart")
	  .def(init<std::string const&>())
	  .def(init<>())
	  .add_property("Name", GetModelPartName, SetModelPartName)
	//  .add_property("ProcessInfo", GetProcessInfo, SetProcessInfo)
		  .add_property("ProcessInfo", make_function(pointer_to_get_process_info,return_internal_reference<>() ), &ModelPart::SetProcessInfo )
	  .def("CreateSolutionStep",&ModelPart::CreateSolutionStep)
	  .def("CloneSolutionStep",&ModelPart::CloneSolutionStep)
	  .def("CreateTimeStep",&ModelPart::CreateTimeStep)
          .def("ReduceTimeStep",&ModelPart::ReduceTimeStep)
	  .def("CloneTimeStep",pointer_to_clone_time_step_1)
	  .def("CloneTimeStep",pointer_to_clone_time_step_2)
// 	  .def("CopySolutionStepData",&ModelPart::CopySolutionStepData)
	  .def("NumberOfNodes",&ModelPart::NumberOfNodes)
	  .def("SetBufferSize",&ModelPart::SetBufferSize)
	  .def("GetBufferSize",&ModelPart::GetBufferSize)
	  .def("NumberOfElements",&ModelPart::NumberOfElements)
	  .def("NumberOfConditions",&ModelPart::NumberOfConditions)
	  .def("GetMesh",ModelPartGetMesh)
	  .add_property("Nodes", ModelPartGetNodes1, ModelPartSetNodes1)
	  .def("GetNodes",ModelPartGetNodes1)
	  .def("SetNodes",ModelPartSetNodes1)
	  .def("GetNodes",ModelPartGetNodes2)
	  .def("SetNodes",ModelPartSetNodes2)
	  .def("NodesArray", &ModelPart::NodesArray, return_internal_reference<>())
	  .add_property("Properties", ModelPartGetProperties1, ModelPartSetProperties1)
	  .def("GetProperties",ModelPartGetProperties1)
	  .def("SetProperties",ModelPartSetProperties1)
	  .def("GetProperties",ModelPartGetProperties2)
	  .def("SetProperties",ModelPartSetProperties2)
	  .def("PropertiesArray", &ModelPart::PropertiesArray, return_internal_reference<>())
	  .add_property("Elements", ModelPartGetElements1, ModelPartSetElements1)
	  .def("GetElements",ModelPartGetElements1)
	  .def("SetElements",ModelPartSetElements1)
	  .def("GetElements",ModelPartGetElements2)
	  .def("SetElements",ModelPartSetElements2)
	  .def("ElementsArray", &ModelPart::ElementsArray, return_internal_reference<>())
	  .add_property("Conditions", ModelPartGetConditions1, ModelPartSetConditions1)
	  .def("GetConditions",ModelPartGetConditions1)
	  .def("SetConditions",ModelPartSetConditions1)
	  .def("GetConditions",ModelPartGetConditions2)
	  .def("SetConditions",ModelPartSetConditions2)
	  .def("ConditionsArray", &ModelPart::ConditionsArray, return_internal_reference<>())
	  .def("AddNodalSolutionStepVariable",AddNodalSolutionStepVariable<bool>)
	  .def("AddNodalSolutionStepVariable",AddNodalSolutionStepVariable<int>)
	  .def("AddNodalSolutionStepVariable",AddNodalSolutionStepVariable<double>)
	  .def("AddNodalSolutionStepVariable",AddNodalSolutionStepVariable<array_1d<double,3> >)
	  .def("AddNodalSolutionStepVariable",AddNodalSolutionStepVariable<Vector>)
	  .def("AddNodalSolutionStepVariable",AddNodalSolutionStepVariable<Matrix>)
	  .def("OverwriteSolutionStepData",&ModelPart::OverwriteSolutionStepData)
	  .def("CreateNewNode",ModelPartCreateNewNode)
	  //.def("",&ModelPart::)
	  .def(self_ns::str(self))
	  ;
  }
	
}  // namespace Python.

} // Namespace Kratos

