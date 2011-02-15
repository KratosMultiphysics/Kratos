/*
==============================================================================
KratosULFApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2009-01-22 17:13:54 $
//   Revision:            $Revision: 1.5 $
//
//
 

// System includes 

// External includes 
#include <boost/python.hpp>

 
// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_processes_to_python.h"

#include "custom_processes/find_nodal_h_process.h" 
//#include "custom_processes/node_erase_process.h" 
#include "custom_processes/pressure_calculate_process.h"
#include "custom_processes/mass_calculate_process.h"  
#include "custom_processes/ulf_apply_bc_process.h" 
#include "custom_processes/ulf_time_step_dec_process.h" 
#include "custom_processes/mark_fluid_process.h"  
#include "custom_processes/mark_close_nodes_process.h"  
#include "custom_processes/mark_outer_nodes_process.h"  
#include "custom_processes/save_structure_model_part_process.h" 
#include "custom_processes/save_structure_conditions_process.h" 
#include "custom_processes/merge_model_parts_process.h"
#include "custom_processes/save_fluid_only_process.h"  
#include "custom_processes/lagrangian_inlet_process.h"  

#include "includes/node.h"

namespace Kratos
{
	
namespace Python
{
  void  AddProcessesToPython()
  {
	using namespace boost::python;

	  class_<FindNodalHProcess, bases<Process> >("FindNodalHProcess",
		 init<ModelPart&>())
		 ;
/*
	  class_<ActOnWallsNodalProcess, bases<Process> >("ActOnWallsNodalProcess",
		 init<ModelPart&>())
		 ;
*/
/*	  class_<MoveMeshProcess, bases<Process> >("MoveMeshProcess",
		 init<ModelPart&>())
		 ;
*/
/*	  class_<LagrangianInletProcess, bases<Process> >("LagrangianInletProcess",
		 init<ModelPart&, double>())
		 ;
*/
/*
	  class_<CoordinateLaplacianSmootherProcess, bases<Process> >("CoordinateLaplacianSmootherProcess",
		 init<ModelPart&, int>())
		 ;
*/      
	/*  
	class_<NodeEraseProcess, bases<Process> >("NodeEraseProcess",
		 init<ModelPart&>())
		 ;
	 */
	class_<PressureCalculateProcess, bases<Process> >("PressureCalculateProcess",
		 init<ModelPart&, unsigned int>())
		 ;
	class_<MassCalculateProcess, bases<Process> >("MassCalculateProcess",
		 init<ModelPart&>())
		 ;

		 class_<UlfApplyBCProcess, bases<Process> >("UlfApplyBCProcess",
		 init<ModelPart&>())
		 ;
	   class_<UlfTimeStepDecProcess, bases<Process> >("UlfTimeStepDecProcess",init<ModelPart&>())
		   .def("EstimateDeltaTime",&UlfTimeStepDecProcess::EstimateDeltaTime)
		 ;
	   class_<MarkOuterNodesProcess, bases<Process> >("MarkOuterNodesProcess",init<ModelPart&>())
		   .def("MarkOuterNodes",&MarkOuterNodesProcess::MarkOuterNodes)
		 ;
	   class_<MarkFluidProcess, bases<Process> >("MarkFluidProcess",
		 init<ModelPart&>())
		 ;
	   class_<MarkCloseNodesProcess, bases<Process> >("MarkCloseNodesProcess", init<ModelPart&>())
		   .def("MarkCloseNodes", &MarkCloseNodesProcess::MarkCloseNodes)
		 ;
	   class_<SaveStructureModelPartProcess, bases<Process> >("SaveStructureModelPartProcess", init<>())
		   .def("SaveStructure", &SaveStructureModelPartProcess::SaveStructure)
		 ;
  	class_<SaveStructureConditionsProcess, bases<Process> >("SaveStructureConditionsProcess", init<>())
		   .def("SaveStructureConditions", &SaveStructureConditionsProcess::SaveStructureConditions)
		 ;
	   class_<MergeModelPartsProcess, bases<Process> >("MergeModelPartsProcess", init<> ())
		   .def("MergeParts", &MergeModelPartsProcess::MergeParts)
		 ;
	 class_<SaveFluidOnlyProcess, bases<Process> >("SaveFluidOnlyProcess", init<> ())
		   .def("SaveFluidOnly", &SaveFluidOnlyProcess::SaveFluidOnly)
		 ;
	class_<LagrangianInletProcess, bases<Process> >("LagrangianInletProcess",
		 init<ModelPart&, double,  array_1d<double,3> >())
		 ;
  }
	
}  // namespace Python.

} // Namespace Kratos


