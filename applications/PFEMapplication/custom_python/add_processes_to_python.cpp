/*
==============================================================================
KratosPFEMApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
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
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-22 17:13:57 $
//   Revision:            $Revision: 1.6 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_processes_to_python.h"
#include "custom_processes/act_on_walls_nodal_process.h" 
#include "custom_processes/find_nodal_h_respecting_distance_process.h"
#include "custom_processes/coordinate_laplacian_smoother_process.h" 
#include "custom_processes/move_mesh_process.h" 
// #include "custom_processes/node_erase_process.h" 
#include "custom_processes/ulf_time_step_dec_process.h" 
//#include "custom_processes/lagrangian_inlet_process.h" 
#include "includes/node.h"
// #include "custom_processes/create_internal_nodes_process.h" 
// #include "custom_processes/particle_to_node_projection_process.h" 
// #include "custom_processes/Lagrangian_move_process.h" 
// #include "custom_processes/ptn_pressure_projection_process.h" 
// #include "custom_processes/ptn_velocity_projection_process.h" 
// #include "custom_processes/find_poor_element_process.h" 



namespace Kratos
{
	
namespace Python
{
  void  AddProcessesToPython()
  {
	using namespace boost::python;

// 	  class_<FindNodalHProcess, bases<Process> >("FindNodalHProcess",
// 		 init<ModelPart&>())
// 		 ;

	  class_<ActOnWallsNodalProcess, bases<Process> >("ActOnWallsNodalProcess",
		 init<ModelPart&>())
		 ;

	  class_<MoveMeshProcess, bases<Process> >("MoveMeshProcess",
		 init<ModelPart&>())
		 ;

/*	  class_<LagrangianInletProcess, bases<Process> >("LagrangianInletProcess",
		 init<ModelPart&, double>())
		 ;
*/
	  class_<CoordinateLaplacianSmootherProcess, bases<Process> >("CoordinateLaplacianSmootherProcess",
		 init<ModelPart&, int, double>())
		 ;

// 	  class_<NodeEraseProcess, bases<Process> >("NodeEraseProcess",
// 		 init<ModelPart&>())
// 		 ;
         class_<UlfTimeStepDecProcess>("UlfTimeStepDecProcess", init<ModelPart&>())
		.def("EstimateDeltaTime",&UlfTimeStepDecProcess::EstimateDeltaTime)
		;


// 	  class_<CreateInternalNodesProcess, bases<Process> >("CreateInternalNodesProcess",
// 		 init<ModelPart&, ModelPart::ElementsContainerType&, int>())
// 		 ;
// 
// 	  class_<ParticleToNodeProjectionProcess, bases<Process> >("ParticleToNodeProjectionProcess",
// 		 init<ModelPart&, double>())
// 		 ;
// 
// // 	  class_<MoveParticlesProcess, bases<Process> >("MoveParticlesProcess",
// // 		 init<ModelPart&>())
// // 		 ;
// 
// 	  class_<LagrangianMoveProcess, bases<Process> >("LagrangianMoveProcess",
// 		 init<ModelPart&, int, double>())
// 		 ;
// 
// 	  class_<PtNPressureProjectionProcess, bases<Process> >("PtNPressureProjectionProcess",
// 		 init<ModelPart&, double>())
// 		 ;
// 
// 
// 	  class_<PtNVelocityProjectionProcess, bases<Process> >("PtNVelocityProjectionProcess",
// 		 init<ModelPart&, double>())
// 		 ;
// 
// 	  class_<FindPoorElementProcess, bases<Process> >("FindPoorElementProcess",
// 		 init<ModelPart::ElementsContainerType&, ModelPart&, int,  double, int>())
// 		 ;

  }
	
}  // namespace Python.

} // Namespace Kratos

