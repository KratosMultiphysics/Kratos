//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes 
#include <boost/python.hpp>

// External includes 

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "processes/process.h"

//Application includes
#include "custom_python/add_custom_processes_to_python.h"

//Processes
#include "custom_processes/rigid_wall_contact_search_process.hpp"

//Set initial mechanical state
#include "custom_processes/set_mechanical_initial_state_process.hpp"

//Modeler Bounding Boxes
#include "custom_bounding/rigid_nose_wall_bounding_box.hpp"
#include "custom_bounding/rigid_circle_wall_bounding_box.hpp"
#include "custom_bounding/rigid_plane_wall_bounding_box.hpp"

//Processes
#include "custom_processes/contact_refine_mesh_boundary_process.hpp"

namespace Kratos
{
	
  namespace Python
  {
 	
    void  AddCustomProcessesToPython()
    {

      using namespace boost::python;
      typedef Process                                         ProcessBaseType;
      typedef RefineMeshBoundaryProcess             RefineMeshProcessBaseType;
      typedef std::vector<SpatialBoundingBox::Pointer>   BoundingBoxContainer;

      
      //**********MESH MODELLER PROCESS*********//

      class_<ContactRefineMeshBoundaryProcess, bases<RefineMeshProcessBaseType>, boost::noncopyable >
	(
	 "ContactRefineMeshBoundary", init<ModelPart&, BoundingBoxContainer&, ModelerUtilities::MeshingParameters&, int, int>()
	 )
	;

      //********WALL CONTACT SEARCH*********//

      class_<RigidWallContactSearchProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "RigidNoseWallContactSearch", init<RigidNoseWallBoundingBox::Pointer, ModelPart&, int>()
      	 )
      	 .def("ExecuteInitializeSolutionStep", &RigidWallContactSearchProcess::ExecuteInitializeSolutionStep)
      	 .def("ExecuteFinalizeSolutionStep", &RigidWallContactSearchProcess::ExecuteFinalizeSolutionStep)
      	;

      class_<RigidWallContactSearchProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "RigidCircleWallContactSearch", init<RigidCircleWallBoundingBox::Pointer, ModelPart&, int>()
      	 )
      	 .def("ExecuteInitializeSolutionStep", &RigidWallContactSearchProcess::ExecuteInitializeSolutionStep)
      	 .def("ExecuteFinalizeSolutionStep", &RigidWallContactSearchProcess::ExecuteFinalizeSolutionStep)
      	;

      class_<RigidWallContactSearchProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "RigidPlaneWallContactSearch", init<RigidPlaneWallBoundingBox::Pointer, ModelPart&, int>()	 
      	 )
      	 .def("ExecuteInitializeSolutionStep", &RigidWallContactSearchProcess::ExecuteInitializeSolutionStep)
      	 .def("ExecuteFinalizeSolutionStep", &RigidWallContactSearchProcess::ExecuteFinalizeSolutionStep)
      	;
      
      // **** SET INITIAL MECHANICAL STATE **** //
      class_<SetMechanicalInitialStateProcess, bases<ProcessBaseType>, boost::noncopyable >
         (
          "SetMechanicalInitialStateProcess", init<ModelPart&, bool, double, double>()
         )
         .def("ExecuteInitialize",           &SetMechanicalInitialStateProcess::ExecuteInitialize)
         .def("ExecuteFinalizeSolutionStep", &SetMechanicalInitialStateProcess::ExecuteFinalizeSolutionStep)
         ;


 

    }
 
  }  // namespace Python.

} // Namespace Kratos

