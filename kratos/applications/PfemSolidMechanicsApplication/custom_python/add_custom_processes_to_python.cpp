//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
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
#include "custom_processes/elemental_neighbours_search_process.hpp"
#include "custom_processes/nodal_neighbours_search_process.hpp"
#include "custom_processes/boundary_skin_build_process.hpp"
#include "custom_processes/rigid_wall_contact_search_process.hpp"
#include "custom_processes/model_volume_calculation_process.hpp"

//Set initial mechanical state
#include "custom_processes/set_mechanical_initial_state_process.hpp"

//Modeler Bounding Boxes
#include "custom_modelers/spatial_bounding_box.hpp"
#include "custom_modelers/rigid_nose_wall_bounding_box.hpp"
#include "custom_modelers/rigid_circle_wall_bounding_box.hpp"
#include "custom_modelers/rigid_plane_wall_bounding_box.hpp"
namespace Kratos
{
	
  namespace Python
  {

  	
    void  AddCustomProcessesToPython()
    {

      using namespace boost::python;
      typedef Process                                         ProcessBaseType;


      //***************NEIGHBOURS**************//
      
      class_<NodalNeighboursSearchProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "NodalNeighboursSearch", init<ModelPart&, int, int, int, int>()
	 )
	.def("CleanNeighbours", &NodalNeighboursSearchProcess::ClearNeighbours)
	;
      
      class_<ElementalNeighboursSearchProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "ElementalNeighboursSearch", init<ModelPart&, int, int, int, int>()
	 )
	.def("CleanNeighbours", &ElementalNeighboursSearchProcess::ClearNeighbours)
	;

 
      //***************BOUNDARY**************//

      class_<BoundarySkinBuildProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "BuildBoundarySkin", init<ModelPart&, int, int, int>()
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


      //********MODEL VOLUME CALCULATION*********//

      class_<ModelVolumeCalculationProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "ModelVolumeCalculation", init<ModelPart&, int, int>()
	 )
	 .def("ExecuteInitializeSolutionStep", &ModelVolumeCalculationProcess::ExecuteInitializeSolutionStep)
	 .def("ExecuteFinalizeSolutionStep", &ModelVolumeCalculationProcess::ExecuteFinalizeSolutionStep)
	;

    }
 
  }  // namespace Python.

} // Namespace Kratos

