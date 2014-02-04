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

//Modeler Bounding Boxes
#include "custom_modelers/rigid_wall_bounding_box.hpp"

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
	 "NodalNeighboursSearch", init<ModelPart&, unsigned int, unsigned int, unsigned int>()
	 )
	.def("CleanNeighbours", &NodalNeighboursSearchProcess::ClearNeighbours)
	;
      
      class_<ElementalNeighboursSearchProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "ElementalNeighboursSearch", init<ModelPart&, unsigned  int, unsigned int, unsigned int>()
	 )
	.def("CleanNeighbours", &ElementalNeighboursSearchProcess::ClearNeighbours)
	;

 
      //***************BOUNDARY**************//

      class_<BoundarySkinBuildProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "BuildBoundarySkin", init<ModelPart&, unsigned int, unsigned int>()
	 )
	;

      
      //********WALL CONTACT SEARCH*********//

      class_<RigidWallContactSearchProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "RigidWallContactSearch", init<RigidWallBoundingBox::Pointer, ModelPart&>()
	 )
	 .def("ExecuteInitializeSolutionStep", &RigidWallContactSearchProcess::ExecuteInitializeSolutionStep)
	 .def("ExecuteFinalizeSolutionStep", &RigidWallContactSearchProcess::ExecuteFinalizeSolutionStep)
	;


      //********MODEL VOLUME CALCULATION*********//

      class_<ModelVolumeCalculationProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "ModelVolumeCalculation", init<ModelPart&, int>()
	 )
	 .def("ExecuteInitializeSolutionStep", &ModelVolumeCalculationProcess::ExecuteInitializeSolutionStep)
	 .def("ExecuteFinalizeSolutionStep", &ModelVolumeCalculationProcess::ExecuteFinalizeSolutionStep)
	;

    }
 
  }  // namespace Python.

} // Namespace Kratos

