//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
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

//PreMeshing processes
#include "includes/model_part.h"

//MiddleMeshing processes
// #include "custom_processes/refine_mesh_elements_on_size_process.hpp"
#include "custom_processes/remove_mesh_nodes_for_fluids_process.hpp"

//PostMeshing processes
#include "custom_processes/select_mesh_elements_for_fluids_process.hpp"
#include "custom_processes/generate_new_nodes_before_meshing_process.hpp"
#include "custom_processes/model_start_end_meshing_for_fluids_process.hpp"

//Processes


namespace Kratos
{
	
  namespace Python
  {


    void  AddCustomProcessesToPython()
    {

      using namespace boost::python;
      typedef Process                                         ProcessBaseType;
      typedef ModelStartEndMeshingProcess     ModelStartEndMeshingProcessType;

      class_<RemoveMeshNodesForFluidsProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "RemoveMeshNodesForFluids", init<ModelPart&, ModelerUtilities::MeshingParameters&, int>()
      	 )
      	;


      class_<GenerateNewNodesBeforeMeshingProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "GenerateNewNodesBeforeMeshing", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>()
      	 )
      	;

      class_<SelectMeshElementsForFluidsProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "SelectMeshElementsForFluids", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>()
	 )
	;

      class_<ModelStartEndMeshingForFluidsProcess, bases<ModelStartEndMeshingProcessType>, boost::noncopyable >
	(
	 "ModelMeshingForFluids", init<ModelPart&, Flags, int>()
	 )
	;

    }	
 
 
  }  // namespace Python.

} // Namespace Kratos

