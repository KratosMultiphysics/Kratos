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


// #include "custom_processes/elemental_neighbours_search_process.hpp"
// #include "custom_processes/nodal_neighbours_search_process.hpp"
// #include "custom_processes/model_volume_calculation_process.hpp"
// #include "custom_processes/refine_mesh_elements_on_threshold_process.hpp" 


#include "custom_processes/build_mesh_boundary_process_for_fluids.hpp"
#include "custom_processes/model_start_end_meshing_process_for_fluids.hpp"
#include "custom_processes/refine_mesh_boundary_process_for_fluids.hpp"
#include "custom_processes/remove_mesh_nodes_process_for_fluids.hpp"
#include "custom_processes/reconstruct_mesh_boundary_process_for_fluids.hpp"

//Processes


namespace Kratos
{
	
  namespace Python
  {

    // typedef Process                        ProcessBaseType;
    // typedef Process::Pointer                ProcessPointer;
    // typedef std::vector<Process::Pointer> ProcessContainer;

    // void Push_Back_Process( ProcessContainer& ThisProcessContainer,
    // 			       ProcessPointer ThisProcess )
    // {
    //    ThisProcessContainer.push_back( ThisProcess );
    // }
  

    void  AddCustomProcessesToPython()
    {

      using namespace boost::python;
      typedef Process                                         ProcessBaseType;

      class_<BuildMeshBoundaryProcessForFluids, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "BuildMeshBoundaryForFluids", init<ModelPart&, int, int, int>()
      	 )
      	;

      class_<ModelStartEndMeshingProcessForFluids, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "ModelMeshingForFluids", init<ModelPart&, Flags, int>()
      	 )
      	;

      class_<RefineMeshBoundaryProcessForFluids, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "RefineMeshBoundaryForFluids", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int, int>()
      	 )
      	; 

      class_<RemoveMeshNodesProcessForFluids, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "RemoveMeshNodesForFluids", init<ModelPart&, ModelerUtilities::MeshingParameters&, int, int>()
      	 )
      	;

       class_<ReconstructMeshBoundaryProcessForFluids, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "ReconstructMeshBoundaryForFluids", init<ModelPart&, ModelerUtilities::MeshingParameters&, int, int>()
      	 )
      	; 




    }	
 
 
  }  // namespace Python.

} // Namespace Kratos

