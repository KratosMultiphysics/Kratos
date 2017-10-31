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
#include "custom_processes/recover_volume_losses_process.hpp"
#include "custom_processes/select_mesh_elements_for_fluids_process.hpp"
#include "custom_processes/generate_new_nodes_before_meshing_process.hpp"
#include "custom_processes/inlet_management_process.hpp"
#include "custom_processes/set_inlet_process.hpp"
#include "custom_processes/model_start_end_meshing_for_fluids_process.hpp"
#include "custom_processes/split_elements_process.hpp"
#include "custom_processes/set_active_flag_process.hpp"
#include "custom_processes/adaptive_time_interval_process.hpp"
#include "custom_processes/transfer_model_part_elements_process.hpp"

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



      class_<RecoverVolumeLossesProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "RecoverVolumeLosses", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>()
	 )
	;


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

      class_<InletManagementProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "InletManagement", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>()
      	 )
      	;

      class_<SetInletProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "SetInlet", init<ModelPart&, int>()
      	 )
      	;

      class_<SplitElementsProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "SplitElementsProcess", init<ModelPart&, int>()
	 )
	;

      class_<SetActiveFlagProcess,  bases<ProcessBaseType>, boost::noncopyable >
	(
	 "SetActiveFlagProcess", init<ModelPart&, bool, bool, int>()
	 )
	;

      class_<AdaptiveTimeIntervalProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "AdaptiveTimeIntervalProcess", init<ModelPart&, int>()
      	 )
      	;

     class_<ModelStartEndMeshingForFluidsProcess, bases<ModelStartEndMeshingProcessType>, boost::noncopyable >
	(
	 "ModelMeshingForFluids", init<ModelPart&, Flags, int>()
	 )
	;

      //**********TRANSFER ELEMENTS TO MODEL PART*********//

      class_<TransferModelPartElementsProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "TransferModelPartElementsProcess", init<ModelPart&, ModelPart&>()
      	)
        .def("Execute", &TransferModelPartElementsProcess::Execute)
      	;
      

    }	
 
 
  }  // namespace Python.

} // Namespace Kratos

