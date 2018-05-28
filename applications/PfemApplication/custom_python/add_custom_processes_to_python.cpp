//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

// System includes 

// External includes 

// Project includes
#include "custom_python/add_custom_processes_to_python.h"


// Processes
#include "custom_processes/adaptive_time_interval_process.hpp"
#include "custom_processes/inlet_management_process.hpp"
#include "custom_processes/set_inlet_process.hpp"

// MeshModeler initialization and finalization processes
#include "custom_processes/model_start_end_meshing_for_fluids_process.hpp"

// PreMeshing processes
#include "custom_processes/split_elements_process.hpp"

// MiddleMeshing processes
#include "custom_processes/remove_mesh_nodes_for_fluids_process.hpp"

// PostMeshing processes
#include "custom_processes/generate_new_nodes_before_meshing_process.hpp"
#include "custom_processes/select_mesh_elements_for_fluids_process.hpp"
#include "custom_processes/recover_volume_losses_process.hpp"
#include "custom_processes/set_active_flag_process.hpp"
#include "custom_processes/transfer_model_part_elements_process.hpp"


namespace Kratos
{
	
namespace Python
{
  	
void  AddCustomProcessesToPython(pybind11::module& m)
{

  using namespace pybind11;

  //**********MESH MODELLER PROCESS*********//
  
  class_<RemoveMeshNodesForFluidsProcess, RemoveMeshNodesForFluidsProcess::Pointer, Process>
      (m, "RemoveMeshNodesForFluids")
      .def(init<ModelPart&, ModelerUtilities::MeshingParameters&, int>());
  
  class_<SplitElementsProcess, SplitElementsProcess::Pointer, Process>
      (m,"SplitElementsProcess")
      .def(init<ModelPart&, int>());

  class_<GenerateNewNodesBeforeMeshingProcess, GenerateNewNodesBeforeMeshingProcess::Pointer, Process>
      (m, "GenerateNewNodesBeforeMeshing")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>());
  
  class_<SelectMeshElementsForFluidsProcess, SelectMeshElementsForFluidsProcess::Pointer, Process>
      (m, "SelectMeshElementsForFluids")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>());
  
  class_<SetActiveFlagProcess, SetActiveFlagProcess::Pointer, Process>
      (m, "SetActiveFlagProcess")
      .def(init<ModelPart&, bool, bool, int>());

  class_<ModelStartEndMeshingForFluidsProcess, ModelStartEndMeshingForFluidsProcess::Pointer, ModelStartEndMeshingProcess>
      (m, "ModelMeshingForFluids")
      .def(init<ModelPart&, Flags, int>());


  //*****TRANSFER ELEMENTS TO MODEL PART****//
  
  class_<TransferModelPartElementsProcess, TransferModelPartElementsProcess::Pointer, Process>
      (m, "TransferModelPartElementsProcess")
      .def(init<ModelPart&, ModelPart&>())
      ;

  //*********ADAPTIVE TIME STEP*************//
  
  class_<AdaptiveTimeIntervalProcess, AdaptiveTimeIntervalProcess::Pointer, Process>
      (m, "AdaptiveTimeIntervalProcess")
      .def(init<ModelPart&, int>());
  
  //***************INLET PROCESS************//
  
  class_<InletManagementProcess, InletManagementProcess::Pointer, Process>
      (m, "InletManagement")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>());

  class_<SetInletProcess, SetInletProcess::Pointer, Process>
      (m, "SetInlet")
      .def(init<ModelPart&, int>());
  
  
  //*********VOLUME RECOVETY PROCESS********//
  
  class_<RecoverVolumeLossesProcess, RecoverVolumeLossesProcess::Pointer, Process>
      (m, "RecoverVolumeLosses")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>());
  
}
 
}  // namespace Python.

} // Namespace Kratos

