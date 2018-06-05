//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

// External includes

// Project includes
#include "includes/node.h"
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


    void  AddCustomProcessesToPython(pybind11::module& m)
    {

      using namespace pybind11;
      typedef Process                                         ProcessBaseType;
      typedef ModelStartEndMeshingProcess     ModelStartEndMeshingProcessType;



      class_<RecoverVolumeLossesProcess, ProcessBaseType>
	(m, "RecoverVolumeLosses")
	.def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>());

      class_<RemoveMeshNodesForFluidsProcess, ProcessBaseType>
      	(m, "RemoveMeshNodesForFluids")
	.def(init<ModelPart&, ModelerUtilities::MeshingParameters&, int>());

      class_<GenerateNewNodesBeforeMeshingProcess, ProcessBaseType>
      	(m, "GenerateNewNodesBeforeMeshing")
	.def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>());

      class_<SelectMeshElementsForFluidsProcess, ProcessBaseType>
	(m, "SelectMeshElementsForFluids")
	.def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>());

      class_<InletManagementProcess, ProcessBaseType>
      	(m, "InletManagement")
	.def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>());

      class_<SetInletProcess, ProcessBaseType>
      	(m, "SetInlet")
	.def(init<ModelPart&, int>());

      class_<SplitElementsProcess, ProcessBaseType>
	(m,"SplitElementsProcess")
	.def(init<ModelPart&, int>());

      class_<SetActiveFlagProcess, ProcessBaseType>
	(m, "SetActiveFlagProcess")
	.def(init<ModelPart&, bool, bool, int>());

      class_<AdaptiveTimeIntervalProcess, ProcessBaseType>
      	(m, "AdaptiveTimeIntervalProcess")
	.def(init<ModelPart&, int>());

     class_<ModelStartEndMeshingForFluidsProcess, ModelStartEndMeshingProcessType>
       (m, "ModelMeshingForFluids")
       .def(init<ModelPart&, Flags, int>());

      //**********TRANSFER ELEMENTS TO MODEL PART*********//

      class_<TransferModelPartElementsProcess, ProcessBaseType>
      	(m, "TransferModelPartElementsProcess")
	  .def(init<ModelPart&, ModelPart&>())
        .def("Execute", &TransferModelPartElementsProcess::Execute)
      	;


    }


  }  // namespace Python.

} // Namespace Kratos

