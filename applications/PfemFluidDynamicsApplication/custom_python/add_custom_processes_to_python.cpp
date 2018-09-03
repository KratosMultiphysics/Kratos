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
#include "custom_processes/set_active_flag_mesher_process.hpp"
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
      typedef SettleModelStructureProcess     ModelStartEndMeshingProcessType;



      class_<RecoverVolumeLossesProcess, RecoverVolumeLossesProcess::Pointer, MesherProcess>
	(m, "RecoverVolumeLosses")
	.def(init<ModelPart&,  MesherUtilities::MeshingParameters&, int>());

      class_<RemoveMeshNodesForFluidsProcess, RemoveMeshNodesForFluidsProcess::Pointer, MesherProcess>
      	(m, "RemoveMeshNodesForFluids")
	.def(init<ModelPart&, MesherUtilities::MeshingParameters&, int>());

      class_<GenerateNewNodesBeforeMeshingProcess, GenerateNewNodesBeforeMeshingProcess::Pointer, MesherProcess>
      	(m, "GenerateNewNodesBeforeMeshing")
	.def(init<ModelPart&,  MesherUtilities::MeshingParameters&, int>());

      class_<SelectMeshElementsForFluidsProcess, SelectMeshElementsForFluidsProcess::Pointer, MesherProcess>
	(m, "SelectMeshElementsForFluids")
	.def(init<ModelPart&,  MesherUtilities::MeshingParameters&, int>());

      class_<InletManagementProcess, InletManagementProcess::Pointer, MesherProcess>
      	(m, "InletManagement")
	.def(init<ModelPart&,  MesherUtilities::MeshingParameters&, int>());

      class_<SetInletProcess, SetInletProcess::Pointer, ProcessBaseType>
      	(m, "SetInlet")
	.def(init<ModelPart&, int>());

      class_<SplitElementsProcess, SplitElementsProcess::Pointer, ProcessBaseType>
	(m,"SplitElementsProcess")
	.def(init<ModelPart&, int>());

      class_<SetActiveFlagProcess, SetActiveFlagProcess::Pointer, MesherProcess>
	(m, "SetActiveFlagProcess")
	.def(init<ModelPart&, bool, bool, int>());

     class_<SetActiveFlagMesherProcess, SetActiveFlagMesherProcess::Pointer, SetActiveFlagProcess>
	(m, "SetActiveFlagMesherProcess")
	.def(init<ModelPart&, bool, bool, int>());

      
      class_<AdaptiveTimeIntervalProcess, AdaptiveTimeIntervalProcess::Pointer, ProcessBaseType>
      	(m, "AdaptiveTimeIntervalProcess")
	.def(init<ModelPart&, int>());

     class_<ModelStartEndMeshingForFluidsProcess, ModelStartEndMeshingForFluidsProcess::Pointer, ModelStartEndMeshingProcessType>
       (m, "ModelMeshingForFluids")
       .def(init<ModelPart&, Flags, int>());

      //**********TRANSFER ELEMENTS TO MODEL PART*********//

      class_<TransferModelPartElementsProcess, TransferModelPartElementsProcess::Pointer, ProcessBaseType>
      	(m, "TransferModelPartElementsProcess")
	  .def(init<ModelPart&, ModelPart&>())
        .def("Execute", &TransferModelPartElementsProcess::Execute)
      	;


    }


  }  // namespace Python.

} // Namespace Kratos

