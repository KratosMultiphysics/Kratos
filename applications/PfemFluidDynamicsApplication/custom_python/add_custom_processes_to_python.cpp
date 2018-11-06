//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2018 $
//   Revision:            $Revision:                     0.0 $
//
//

// External includes


//Application includes
#include "custom_python/add_custom_processes_to_python.h"

// Project includes
#include "includes/node.h"
#include "processes/process.h"

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
#include "custom_processes/model_start_end_meshing_with_conditions_for_fluids_process.hpp"
#include "custom_processes/split_elements_process.hpp"
#include "custom_processes/set_active_flag_process.hpp"
#include "custom_processes/set_active_flag_mesher_process.hpp"
#include "custom_processes/set_material_properties_to_fluid_nodes_process.hpp"
#include "custom_processes/set_material_properties_from_fluid_to_rigid_nodes_process.hpp"
#include "custom_processes/set_material_properties_to_solid_nodes_process.hpp"
#include "custom_processes/adaptive_time_interval_process.hpp"
#include "custom_processes/transfer_model_part_elements_process.hpp"
#include "custom_processes/build_mesh_boundary_for_fluids_process.hpp"
#include "custom_processes/build_model_part_boundary_for_fluids_process.hpp"
#include "custom_processes/generate_new_conditions_mesher_for_fluids_process.hpp"

//Processes


namespace Kratos
{

  namespace Python
  {


    void  AddCustomProcessesToPython(pybind11::module& m)
    {

      namespace py = pybind11;
      typedef Process                                         ProcessBaseType;
      typedef SettleModelStructureProcess     ModelStartEndMeshingProcessType;



      py::class_<RecoverVolumeLossesProcess, RecoverVolumeLossesProcess::Pointer, MesherProcess>
	(m, "RecoverVolumeLosses")
	.def(py::init<ModelPart&,  MesherUtilities::MeshingParameters&, int>());

      py::class_<RemoveMeshNodesForFluidsProcess, RemoveMeshNodesForFluidsProcess::Pointer, MesherProcess>
      	(m, "RemoveMeshNodesForFluids")
	.def(py::init<ModelPart&, MesherUtilities::MeshingParameters&, int>());

      py::class_<GenerateNewNodesBeforeMeshingProcess, GenerateNewNodesBeforeMeshingProcess::Pointer, MesherProcess>
      	(m, "GenerateNewNodesBeforeMeshing")
	.def(py::init<ModelPart&,  MesherUtilities::MeshingParameters&, int>());

      py::class_<SelectMeshElementsForFluidsProcess, SelectMeshElementsForFluidsProcess::Pointer, MesherProcess>
	(m, "SelectMeshElementsForFluids")
	.def(py::init<ModelPart&,  MesherUtilities::MeshingParameters&, int>());

      py::class_<InletManagementProcess, InletManagementProcess::Pointer, MesherProcess>
      	(m, "InletManagement")
	.def(py::init<ModelPart&,  MesherUtilities::MeshingParameters&, int>());

      py::class_<SetInletProcess, SetInletProcess::Pointer, ProcessBaseType>
      	(m, "SetInlet")
	.def(py::init<ModelPart&, int>());

      py::class_<SplitElementsProcess, SplitElementsProcess::Pointer, ProcessBaseType>
	(m,"SplitElementsProcess")
	.def(py::init<ModelPart&, int>());

      py::class_<SetActiveFlagProcess, SetActiveFlagProcess::Pointer, MesherProcess>
	(m, "SetActiveFlagProcess")
	.def(py::init<ModelPart&, bool, bool, int>());

      py::class_<SetMaterialPropertiesToFluidNodesProcess, SetMaterialPropertiesToFluidNodesProcess::Pointer, MesherProcess>
	(m, "SetMaterialPropertiesToFluidNodes")
	.def(py::init<ModelPart&>());

      py::class_<SetMaterialPropertiesFromFluidToRigidNodesProcess, SetMaterialPropertiesFromFluidToRigidNodesProcess::Pointer, MesherProcess>
	(m, "SetMaterialPropertiesFromFluidToRigidNodes")
	.def(py::init<ModelPart&, ModelPart&>());

      py::class_<SetMaterialPropertiesToSolidNodesProcess, SetMaterialPropertiesToSolidNodesProcess::Pointer, MesherProcess>
	(m, "SetMaterialPropertiesToSolidNodes")
	.def(py::init<ModelPart&>());

     py::class_<SetActiveFlagMesherProcess, SetActiveFlagMesherProcess::Pointer, SetActiveFlagProcess>
	(m, "SetActiveFlagMesherProcess")
	.def(py::init<ModelPart&, bool, bool, int>());


      py::class_<AdaptiveTimeIntervalProcess, AdaptiveTimeIntervalProcess::Pointer, ProcessBaseType>
      	(m, "AdaptiveTimeIntervalProcess")
	.def(py::init<ModelPart&, int>());

     py::class_<ModelStartEndMeshingWithConditionsForFluidsProcess, ModelStartEndMeshingWithConditionsForFluidsProcess::Pointer, ModelStartEndMeshingProcessType>
       (m, "ModelMeshingWithConditionsForFluids")
       .def(py::init<ModelPart&, Flags, int>());

     py::class_<ModelStartEndMeshingForFluidsProcess, ModelStartEndMeshingForFluidsProcess::Pointer, ModelStartEndMeshingProcessType>
       (m, "ModelMeshingForFluids")
       .def(py::init<ModelPart&, Flags, int>());

     py::class_<BuildMeshBoundaryForFluidsProcess, BuildMeshBoundaryForFluidsProcess::Pointer, MesherProcess>
	(m, "BuildMeshBoundaryForFluids")
	.def(py::init<ModelPart&, MesherUtilities::MeshingParameters&, int>());

     py::class_<BuildModelPartBoundaryForFluidsProcess, BuildModelPartBoundaryForFluidsProcess::Pointer, MesherProcess>
	(m, "BuildModelPartBoundaryForFluids")
       .def(py::init<ModelPart&, std::string, int>())
       .def("SearchConditionMasters", &BuildModelPartBoundaryForFluidsProcess::SearchConditionMasters)
       ;
      //**********TRANSFER ELEMENTS TO MODEL PART*********//

      py::class_<TransferModelPartElementsProcess, TransferModelPartElementsProcess::Pointer, ProcessBaseType>
      	(m, "TransferModelPartElementsProcess")
	  .def(py::init<ModelPart&, ModelPart&>())
        .def("Execute", &TransferModelPartElementsProcess::Execute)
      	;

        py::class_<GenerateNewConditionsMesherForFluidsProcess, GenerateNewConditionsMesherForFluidsProcess::Pointer, BuildModelPartBoundaryProcess>
      (m,"GenerateNewConditionsForFluids")
      .def(py::init<ModelPart&, MesherUtilities::MeshingParameters&, int>())
      ;



    }


  }  // namespace Python.

} // Namespace Kratos

