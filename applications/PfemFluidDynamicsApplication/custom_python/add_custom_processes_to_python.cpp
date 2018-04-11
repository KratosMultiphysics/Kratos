//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes 

// External includes 

// Project includes
#include "includes/model_part.h"
#include "custom_python/add_custom_processes_to_python.h"

// Processes

// Middle meshing processes
#include "custom_processes/remove_mesh_nodes_for_fluids_process.hpp"

// Post meshing processes
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

namespace Kratos
{
	
namespace Python
{


void  AddCustomProcessesToPython(pybind11::module& m)
{

  using namespace pybind11;    

  class_<RecoverVolumeLossesProcess, Process>(m,"RecoverVolumeLosses")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>())
      ;

  class_<RemoveMeshNodesForFluidsProcess, Process>(m,"RemoveMeshNodesForFluids")
      .def(init<ModelPart&, ModelerUtilities::MeshingParameters&, int>())
      ;

  class_<GenerateNewNodesBeforeMeshingProcess, Process>(m,"GenerateNewNodesBeforeMeshing")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>())
      ;

  class_<SelectMeshElementsForFluidsProcess, Process>(m,"SelectMeshElementsForFluids")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>())
      ;

  class_<InletManagementProcess, Process>(m,"InletManagement")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>())
      ;

  class_<SetInletProcess, Process>(m,"SetInlet")
      .def(init<ModelPart&, int>())
      ;

  class_<SplitElementsProcess, Process>(m,"SplitElementsProcess")
      .def(init<ModelPart&, int>())
      ;

  class_<SetActiveFlagProcess,  Process>(m,"SetActiveFlagProcess")
      .def(init<ModelPart&, bool, bool, int>())
      ;

  class_<AdaptiveTimeIntervalProcess, Process>(m,"AdaptiveTimeIntervalProcess")
      .def(init<ModelPart&, int>())
      ;

  class_<ModelStartEndMeshingForFluidsProcess, ModelStartEndMeshingProcess>(m,"ModelMeshingForFluids")
      .def(init<ModelPart&, Flags, int>())
      ;

  //**********TRANSFER ELEMENTS TO MODEL PART*********//

  class_<TransferModelPartElementsProcess, Process>(m,"TransferModelPartElementsProcess")
      .def(init<ModelPart&, ModelPart&>())
      .def("Execute", &TransferModelPartElementsProcess::Execute)
      ;
      
}	
 
 
}  // namespace Python.

} // Namespace Kratos

