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

// Properties
#include "includes/properties.h"

// Processes
//#include "custom_processes/adaptive_time_interval_process.hpp"
#include "custom_processes/assign_properties_to_nodes_process.hpp"
#include "custom_processes/manage_isolated_nodes_process.hpp"
#include "custom_processes/manage_selected_elements_process.hpp"
#include "custom_processes/recover_volume_losses_process.hpp"

// PreMeshing processes
#include "custom_processes/inlet_mesher_process.hpp"
#include "custom_processes/insert_fluid_nodes_mesher_process.hpp"
#include "custom_processes/remove_fluid_nodes_mesher_process.hpp"
#include "custom_processes/refine_fluid_elements_in_edges_mesher_process.hpp"

// MiddleMeshing processes

// PostMeshing processes


namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{

  namespace py = pybind11;

  //**********MODEL PROPERTIES*********//

  /// Properties container. A vector set of properties with their Id's as key.
  typedef PointerVectorSet<Properties, IndexedObject> PropertiesContainerType;
  typedef typename PropertiesContainerType::Pointer   PropertiesContainerPointerType;

  //to define it as a variable
  py::class_<Variable<PropertiesContainerPointerType>, VariableData>(m,"PropertiesVectorPointerVariable")
      .def( "__repr__", &Variable<PropertiesContainerPointerType>::Info )
      ;


  //**********MESHER PROCESSES*********//

  py::class_<RemoveFluidNodesMesherProcess, RemoveFluidNodesMesherProcess::Pointer, RemoveNodesMesherProcess>
      (m, "RemoveFluidNodes")
      .def(py::init<ModelPart&, MesherUtilities::MeshingParameters&, int>());

  py::class_<InsertFluidNodesMesherProcess, InsertFluidNodesMesherProcess::Pointer, MesherProcess>
      (m, "InsertFluidNodes")
      .def(py::init<ModelPart&, MesherUtilities::MeshingParameters&, int>());

  py::class_<InletMesherProcess, InletMesherProcess::Pointer, MesherProcess>
      (m, "InsertInlet")
      .def(py::init<ModelPart&, MesherUtilities::MeshingParameters&, int>());

  py::class_<RefineFluidElementsInEdgesMesherProcess, RefineFluidElementsInEdgesMesherProcess::Pointer, RefineElementsInEdgesMesherProcess>
      (m,"RefineFluidElementsInEdges")
      .def(py::init<ModelPart&, MesherUtilities::MeshingParameters&, int>())
      ;

  //*********SET SOLVER PROCESSES*************//
  py::class_<AssignPropertiesToNodesProcess, AssignPropertiesToNodesProcess::Pointer, Process>
      (m, "AssignPropertiesToNodes")
      .def(py::init<ModelPart&, Parameters>())
      .def(py::init<ModelPart&, Parameters&>());

  //*********ADAPTIVE TIME STEP*************//
  // py::class_<AdaptiveTimeIntervalProcess, AdaptiveTimeIntervalProcess::Pointer, Process>
  //     (m, "AdaptiveTimeIntervalProcess")
  //     .def(py::init<ModelPart&, int>());

  //*********VOLUME RECOVERY PROCESS********//
  py::class_<RecoverVolumeLossesProcess, RecoverVolumeLossesProcess::Pointer, Process>
      (m, "RecoverVolumeLosses")
      .def(py::init<ModelPart&, int>());

  //*********MANAGE PARTICULAR ENTITIES PROCESS********//
  py::class_<ManageIsolatedNodesProcess, ManageIsolatedNodesProcess::Pointer, Process>
      (m, "ManageIsolatedNodesProcess")
      .def(py::init<ModelPart&>());

  py::class_<ManageSelectedElementsProcess, ManageSelectedElementsProcess::Pointer, Process>
      (m, "ManageSelectedElementsProcess")
      .def(py::init<ModelPart&>());

}

}  // namespace Python.

} // Namespace Kratos
