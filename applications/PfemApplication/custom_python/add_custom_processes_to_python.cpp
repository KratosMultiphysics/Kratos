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
#include "custom_processes/elemental_neighbours_search_process.hpp"
#include "custom_processes/nodal_neighbours_search_process.hpp"
#include "custom_processes/build_model_part_boundary_process.hpp"
#include "custom_processes/model_volume_calculation_process.hpp"

// MeshModeler initialization and finalization processes
#include "custom_processes/model_start_end_meshing_process.hpp"

// PreMeshing processes
#include "custom_processes/refine_mesh_elements_on_threshold_process.hpp"
#include "custom_processes/refine_mesh_elements_in_edges_process.hpp"
#include "custom_processes/refine_mesh_boundary_process.hpp"
#include "custom_processes/remove_mesh_nodes_process.hpp"

// MiddleMeshing processes
#include "custom_processes/refine_mesh_elements_on_size_process.hpp"
#include "custom_processes/print_output_mesh_process.hpp"

// PostMeshing processes
#include "custom_processes/generate_new_nodes_process.hpp"
#include "custom_processes/select_mesh_elements_process.hpp"
#include "custom_processes/build_mesh_elements_process.hpp"
#include "custom_processes/build_mesh_boundary_process.hpp"

// Kinematics
#include "custom_processes/constant_rotation_process.h"


namespace Kratos
{
	
namespace Python
{

typedef Process                        ProcessBaseType;
typedef Process::Pointer                ProcessPointer;
typedef std::vector<Process::Pointer> ProcessContainer;

void Push_Back_Process( ProcessContainer& ThisProcessContainer,
                        ProcessPointer ThisProcess )
{
  ThisProcessContainer.push_back( ThisProcess );
}
  	
void  AddCustomProcessesToPython(pybind11::module& m)
{

  using namespace pybind11;

  //process container
  class_<ProcessContainer>(m,"ProcessContainer")
      .def( init<>() )
      .def( "PushBack", Push_Back_Process )
      ;


  //***************NEIGHBOURS**************//
      
  class_<NodalNeighboursSearchProcess, typename NodalNeighboursSearchProcess::Pointer, Process>
      (m,"NodalNeighboursSearch")
      .def(init<ModelPart&, int, int, int>())
      .def("CleanNeighbours", &NodalNeighboursSearchProcess::ClearNeighbours)
      ;
      
  class_<ElementalNeighboursSearchProcess, typename ElementalNeighboursSearchProcess::Pointer, Process>
      (m,"ElementalNeighboursSearch")
      .def(init<ModelPart&, int, int, int>())
      .def("CleanNeighbours", &ElementalNeighboursSearchProcess::ClearNeighbours)
      ;


  //***************BOUNDARY**************//

  class_<BuildModelPartBoundaryProcess, typename BuildModelPartBoundaryProcess::Pointer, Process>
      (m,"BuildModelPartBoundary")
      .def(init<ModelPart&, std::string, int>())
      .def("SearchConditionMasters", &BuildModelPartBoundaryProcess::SearchConditionMasters)
      ;


  //**********MESH MODELLER PROCESS*********//

  class_<ModelStartEndMeshingProcess, typename ModelStartEndMeshingProcess::Pointer, Process>
      (m,"ModelMeshing")
      .def(init<ModelPart&, Flags, int>())
      ;


  class_<RefineMeshElementsOnThresholdProcess, typename RefineMeshElementsOnThresholdProcess::Pointer, Process>
      (m,"SetElementNodesToRefineOnThreshold")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>())
      ;

  class_<RefineMeshElementsInEdgesProcess, typename RefineMeshElementsInEdgesProcess::Pointer, Process>
      (m,"SetElementEdgesToRefine")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>())
      ;
      
  class_<RefineMeshElementsOnSizeProcess, typename RefineMeshElementsOnSizeProcess::Pointer, Process>
      (m,"SetElementsToRefineOnSize")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>())
      ;

  class_<RefineMeshBoundaryProcess, typename RefineMeshBoundaryProcess::Pointer, Process>
      (m,"RefineMeshBoundary")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>())
      ;

  class_<RemoveMeshNodesProcess, typename RemoveMeshNodesProcess::Pointer, Process>
      (m,"RemoveMeshNodes")
      .def(init<ModelPart&, ModelerUtilities::MeshingParameters&, int>())
      ;


  class_<GenerateNewNodesProcess, typename GenerateNewNodesProcess::Pointer, Process>
      (m,"GenerateNewNodes")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>())
      ;

  class_<SelectMeshElementsProcess, typename SelectMeshElementsProcess::Pointer, Process>
      (m,"SelectMeshElements")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>())
      ;

  class_<BuildMeshElementsProcess, typename BuildMeshElementsProcess::Pointer, Process>
      (m,"BuildMeshElements")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>())
      ;


  class_<BuildMeshBoundaryProcess, typename BuildMeshBoundaryProcess::Pointer, BuildModelPartBoundaryProcess>
      (m,"BuildMeshBoundary")
      .def(init<ModelPart&, ModelerUtilities::MeshingParameters&, int>())
      ;


  class_<PrintOutputMeshProcess, typename PrintOutputMeshProcess::Pointer, Process>
      (m,"PrintOutputMeshProcess")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, std::string, int>())
      ;
      

  //********MODEL VOLUME CALCULATION*********//

  class_<ModelVolumeCalculationProcess, typename ModelVolumeCalculationProcess::Pointer, Process>
      (m,"ModelVolumeCalculation")
      .def(init<ModelPart&, bool, int>())
      .def("ExecuteInitializeSolutionStep", &ModelVolumeCalculationProcess::ExecuteInitializeSolutionStep)
      .def("ExecuteFinalizeSolutionStep", &ModelVolumeCalculationProcess::ExecuteFinalizeSolutionStep)
      ;
      
  //********MODEL VOLUME CALCULATION*********//

  class_<ConstantRotationProcess, typename ConstantRotationProcess::Pointer, Process>
      (m,"ConstantRotationProcess")
      .def(init<ModelPart&, const double, const double, const double, const double, const double, const double>())	 
      .def(init< ModelPart&, Parameters& >())
      ;
      
      
}
 
}  // namespace Python.

} // Namespace Kratos

