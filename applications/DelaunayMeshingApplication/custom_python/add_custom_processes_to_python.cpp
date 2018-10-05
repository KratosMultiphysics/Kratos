//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
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
#include "custom_processes/constant_rotation_process.hpp"

// Mesher initialization and finalization processes
#include "custom_processes/settle_model_structure_process.hpp"

// Mesher processes:

// PreMeshing processes
#include "custom_processes/refine_elements_on_threshold_mesher_process.hpp"
#include "custom_processes/refine_elements_in_edges_mesher_process.hpp"
#include "custom_processes/refine_conditions_mesher_process.hpp"
#include "custom_processes/remove_nodes_mesher_process.hpp"

// MiddleMeshing processes
#include "custom_processes/refine_elements_on_size_mesher_process.hpp"
#include "custom_processes/print_mesh_output_mesher_process.hpp"

// PostMeshing processes
#include "custom_processes/generate_new_nodes_mesher_process.hpp"
#include "custom_processes/select_elements_mesher_process.hpp"
#include "custom_processes/generate_new_elements_mesher_process.hpp"
#include "custom_processes/generate_new_conditions_mesher_process.hpp"


namespace Kratos
{

namespace Python
{

typedef Process::Pointer                           ProcessPointer;
typedef MesherProcess::Pointer               MesherProcessPointer;
typedef std::vector<MesherProcessPointer>  MesherProcessContainer;

void Push_Back_Process( MesherProcessContainer& ThisProcessContainer,
                        MesherProcessPointer ThisProcess )
{
  ThisProcessContainer.push_back( ThisProcess );
}

void  AddCustomProcessesToPython(pybind11::module& m)
{

  using namespace pybind11;

  //**********MESHER PROCESS*********//

  //mesher process container
  class_<MesherProcessContainer>(m,"MesherProcessContainer")
      .def( init<>() )
      .def( "PushBack", Push_Back_Process )
      ;

  class_<MesherProcess, MesherProcess::Pointer, Process>(m,"MesherProcess")
      .def(init<>())
      ;

  //***************NEIGHBOURS**************//

  class_<NodalNeighboursSearchProcess, NodalNeighboursSearchProcess::Pointer, MesherProcess>
      (m,"NodalNeighboursSearch")
      .def(init<ModelPart&, int, int, int>())
      .def("CleanNeighbours", &NodalNeighboursSearchProcess::ClearNeighbours)
      ;

  class_<ElementalNeighboursSearchProcess, ElementalNeighboursSearchProcess::Pointer, MesherProcess>
      (m,"ElementalNeighboursSearch")
      .def(init<ModelPart&, int, int, int>())
      .def("CleanNeighbours", &ElementalNeighboursSearchProcess::ClearNeighbours)
      ;


  //***************BOUNDARY**************//

  class_<BuildModelPartBoundaryProcess, BuildModelPartBoundaryProcess::Pointer, MesherProcess>
      (m,"BuildModelPartBoundary")
      .def(init<ModelPart&, std::string, int>())
      .def("SearchConditionMasters", &BuildModelPartBoundaryProcess::SearchConditionMasters)
      ;


  //**********MODEL STRUCTURE*********//

  class_<SettleModelStructureProcess, SettleModelStructureProcess::Pointer, Process>
      (m,"ModelStructure")
      .def(init<ModelPart&, Flags, int>())
      .def("ExecuteInitialize", &SettleModelStructureProcess::ExecuteInitialize)
      .def("ExecuteFinalize", &SettleModelStructureProcess::ExecuteFinalize)
      ;


  //**********MESHER PROCESSES*********//


  class_<RefineElementsOnThresholdMesherProcess, RefineElementsOnThresholdMesherProcess::Pointer, MesherProcess>
      (m,"RefineElementsOnThreshold")
      .def(init<ModelPart&, MesherUtilities::MeshingParameters&, int>())
      ;

  class_<RefineElementsOnSizeMesherProcess, RefineElementsOnSizeMesherProcess::Pointer, MesherProcess>
      (m,"RefineElementsOnSize")
      .def(init<ModelPart&,  MesherUtilities::MeshingParameters&, int>())
      ;

  class_<RefineElementsInEdgesMesherProcess, RefineElementsInEdgesMesherProcess::Pointer, MesherProcess>
      (m,"RefineElementsInEdges")
      .def(init<ModelPart&, MesherUtilities::MeshingParameters&, int>())
      ;

  class_<RefineConditionsMesherProcess, RefineConditionsMesherProcess::Pointer, MesherProcess>
      (m,"RefineConditions")
      .def(init<ModelPart&,  MesherUtilities::MeshingParameters&, int>())
      ;

  class_<RemoveNodesMesherProcess, RemoveNodesMesherProcess::Pointer, MesherProcess>
      (m,"RemoveNodes")
      .def(init<ModelPart&, MesherUtilities::MeshingParameters&, int>())
      ;

  class_<GenerateNewNodesMesherProcess, GenerateNewNodesMesherProcess::Pointer, MesherProcess>
      (m,"GenerateNewNodes")
      .def(init<ModelPart&,  MesherUtilities::MeshingParameters&, int>())
      ;

  class_<SelectElementsMesherProcess, SelectElementsMesherProcess::Pointer, MesherProcess>
      (m,"SelectElements")
      .def(init<ModelPart&,  MesherUtilities::MeshingParameters&, int>())
      ;

  class_<GenerateNewElementsMesherProcess, GenerateNewElementsMesherProcess::Pointer, MesherProcess>
      (m,"GenerateNewElements")
      .def(init<ModelPart&,  MesherUtilities::MeshingParameters&, int>())
      ;

  class_<GenerateNewConditionsMesherProcess, GenerateNewConditionsMesherProcess::Pointer, BuildModelPartBoundaryProcess>
      (m,"GenerateNewConditions")
      .def(init<ModelPart&, MesherUtilities::MeshingParameters&, int>())
      ;

  class_<PrintMeshOutputMesherProcess, PrintMeshOutputMesherProcess::Pointer, MesherProcess>
      (m,"PrintMeshOutput")
      .def(init<ModelPart&,  MesherUtilities::MeshingParameters&, std::string, int>())
      ;


  //********MODEL VOLUME CALCULATION*********//

  class_<ModelVolumeCalculationProcess, ModelVolumeCalculationProcess::Pointer, Process>
      (m,"ModelVolumeCalculation")
      .def(init<ModelPart&, bool, int>())
      .def("ExecuteInitializeSolutionStep", &ModelVolumeCalculationProcess::ExecuteInitializeSolutionStep)
      .def("ExecuteFinalizeSolutionStep", &ModelVolumeCalculationProcess::ExecuteFinalizeSolutionStep)
      ;

  //********CONSTANT ROTATION CALCULATION*********//

  class_<ConstantRotationProcess, ConstantRotationProcess::Pointer, Process>
      (m,"ConstantRotationProcess")
      .def(init<ModelPart&, const double, const double, const double, const double, const double, const double>())
      .def(init< ModelPart&, Parameters& >())
      ;


}

}  // namespace Python.

} // Namespace Kratos
