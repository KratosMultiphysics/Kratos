//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
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

//General model processes

//Processes
#include "custom_processes/elemental_neighbours_search_process.hpp"
#include "custom_processes/nodal_neighbours_search_process.hpp"
#include "custom_processes/build_model_part_boundary_process.hpp"
#include "custom_processes/model_volume_calculation_process.hpp"

//MeshModeler initialization and finalization processes
#include "custom_processes/model_start_end_meshing_process.hpp"

//PreMeshing processes
#include "custom_processes/add_nodes_of_interest_process.hpp"
#include "custom_processes/refine_mesh_elements_on_threshold_process.hpp"
#include "custom_processes/refine_mesh_elements_in_edges_process.hpp"
#include "custom_processes/refine_mesh_boundary_process.hpp"
#include "custom_processes/remove_mesh_nodes_process.hpp"

//MiddleMeshing processes
#include "custom_processes/refine_mesh_elements_on_size_process.hpp"
#include "custom_processes/cpt_refine_mesh_elements_on_size_process.hpp"
#include "custom_processes/print_output_mesh_process.hpp"

//PostMeshing processes
#include "custom_processes/generate_new_nodes_process.hpp"
#include "custom_processes/select_mesh_elements_process.hpp"
#include "custom_processes/build_mesh_elements_process.hpp"
#include "custom_processes/build_mesh_boundary_process.hpp"

//Kinematics
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
  	
    void  AddCustomProcessesToPython()
    {

      using namespace boost::python;
      typedef Process                                         ProcessBaseType;


      //process container
      class_< ProcessContainer >( "ProcessContainer", init<>() )
        .def( "PushBack", Push_Back_Process )
      ;


     //***************NEIGHBOURS**************//
      
      class_<NodalNeighboursSearchProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "NodalNeighboursSearch", init<ModelPart&, int, int, int>()
	 )
	.def("CleanNeighbours", &NodalNeighboursSearchProcess::ClearNeighbours)
	;
      
      class_<ElementalNeighboursSearchProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "ElementalNeighboursSearch", init<ModelPart&, int, int, int>()
	 )
	.def("CleanNeighbours", &ElementalNeighboursSearchProcess::ClearNeighbours)
	;


      //***************BOUNDARY**************//

      class_<BuildModelPartBoundaryProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "BuildModelPartBoundary", init<ModelPart&, std::string, int>()
	 )
	.def("SearchConditionMasters", &BuildModelPartBoundaryProcess::SearchConditionMasters)
	;


      //**********MESH MODELLER PROCESS*********//

      class_<ModelStartEndMeshingProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "ModelMeshing", init<ModelPart&, Flags, int>()
	 )
	;

      class_<AddNodesOfInterestProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "AddNodesOfInterest", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>()
	 )
	;

      class_<RefineMeshElementsOnThresholdProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "SetElementNodesToRefineOnThreshold", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>()
	 )
	;

      class_<RefineMeshElementsInEdgesProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "SetElementEdgesToRefine", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>()
	 )
	;
      
      class_<RefineMeshElementsOnSizeProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "SetElementsToRefineOnSize", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>()
	 )
	;
      
      class_<CptRefineMeshElementsOnSizeProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "SetElementsToRefineOnSizeCPT", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>()
	 )
	;
      class_<RefineMeshBoundaryProcess, bases<ProcessBaseType>, boost::noncopyable >
      	(
      	 "RefineMeshBoundary", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>()
      	 )
      	;

      class_<RemoveMeshNodesProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "RemoveMeshNodes", init<ModelPart&, ModelerUtilities::MeshingParameters&, int>()
	 )
	;


      class_<GenerateNewNodesProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "GenerateNewNodes", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>()
	 )
	;

      class_<SelectMeshElementsProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "SelectMeshElements", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>()
	 )
	;

      class_<BuildMeshElementsProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "BuildMeshElements", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>()
	 )
	;


      class_<BuildMeshBoundaryProcess, bases<BuildModelPartBoundaryProcess>, boost::noncopyable >
	(
	 "BuildMeshBoundary", init<ModelPart&, ModelerUtilities::MeshingParameters&, int>()
	 )
	;


      class_<PrintOutputMeshProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "PrintOutputMeshProcess", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>()
	 )
	;
      

      //********MODEL VOLUME CALCULATION*********//

      class_<ModelVolumeCalculationProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "ModelVolumeCalculation", init<ModelPart&, bool, int>()
	 )
	 .def("ExecuteInitializeSolutionStep", &ModelVolumeCalculationProcess::ExecuteInitializeSolutionStep)
	 .def("ExecuteFinalizeSolutionStep", &ModelVolumeCalculationProcess::ExecuteFinalizeSolutionStep)
	;
      
      //********MODEL VOLUME CALCULATION*********//

      class_<ConstantRotationProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "ConstantRotationProcess", init<ModelPart&, const double, const double, const double, const double, const double, const double>()
	 )	 
         .def(init< ModelPart&, Parameters& >())
	;
      
      
    }
 
  }  // namespace Python.

} // Namespace Kratos

