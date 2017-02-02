//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <vector>

// External includes 

// Project includes
#include "includes/node.h"
#include "processes/process.h"

//Application includes
#include "custom_python/add_custom_processes_to_python.h"

//Processes
#include "custom_processes/contact_model_start_end_meshing_process.hpp"
#include "custom_processes/parametric_wall_contact_search_process.hpp"
#include "custom_processes/build_contact_model_part_process.hpp"
#include "custom_processes/clear_contact_conditions_process.hpp"
#include "custom_processes/clear_point_contact_conditions_process.hpp"
#include "custom_processes/build_contact_conditions_process.hpp"


namespace Kratos
{
	
  namespace Python
  {
    
    void Push_Back_String( std::vector<std::string>& ThisStringVector, std::string ThisString)
    {
      ThisStringVector.push_back(ThisString);
    }


    void  AddCustomProcessesToPython()
    {

      using namespace boost::python;
      typedef Process                                              ProcessBaseType;
      typedef ModelStartEndMeshingProcess      ModelStartEndMeshingProcessBaseType;


      class_< std::vector<std::string> >("StringVector",init<>())
	.def("PushBack", Push_Back_String)
	;

      //**********MESH MODELLER PROCESS*********//

      class_<ContactModelStartEndMeshingProcess, bases< ModelStartEndMeshingProcessBaseType >, boost::noncopyable >
	("ContactModelMeshing", init<ModelPart&, Flags, int>())
	;

      class_<ParametricWallContactSearchProcess, bases< ProcessBaseType >, boost::noncopyable >
	("ParametricWallContactSearch", init<ModelPart&, std::string, SpatialBoundingBox::Pointer, Parameters>())
	;

      class_<BuildContactModelPartProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "BuildContactModelPart", init<ModelPart&, ModelerUtilities::MeshingParameters&, std::vector<std::string>&, int>()
	 )
	;

      class_<ClearPointContactConditionsProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "ClearPointContactConditions", init<ModelPart&, int>()
	 )
	;

      class_<ClearContactConditionsProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "ClearContactConditions", init<ModelPart&, int>()
	 )
	;

      class_<BuildContactConditionsProcess, bases<ProcessBaseType>, boost::noncopyable >
	(
	 "BuildContactConditions", init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>()
	 )
	;

    }
 
  }  // namespace Python.

} // Namespace Kratos

