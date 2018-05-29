//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes 

// External includes 

// Project includes
#include "custom_python/add_custom_processes_to_python.h"

// Processes
#include "custom_processes/contact_model_start_end_meshing_process.hpp"
#include "custom_processes/parametric_wall_contact_search_process.hpp"
#include "custom_processes/build_contact_model_part_process.hpp"
#include "custom_processes/clear_contact_conditions_process.hpp"
#include "custom_processes/clear_point_contact_conditions_process.hpp"
#include "custom_processes/build_contact_conditions_process.hpp"

#include "custom_processes/hm_parametric_wall_contact_search_process.hpp"

namespace Kratos
{
	
namespace Python
{
    
void Push_Back_String( std::vector<std::string>& ThisStringVector, std::string ThisString)
{
  ThisStringVector.push_back(ThisString);
}


void  AddCustomProcessesToPython(pybind11::module& m)
{

  using namespace pybind11;

  class_< std::vector<std::string> >(m,"StringVector")
      .def(init<>())
      .def("PushBack", Push_Back_String)
      ;

  //**********MESH MODELLER PROCESS*********//

  class_<ContactModelStartEndMeshingProcess, ContactModelStartEndMeshingProcess::Pointer, ModelStartEndMeshingProcess>
      (m, "ContactModelMeshing")
      .def(init<ModelPart&, Flags, int>())
      ;

  class_<ParametricWallContactSearchProcess, ParametricWallContactSearchProcess::Pointer, Process>
      (m,"ParametricWallContactSearch")
      .def(init<ModelPart&, std::string, SpatialBoundingBox::Pointer, Parameters>())
      ;

  class_<HMParametricWallContactSearchProcess, HMParametricWallContactSearchProcess::Pointer, Process> 
	   (m,"HMParametricWallContactSearch")
      .def(init<ModelPart&, std::string, SpatialBoundingBox::Pointer, Parameters>())
	   ;

  class_<BuildContactModelPartProcess, BuildContactModelPartProcess::Pointer, Process>
      (m,"BuildContactModelPart")
      .def(init<ModelPart&, ModelerUtilities::MeshingParameters&, std::vector<std::string>&, int>())
      ;

  class_<ClearPointContactConditionsProcess, ClearPointContactConditionsProcess::Pointer, Process>
      (m,"ClearPointContactConditions")
      .def(init<ModelPart&, int>())
      ;

  class_<ClearContactConditionsProcess, ClearContactConditionsProcess::Pointer, Process>
      (m,"ClearContactConditions")
      .def(init<ModelPart&, int>())
      ;

  class_<BuildContactConditionsProcess, BuildContactConditionsProcess::Pointer, Process>
      (m,"BuildContactConditions")
      .def(init<ModelPart&,  ModelerUtilities::MeshingParameters&, int>())
      ;

}
 
}  // namespace Python.

} // Namespace Kratos

