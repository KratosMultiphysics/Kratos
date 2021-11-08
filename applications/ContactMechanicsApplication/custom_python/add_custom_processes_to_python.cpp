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
#include "custom_processes/settle_contact_model_structure_process.hpp"
#include "custom_processes/build_contact_model_part_process.hpp"

#include "custom_processes/parametric_wall_contact_search_process.hpp"
#include "custom_processes/hm_parametric_wall_contact_search_process.hpp"
#include "custom_processes/clear_point_contact_conditions_process.hpp"

// Mesher processes:
#include "custom_processes/clear_contact_conditions_mesher_process.hpp"
#include "custom_processes/generate_new_contact_conditions_mesher_process.hpp"



namespace Kratos
{

namespace Python
{


void  AddCustomProcessesToPython(pybind11::module& m)
{

  namespace py = pybind11;


  //**********MODEL STRUCTURE*********//

  py::class_<SettleContactModelStructureProcess, SettleContactModelStructureProcess::Pointer, SettleModelStructureProcess>
      (m, "ContactModelStructure")
      .def(py::init<ModelPart&, Flags, int>())
      ;

  py::class_<BuildContactModelPartProcess, BuildContactModelPartProcess::Pointer, Process>
      (m,"BuildContactModelPart")
      .def(py::init<ModelPart&, MesherUtilities::MeshingParameters&, std::vector<std::string>&, int>())
      ;

  //**********CONTACT WITH PARAMETRIC WALLS*********//

  py::class_<ParametricWallContactSearchProcess, ParametricWallContactSearchProcess::Pointer, Process>
      (m,"ParametricWallContactSearch")
      .def(py::init<ModelPart&, std::string, SpatialBoundingBox::Pointer, Parameters>())
      ;

  py::class_<HMParametricWallContactSearchProcess, HMParametricWallContactSearchProcess::Pointer, Process>
      (m,"HMParametricWallContactSearch")
      .def(py::init<ModelPart&, std::string, SpatialBoundingBox::Pointer, Parameters>())
      ;

  py::class_<ClearPointContactConditionsProcess, ClearPointContactConditionsProcess::Pointer, Process>
      (m,"ClearPointContactConditions")
      .def(py::init<ModelPart&, int>())
      ;

  //**********MESHER PROCESSES*********//

  py::class_<ClearContactConditionsMesherProcess, ClearContactConditionsMesherProcess::Pointer, MesherProcess>
      (m,"ClearContactConditions")
      .def(py::init<ModelPart&, int>())
      ;

  py::class_<GenerateNewContactConditionsMesherProcess, GenerateNewContactConditionsMesherProcess::Pointer, MesherProcess>
      (m,"BuildContactConditions")
      .def(py::init<ModelPart&, MesherUtilities::MeshingParameters&, int>())
      ;

}

}  // namespace Python.

} // Namespace Kratos
