//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_processes_to_python.h"

//Processes
#include "custom_processes/refine_conditions_in_contact_mesher_process.hpp"
#include "custom_processes/set_mechanical_initial_state_process.hpp"

namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{

  namespace py = pybind11;

  typedef std::vector<SpatialBoundingBox::Pointer>   BoundingBoxContainer;

  // Mesher process
  py::class_<RefineConditionsInContactMesherProcess, RefineConditionsInContactMesherProcess::Pointer, RefineConditionsMesherProcess>(m,"RefineConditionsInContact")
      .def(py::init<ModelPart&, BoundingBoxContainer&, MesherUtilities::MeshingParameters&, int>())
      ;

  // Set initial mechanical state
  py::class_<SetMechanicalInitialStateProcess, SetMechanicalInitialStateProcess::Pointer, Process>(m,"SetMechanicalInitialStateProcess")
      .def(py::init<ModelPart&, Parameters>())
      .def(py::init<ModelPart&, Parameters>())
      .def("Execute", &SetMechanicalInitialStateProcess::Execute)
      ;
}

}  // namespace Python.

} // Namespace Kratos
