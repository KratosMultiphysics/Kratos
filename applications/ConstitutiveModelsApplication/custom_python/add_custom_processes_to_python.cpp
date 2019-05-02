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
#include "custom_processes/non_local_plasticity_process.hpp"

namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{

  namespace py = pybind11;


  // Set initial mechanical state
  py::class_<NonLocalPlasticityProcess, NonLocalPlasticityProcess::Pointer, Process>(m,"NonLocalPlasticityProcess")
      .def(py::init<ModelPart&, Parameters>())
      .def(py::init<ModelPart&, Parameters>())
      .def("Execute", &NonLocalPlasticityProcess::Execute)
      ;
}

}  // namespace Python.

} // Namespace Kratos

