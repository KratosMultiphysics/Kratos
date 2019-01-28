//
//   Project Name:        KratosSolversApplication $
//   Created by:          $Author:     JMCarbonell $
//   Last modified by:    $Co-Author:              $
//   Date:                $Date:      January 2019 $
//   Revision:            $Revision:           0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_processes_to_python.h"

// Processes
#include "custom_processes/add_dofs_process.hpp"

// Solver Processes
#include "custom_processes/solver_process.hpp"

namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{

  namespace py = pybind11;


  py::class_<AddDofsProcess, AddDofsProcess::Pointer, Process>(m,"AddDofsProcess")
      .def(py::init<ModelPart&, Parameters>())
      .def(py::init<ModelPart&, Parameters&>())
      .def(py::init<ModelPart&, const pybind11::list&, const pybind11::list&>())
      .def("Execute", &AddDofsProcess::Execute)
      ;

  py::class_<SolverProcess, SolverProcess::Pointer, Process>(m,"SolverProcess")
      .def(py::init<>())
      ;


}

}  // namespace Python.

} // Namespace Kratos
