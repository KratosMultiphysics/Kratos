//
//   Project Name:        KratosSolversApplication $
//   Developed by:        $Developer:  JMCarbonell $
//   Maintained by:       $Maintainer:        JMC  $
//   Date:                $Date:      January 2019 $
//
//

// System includes

#if defined(KRATOS_PYTHON)

// External includes

// Project includes

#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "solvers_application_variables.h"
#include "solvers_application.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

PYBIND11_MODULE(KratosSolversApplication,m)
{
  py::class_<KratosSolversApplication,
             KratosSolversApplication::Pointer,
             KratosApplication>(m, "KratosSolversApplication")
      .def(py::init<>())
      ;

  AddCustomStrategiesToPython(m);
  AddCustomProcessesToPython(m);
  AddCustomUtilitiesToPython(m);

  // Register python variables:

  // time settings
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MESHING_STEP_TIME )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_STEP_TIME )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RESTART_STEP_TIME )

  // time integration methods
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VECTOR_TIME_INTEGRATION_METHODS)
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COMPONENT_TIME_INTEGRATION_METHODS)
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_TIME_INTEGRATION_METHODS)

  // implicit solver
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TIME_INTEGRATION_ORDER)
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RAYLEIGH_ALPHA)
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RAYLEIGH_BETA)
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONVERGENCE_ACHIEVED)
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COMPUTE_CONSISTENT_MASS_MATRIX)

  // eigenvalue solver
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BUILD_LEVEL)
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EIGENVALUE_VECTOR)
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EIGENVECTOR_MATRIX)

}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
