//
//   Project Name:        KratosSolversApplication    $
//   Created by:          $Author:        JMCarbonell $
//   Last modified by:    $Co-Author:                 $
//   Date:                $Date:         January 2019 $
//   Revision:            $Revision:              0.0 $
//
//

// System includes

#if defined(KRATOS_PYTHON)

// External includes

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"

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

}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
