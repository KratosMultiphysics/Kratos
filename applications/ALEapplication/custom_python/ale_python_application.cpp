//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "ale_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/define.h"

namespace Kratos {

namespace Python {

using namespace pybind11;

PYBIND11_MODULE(KratosALEApplication,m) {

  class_<KratosALEApplication, KratosALEApplication::Pointer,
         KratosApplication>(m,"KratosALEApplication")
         .def(init<>());

  AddCustomStrategiesToPython(m);
  AddCustomUtilitiesToPython(m);

}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
