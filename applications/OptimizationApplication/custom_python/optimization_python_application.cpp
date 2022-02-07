// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// =================================================================================

#if defined(KRATOS_PYTHON)

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <pybind11/pybind11.h>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define_python.h"
#include "optimization_application.h"
#include "optimization_application_variables.h"

// ==============================================================================

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosOptimizationApplication, m)
{
    namespace py = pybind11;

    py::class_<KratosOptimizationApplication,
        KratosOptimizationApplication::Pointer,
        KratosApplication >(m, "KratosOptimizationApplication")
        .def(py::init<>())
        ;

  }

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
