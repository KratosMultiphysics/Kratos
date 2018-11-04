//
//   Project Name:        Kratos
//   Last modified by:    $Author: janosch $
//   Date:                $Date: 2008-07-23 14:47:45 $
//   Revision:            $Revision: 1.3 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "externalsolvers_application.h"
#include "custom_python/add_linear_solvers_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosExternalSolversApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosExternalSolversApplication,
        KratosExternalSolversApplication::Pointer,
        KratosApplication >(m,"KratosExternalSolversApplication")
        .def(py::init<>())
        ;

    AddLinearSolversToPython(m);
}

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
