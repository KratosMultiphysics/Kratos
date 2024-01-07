/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Author: Thomas Oberbichler
*/

#if defined(KRATOS_PYTHON)

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define_python.h"
#include "linear_solvers_application.h"
#include "custom_python/add_custom_solvers_to_python.h"
#include "custom_python/add_custom_decompositions_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosLinearSolversApplication, m)
{
    namespace py = pybind11;

    py::class_<KratosLinearSolversApplication,
           KratosLinearSolversApplication::Pointer,
           KratosApplication>(m, "KratosLinearSolversApplication")
        .def(py::init<>())
        ;

    AddCustomSolversToPython(m);
    AddCustomDecompositionsToPython(m);
    AddCustomProcessesToPython(m);

    m.def("HasMKL", []() {
#if defined(USE_EIGEN_MKL)
        return true;
#else
        return false;
#endif
        });

    m.def("HasFEAST", []() {
#if defined(USE_EIGEN_FEAST)
        return true;
#else
        return false;
#endif
        });
}

} // namespace Python
} // namespace Kratos

#endif // defined(KRATOS_PYTHON)
