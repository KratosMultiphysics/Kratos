/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Thomas Oberbichler
*/

#if defined(KRATOS_PYTHON)

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define_python.h"
#include "eigen_solvers_application.h"
#include "custom_python/add_custom_solvers_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

PYBIND11_MODULE(KratosEigenSolversApplication, m)
{

    class_<KratosEigenSolversApplication,
           KratosEigenSolversApplication::Pointer,
           KratosApplication>(m, "KratosEigenSolversApplication")
        .def(init<>())
        ;

    AddCustomSolversToPython(m);
}

} // namespace Python

} // namespace Kratos

#endif // defined(KRATOS_PYTHON)
