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

#if !defined(KRATOS_ADD_SOLVERS_TO_PYTHON_H_INCLUDED)
#define KRATOS_ADD_SOLVERS_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

// External includes
#if defined USE_FEAST4
#define _MKL_SOLVERS_EE_H //exclude MKL's FEAST implementation
#endif

// Project includes

namespace Kratos {
namespace Python {

void AddCustomSolversToPython(pybind11::module& m);

} // namespace Python
} // namespace Kratos

#endif // defined(KRATOS_ADD_SOLVERS_TO_PYTHON_H_INCLUDED)
