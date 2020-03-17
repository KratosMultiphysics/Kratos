/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Armin Geiser
*/

#if !defined(KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED)
#define KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes

namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(pybind11::module& m);

} // namespace Python
} // namespace Kratos

#endif // defined(KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED)
