/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Manuel Messmer
*/

#if !defined(KRATOS_ADD_SOLVERS_TO_PYTHON_H_INCLUDED)
#define KRATOS_ADD_SOLVERS_TO_PYTHON_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define_python.h"

namespace Kratos {
namespace Python {

void  AddCustomSolversToPython(pybind11::module& m);

} // namespace Python
} // namespace Kratos

#endif // KRATOS_ADD_SOLVERS_TO_PYTHON_H_INCLUDED defined