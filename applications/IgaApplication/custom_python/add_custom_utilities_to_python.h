/*
//  KRATOS .___  ________    _____
//         |   |/  _____/   /  _  \
//         |   /   \  ___  /  /_\  \
//         |   \    \_\  \/    |    \
//         |___|\______  /\____|__  /
//                     \/         \/  Application
//
//  License: BSD License
//           Kratos default license: kratos/license.txt
*/

#if !defined(KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED)
#define KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define_python.h"

namespace Kratos {
namespace Python {

void  AddCustomUtilitiesToPython(pybind11::module& m);

} // namespace Python
} // namespace Kratos

#endif // !defined(KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED)
