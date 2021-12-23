//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//

#if !defined(KRATOS_ADD_PROCESSES_TO_PYTHON_H_INCLUDED)
#define KRATOS_ADD_PROCESSES_TO_PYTHON_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define_python.h"

namespace Kratos {
namespace Python {

void  AddCustomProcessesToPython(pybind11::module& m);

} // namespace Python
} // namespace Kratos

#endif // !defined(KRATOS_ADD_PROCESSES_TO_PYTHON_H_INCLUDED)
