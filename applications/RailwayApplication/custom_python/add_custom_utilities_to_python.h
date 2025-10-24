

#if !defined(KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED)
#define KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define_python.h"

namespace Kratos::Python
{

void AddCustomUtilitiesToPython(const pybind11::module& rModule);

} // namespace Kratos::Python.

#endif // KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED  defined
