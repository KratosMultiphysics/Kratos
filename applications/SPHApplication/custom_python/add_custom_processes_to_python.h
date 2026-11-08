#pragma once

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

namespace Kratos::Python
{
    void  AddCustomProcessesToPython(pybind11::module& m);
}



