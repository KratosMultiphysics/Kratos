#pragma once

// --- External Includes ---
#include <pybind11/pybind11.h>

// --- Core Includes ---
#include "includes/define_python.h"


namespace Kratos::Python {


void AddCustomOperationsToPython(pybind11::module& rModule);


} // namespace Kratos::Python
