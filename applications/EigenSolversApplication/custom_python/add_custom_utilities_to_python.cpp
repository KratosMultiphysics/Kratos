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

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/make_consistent.h"

namespace Kratos {
namespace Python {

// ==============================================================================
void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;
    m.def("MakeConsistent", &MakeConsistent);
}

} // namespace Python

} // namespace Kratos
