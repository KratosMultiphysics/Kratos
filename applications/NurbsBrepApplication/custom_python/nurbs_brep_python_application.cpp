//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos {
namespace Python {

using namespace pybind11;

PYBIND11_MODULE(KratosNurbsBrepApplication,m)
{

    class_<KratosNurbsBrepApplication, KratosNurbsBrepApplication::Pointer, 
        KratosApplication>(m, "KratosNurbsBrepApplication")
        .def(init<>())
    ;

    AddCustomUtilitiesToPython(m);
}

} // namespace Python
} // namespace Kratos

#endif // defined(KRATOS_PYTHON)
