/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "iga_application_variables.h"

#include "spaces/ublas_space.h"
#include "custom_utilities/director_utilities.h"


namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(
    pybind11::module& m)
{
    pybind11::class_< DirectorUtilities >(m, "DirectorUtilities")
        .def(pybind11::init<ModelPart&, Parameters>())
        .def("ComputeDirectors",
            &DirectorUtilities::ComputeDirectors)
        ;
}

} // namespace Python
} // Namespace Kratos
