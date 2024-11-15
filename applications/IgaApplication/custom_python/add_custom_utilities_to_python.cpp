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
#include "custom_utilities/iga_flags.h"
#include "custom_utilities/array_1d_max.h"

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

    pybind11::class_< IgaFlags > iga_flags = pybind11::class_< IgaFlags >(m, "IgaFlags")
        .def(pybind11::init<>())
        ;
    pybind11::class_<Array1DMax>(m, "Array1DMax")
            .def_static("FindMaxInArray1D", &Array1DMax::FindMaxInArray1D, "Find the maximum value in an array_1d<double, 3>");     

    iga_flags.attr("FIX_DISPLACEMENT_X") = IgaFlags::FIX_DISPLACEMENT_X;
    iga_flags.attr("FIX_DISPLACEMENT_Y") = IgaFlags::FIX_DISPLACEMENT_Y;
    iga_flags.attr("FIX_DISPLACEMENT_Z") = IgaFlags::FIX_DISPLACEMENT_Z;
    iga_flags.attr("FIX_ROTATION_X") = IgaFlags::FIX_ROTATION_X;
    iga_flags.attr("FIX_ROTATION_Y") = IgaFlags::FIX_ROTATION_Y;
    iga_flags.attr("FIX_ROTATION_Z") = IgaFlags::FIX_ROTATION_Z;
}

} // namespace Python
} // Namespace Kratos
