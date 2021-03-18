//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//


// System includes

// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "custom_utilities/drag_and_moment_utilities.h"


namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Add drag and moment utilities
    py::class_< DragAndMomentUtilities> (m,"DragAndMomentUtilities")
        .def(py::init<>())
        .def("CalculateBodyFittedDragAndMoment", &DragAndMomentUtilities::CalculateBodyFittedDragAndMoment)
        .def("CalculateEmbeddedDragAndMoment", &DragAndMomentUtilities::CalculateEmbeddedDragAndMoment)
        ;

}

} // namespace Python.
} // Namespace Kratos
