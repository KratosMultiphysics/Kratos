//
//  Main authors:    Miguel Angel Celigueta   maceli@cimne.upc.edu
//
//

// System includes

// External includes

// Project includes

#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/post_utilities.h"

#include "includes/define.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

void AddCustomUtilitiesToPython(pybind11::module &m)
{

    py::class_<PostUtilities, PostUtilities::Pointer>(m, "PostUtilities", py::module_local())
        .def(py::init<>())
        // .def("ComputeCurvatureOfBeamSolids", &PostUtilities::ComputeCurvatureOfBeamSolids)
        // .def("ComputeCurvatureOfBeam", &PostUtilities::ComputeCurvatureOfBeam)
        // .def("CreateSkinForBeam", &PostUtilities::CreateSkinForBeam)
        ;
}

} // namespace Python.

} // Namespace Kratos
