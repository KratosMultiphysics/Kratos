//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// // ------------------------------------------------------------------------------
// // External includes
// // ------------------------------------------------------------------------------
#include <pybind11/stl.h>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/geometrical/symmetry_utility.h"

// ==============================================================================

namespace Kratos {
namespace Python {



// ==============================================================================
void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // ================================================================
    // 
    // ================================================================
    py::class_<SymmetryUtility >(m, "SymmetryUtility")
        .def(py::init<std::string, ModelPart&, Parameters>())
        .def("Initialize", &SymmetryUtility::Initialize)
        .def("Update", &SymmetryUtility::Update)
        .def("ApplyOnVectorField", &SymmetryUtility::ApplyOnVectorField)
        .def("ApplyOnScalarField", &SymmetryUtility::ApplyOnScalarField)
        ;
}

}  // namespace Python.
} // Namespace Kratos

