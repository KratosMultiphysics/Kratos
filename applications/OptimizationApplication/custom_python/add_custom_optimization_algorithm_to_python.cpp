// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// =================================================================================

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
#include "custom_python/add_custom_optimization_algorithm_to_python.h"
#include "custom_algorithms/algorithm_steepest_descent.h"


// ==============================================================================

namespace Kratos {
namespace Python {



// ==============================================================================
void  AddCustomOptimizationAlgorithmToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // ================================================================
    // 
    // ================================================================
    py::class_<AlgorithmSteepestDescent >(m, "AlgorithmSteepestDescent")
        .def(py::init<std::string, Model&, Parameters& >())
        .def("Initialize", &AlgorithmSteepestDescent::Initialize)
        ;
 
}

}  // namespace Python.
} // Namespace Kratos

