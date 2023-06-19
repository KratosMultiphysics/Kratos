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
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_python/add_custom_optimization_algorithm_to_python.h"
#include "custom_algorithms/algorithm_gradient_projection.h"


// ==============================================================================

namespace Kratos {
namespace Python {



// ==============================================================================
void  AddCustomOptimizationAlgorithmToPython(pybind11::module& m)
{
    namespace py = pybind11;
    typedef UblasSpace<double, Matrix, Vector> DenseSpace;

    // ================================================================
    // 
    // ================================================================
    py::class_<AlgorithmGradientProjection >(m, "AlgorithmGradientProjection")
        .def(py::init<std::string, Model&, LinearSolver<DenseSpace, DenseSpace>&, Parameters& >())
        .def("Initialize", &AlgorithmGradientProjection::Initialize)
        .def("CalculateSolutionStep", &AlgorithmGradientProjection::CalculateSolutionStep)
        ;  
 
}

}  // namespace Python.
} // Namespace Kratos

