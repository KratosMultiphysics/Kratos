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
#include "custom_algorithms/algorithm_steepest_descent.h"
#include "custom_algorithms/algorithm_gradient_projection.h"
#include "custom_algorithms/nlopt_optimizer.h"

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include "custom_external_libraries/nlopt/src/api/nlopt.h"


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
    py::class_<AlgorithmSteepestDescent >(m, "AlgorithmSteepestDescent")
        .def(py::init<std::string, Model&, Parameters& >())
        .def("Initialize", &AlgorithmSteepestDescent::Initialize)
        .def("CalculateSolutionStep", &AlgorithmSteepestDescent::CalculateSolutionStep)
        ;  

    py::class_<AlgorithmGradientProjection >(m, "AlgorithmGradientProjection")
        .def(py::init<std::string, Model&, LinearSolver<DenseSpace, DenseSpace>&, Parameters& >())
        .def("Initialize", &AlgorithmGradientProjection::Initialize)
        .def("CalculateSolutionStep", &AlgorithmGradientProjection::CalculateSolutionStep)
        ;  

    // Expose nlopt_algorithm enum
    py::enum_<nlopt_algorithm>(m, "nlopt_algorithm")
        .value("NLOPT_GN_DIRECT", nlopt_algorithm::NLOPT_GN_DIRECT)
        .value("NLOPT_GN_DIRECT_L", nlopt_algorithm::NLOPT_GN_DIRECT_L)
        .value("NLOPT_GN_DIRECT_L_RAND", nlopt_algorithm::NLOPT_GN_DIRECT_L_RAND)
        .value("NLOPT_GN_DIRECT_NOSCAL", nlopt_algorithm::NLOPT_GN_DIRECT_NOSCAL)
        .value("NLOPT_GN_DIRECT_L_NOSCAL", nlopt_algorithm::NLOPT_GN_DIRECT_L_NOSCAL)
        .value("NLOPT_GN_DIRECT_L_RAND_NOSCAL", nlopt_algorithm::NLOPT_GN_DIRECT_L_RAND_NOSCAL)
        .value("NLOPT_GN_ORIG_DIRECT", nlopt_algorithm::NLOPT_GN_ORIG_DIRECT)
        .value("NLOPT_GN_ORIG_DIRECT_L", nlopt_algorithm::NLOPT_GN_ORIG_DIRECT_L)
        .value("NLOPT_GD_STOGO", nlopt_algorithm::NLOPT_GD_STOGO)
        .value("NLOPT_GD_STOGO_RAND", nlopt_algorithm::NLOPT_GD_STOGO_RAND)
        .value("NLOPT_LD_LBFGS_NOCEDAL", nlopt_algorithm::NLOPT_LD_LBFGS_NOCEDAL)
        .value("NLOPT_LD_LBFGS", nlopt_algorithm::NLOPT_LD_LBFGS)
        .value("NLOPT_LN_PRAXIS", nlopt_algorithm::NLOPT_LN_PRAXIS)
        .value("NLOPT_LD_VAR1", nlopt_algorithm::NLOPT_LD_VAR1)
        .value("NLOPT_LD_VAR2", nlopt_algorithm::NLOPT_LD_VAR2)
        .value("NLOPT_LD_TNEWTON", nlopt_algorithm::NLOPT_LD_TNEWTON)
        .value("NLOPT_LD_TNEWTON_RESTART", nlopt_algorithm::NLOPT_LD_TNEWTON_RESTART)
        .value("NLOPT_LD_TNEWTON_PRECOND", nlopt_algorithm::NLOPT_LD_TNEWTON_PRECOND)
        .value("NLOPT_LD_TNEWTON_PRECOND_RESTART", nlopt_algorithm::NLOPT_LD_TNEWTON_PRECOND_RESTART)
        .value("NLOPT_GN_CRS2_LM", nlopt_algorithm::NLOPT_GN_CRS2_LM)
        .value("NLOPT_GN_MLSL", nlopt_algorithm::NLOPT_GN_MLSL)
        .value("NLOPT_GD_MLSL", nlopt_algorithm::NLOPT_GD_MLSL)
        .value("NLOPT_GN_MLSL_LDS", nlopt_algorithm::NLOPT_GN_MLSL_LDS)
        .value("NLOPT_GD_MLSL_LDS", nlopt_algorithm::NLOPT_GD_MLSL_LDS)
        .value("NLOPT_LD_MMA", nlopt_algorithm::NLOPT_LD_MMA)
        .value("NLOPT_LN_COBYLA", nlopt_algorithm::NLOPT_LN_COBYLA)
        .value("NLOPT_LN_NEWUOA", nlopt_algorithm::NLOPT_LN_NEWUOA)
        .value("NLOPT_LN_NEWUOA_BOUND", nlopt_algorithm::NLOPT_LN_NEWUOA_BOUND)
        .value("NLOPT_LN_NELDERMEAD", nlopt_algorithm::NLOPT_LN_NELDERMEAD)
        .value("NLOPT_LN_SBPLX", nlopt_algorithm::NLOPT_LN_SBPLX)
        .value("NLOPT_LN_AUGLAG", nlopt_algorithm::NLOPT_LN_AUGLAG)
        .value("NLOPT_LD_AUGLAG", nlopt_algorithm::NLOPT_LD_AUGLAG)
        .value("NLOPT_LN_AUGLAG_EQ", nlopt_algorithm::NLOPT_LN_AUGLAG_EQ)
        .value("NLOPT_LD_AUGLAG_EQ", nlopt_algorithm::NLOPT_LD_AUGLAG_EQ)
        .value("NLOPT_G_MLSL", nlopt_algorithm::NLOPT_G_MLSL)
        .value("NLOPT_G_MLSL_LDS", nlopt_algorithm::NLOPT_G_MLSL_LDS)
        .value("NLOPT_LD_SLSQP", nlopt_algorithm::NLOPT_LD_SLSQP)
        .value("NLOPT_LD_CCSAQ", nlopt_algorithm::NLOPT_LD_CCSAQ)
        .value("NLOPT_GN_ESCH", nlopt_algorithm::NLOPT_GN_ESCH)
        .value("NLOPT_GN_AGS", nlopt_algorithm::NLOPT_GN_AGS)
        .value("NLOPT_NUM_ALGORITHMS", nlopt_algorithm::NLOPT_NUM_ALGORITHMS);        

    py::class_<NLOptOptimizer>(m, "NLOptOptimizer")
        .def(py::init<nlopt_algorithm, unsigned int>())
        .def("set_objective_function", &NLOptOptimizer::SetObjectiveFunction)
        .def("set_gradient", &NLOptOptimizer::SetGradient)
        .def("set_lower_bounds", &NLOptOptimizer::SetLowerBounds)
        .def("set_upper_bounds", &NLOptOptimizer::SetUpperBounds)
        .def("set_initial_guess", &NLOptOptimizer::SetInitialGuess)
        .def("set_relative_tolerance", &NLOptOptimizer::SetRelativeTolerance)
        .def("set_absolute_tolerance", &NLOptOptimizer::SetAbsoluteTolerance)
        .def("set_max_iterations", &NLOptOptimizer::SetMaxIterations)
        .def("optimize", &NLOptOptimizer::Optimize);        
 
}

}  // namespace Python.
} // Namespace Kratos

