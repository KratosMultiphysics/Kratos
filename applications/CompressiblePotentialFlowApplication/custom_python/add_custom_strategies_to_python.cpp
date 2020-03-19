//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// External includes
#include "boost/numeric/ublas/vector.hpp"

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "spaces/ublas_space.h"

// Strategies
#include "custom_strategies/upwind_residualbased_newton_raphson_strategy.hpp"
#include "custom_strategies/upwind_line_search_strategy.hpp"

// Linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos {
namespace Python {

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    // Base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef LinearSolverType::Pointer LinearSolverPointer;
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
    typedef ConvergenceCriteriaType::Pointer ConvergenceCriteriaPointer;
    typedef ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonStrategyType;

    // Custom strategy types
    typedef UpwindResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > UpwindResidualBasedNewtonRaphsonStrategyType;
    typedef UpwindLineSearchStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > UpwindLineSearchStrategyType;

    py::class_<UpwindResidualBasedNewtonRaphsonStrategyType,
               typename UpwindResidualBasedNewtonRaphsonStrategyType::Pointer, ResidualBasedNewtonRaphsonStrategyType>(
        m, "UpwindResidualBasedNewtonRaphsonStrategy")
        .def(py::init<ModelPart&, BaseSchemeType::Pointer, LinearSolverPointer,
                      ConvergenceCriteriaPointer, int, bool, bool, bool>());

    py::class_<UpwindLineSearchStrategyType,
               typename UpwindLineSearchStrategyType::Pointer, UpwindResidualBasedNewtonRaphsonStrategyType>(
        m, "UpwindLineSearchStrategy")
        .def(py::init<ModelPart&, BaseSchemeType::Pointer, LinearSolverPointer,
                      ConvergenceCriteriaPointer, int, bool, bool, bool>());
}

}  // namespace Python.
} // Namespace Kratos

