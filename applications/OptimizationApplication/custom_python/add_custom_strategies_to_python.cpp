//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


// System includes


// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "custom_strategies/strategies/helmholtz_strategy.h"

//linear solvers
#include "linear_solvers/linear_solver.h"


namespace Kratos {
namespace Python {

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef ImplicitSolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;


    py::class_<HelmholtzStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType>,
        HelmholtzStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType>::Pointer,
        BaseSolvingStrategyType>(m,"HelmholtzStrategy")
        .def(py::init<ModelPart &, LinearSolverType::Pointer, bool, bool, int, double>());
}

} // namespace Python.
} // Namespace Kratos
