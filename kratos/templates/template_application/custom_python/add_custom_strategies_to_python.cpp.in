//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/implicit_solving_strategy.h"

//linear solvers
#include "linear_solvers/linear_solver.h"


namespace Kratos::Python {

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef ImplicitSolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    //********************************************************************
    //********************************************************************
    //  py::class_< TestStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType>,
    //  		TestStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType>::Pointer,
    //          BaseSolvingStrategyType>(m, "TestStrategy")
    //  	.def(py::init<ModelPart&, LinearSolverType::Pointer, int, int, bool >() )
    //  	.def("MoveNodes",&TestStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::MoveNodes)
    //  	;

}

} // namespace Kratos::Python
