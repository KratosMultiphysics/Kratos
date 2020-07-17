//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Raul Bravo
//                  Altug Emiroglu, http://github.com/emiroglu
//
//


// System includes


// External includes

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "spaces/ublas_space.h"

// Strategies
#include "custom_strategies/custom_strategies/modal_derivative_strategy.hpp"

// Schemes
#include "custom_strategies/custom_schemes/modal_derivative_scheme.hpp"

// Builders and solvers
#include "custom_strategies/rom_builder_and_solver.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos {
namespace Python {

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;

    //********************************************************************
    //********************************************************************

    // Custom strategy types
    typedef ModalDerivativeStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ModalDerivativeStrategyType;

    // Custom scheme types
    typedef ModalDerivativeScheme< SparseSpaceType, LocalSpaceType >  ModalDerivativeSchemeType;

    // Custom builder and solvers types
    typedef ROMBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> ROMBuilderAndSolverType;

    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************

    py::class_< ModalDerivativeStrategyType, typename ModalDerivativeStrategyType::Pointer,BaseSolvingStrategyType >(m,"ModalDerivativeStrategy")
        .def(py::init<ModelPart&,
             BaseSchemeType::Pointer,
             BuilderAndSolverType::Pointer,
             Parameters>())
        ;

    //********************************************************************
    //*************************SCHEME CLASSES*****************************
    //********************************************************************

    // Modal Derivative scheme
    py::class_< ModalDerivativeSchemeType,typename ModalDerivativeSchemeType::Pointer, BaseSchemeType>(m,"ModalDerivativeScheme")
        .def(py::init<Variable<double>&,
        Parameters>())
        ;

    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************

    // ROM builder and solver
    py::class_<ROMBuilderAndSolverType, typename ROMBuilderAndSolverType::Pointer, BuilderAndSolverType>(m, "ROMBuilderAndSolver")
        .def(py::init< LinearSolverType::Pointer, Parameters>() )
        ;

}

} // namespace Python.
} // Namespace Kratos

