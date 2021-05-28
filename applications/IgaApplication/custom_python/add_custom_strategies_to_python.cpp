// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes
#include "boost/numeric/ublas/vector.hpp"

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

// Strategies
#include "custom_strategies/custom_strategies/eigensolver_nitsche_stabilization_strategy.hpp"
// Schemes
#include "custom_strategies/custom_schemes/eigensolver_nitsche_stabilization_scheme.hpp"

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
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
    typedef BuilderAndSolverType::Pointer BuilderAndSolverPointer;
    
    // Custom strategy types
    typedef EigensolverNitscheStabilizationStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > EigensolverNitscheStabilizationStrategyType;

    // Custom scheme types
    typedef EigensolverNitscheStabilizationScheme< SparseSpaceType, LocalSpaceType > EigensolverNitscheStabilizationSchemeType;

    // Eigensolver Strategy
    py::class_< EigensolverNitscheStabilizationStrategyType, typename EigensolverNitscheStabilizationStrategyType::Pointer,BaseSolvingStrategyType >(m,"EigensolverNitscheStabilizationStrategy")
        .def(py::init<ModelPart&,
             BaseSchemeType::Pointer,
             BuilderAndSolverPointer>(),
                py::arg("model_part"),
                py::arg("scheme"),
                py::arg("builder_and_solver"))
        ;

    // Eigensolver Scheme Type
    py::class_< EigensolverNitscheStabilizationSchemeType,typename EigensolverNitscheStabilizationSchemeType::Pointer, BaseSchemeType>(m,"EigensolverNitscheStabilizationScheme")
        .def(py::init<>() )
        ;
}

}  // namespace Python.
} // Namespace Kratos


