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
#include "custom_strategies/custom_schemes/iga_contact_scheme.hpp"
// Criterias
#include "custom_strategies/custom_convergence_criteria/active_set_criteria.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos {
namespace Python {

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    // Custom convergence criterion types
    // using ActiveSetCriteriaType = ActiveSetCriteriaCriteria<SparseSpaceType, LocalSpaceType>;

    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    // Base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef ImplicitSolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseImplicitSolvingStrategyType;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
    typedef BuilderAndSolverType::Pointer BuilderAndSolverPointer;

    // Custom strategy types
    typedef EigensolverNitscheStabilizationStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > EigensolverNitscheStabilizationStrategyType;

    // Custom scheme types
    typedef EigensolverNitscheStabilizationScheme< SparseSpaceType, LocalSpaceType > EigensolverNitscheStabilizationSchemeType;
    typedef IgaContactScheme< SparseSpaceType, LocalSpaceType > IgaContactSchemeType;

    // Custom criteria types

    using ConvergenceCriteriaType = ConvergenceCriteria<SparseSpaceType, LocalSpaceType>;
    using BaseSolvingStrategyType = ImplicitSolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    typedef ActiveSetCriteria< SparseSpaceType, LocalSpaceType > ActiveSetCriteriaType;

    // ********************************************************************************
    // STRATEGIES
    // ********************************************************************************

    // Eigensolver Strategy
    py::class_< EigensolverNitscheStabilizationStrategyType, typename EigensolverNitscheStabilizationStrategyType::Pointer,BaseImplicitSolvingStrategyType >(m,"EigensolverNitscheStabilizationStrategy")
        .def(py::init<ModelPart&,
             BaseSchemeType::Pointer,
             BuilderAndSolverPointer>(),
                py::arg("model_part"),
                py::arg("scheme"),
                py::arg("builder_and_solver"))
        ;

    // ********************************************************************************
    //  SCHEMES
    // ********************************************************************************

    // Eigensolver Scheme Type
    py::class_< EigensolverNitscheStabilizationSchemeType,typename EigensolverNitscheStabilizationSchemeType::Pointer, BaseSchemeType>(m,"EigensolverNitscheStabilizationScheme")
        .def(py::init<>() )
        ;

    py::class_< IgaContactSchemeType,typename IgaContactSchemeType::Pointer, BaseSchemeType>(m,"IgaContactScheme")
        .def(py::init<>() )
        ;

    
    // ********************************************************************************
    // CRITERIA
    // ********************************************************************************


    // Custom mortar and criteria
    // py::class_< MortarAndConvergenceCriteriaType, typename MortarAndConvergenceCriteriaType::Pointer,
    //     ConvergenceCriteriaType >
    //     (m, "MortarAndConvergenceCriteria")
    //     .def(py::init<Parameters>())
    //     .def(py::init<ConvergenceCriteriaPointer, ConvergenceCriteriaPointer>())
    //     .def(py::init<ConvergenceCriteriaPointer, ConvergenceCriteriaPointer, bool>())
    //     ;

    // Displacement and lagrange multiplier Convergence Criterion
    py::class_< ActiveSetCriteriaType, typename ActiveSetCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "ActiveSetCriteria")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;
}

}  // namespace Python.
} // Namespace Kratos


