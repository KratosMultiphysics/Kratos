//
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// External includes
#include "spaces/ublas_space.h"

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "includes/kratos_parameters.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/poromechanics_newton_raphson_strategy.hpp"
#include "custom_strategies/strategies/poromechanics_ramm_arc_length_strategy.hpp"
#include "custom_strategies/strategies/poromechanics_newton_raphson_nonlocal_strategy.hpp"
#include "custom_strategies/strategies/poromechanics_ramm_arc_length_nonlocal_strategy.hpp"

//builders and solvers

//schemes
#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/newmark_quasistatic_damped_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/newmark_dynamic_U_Pw_scheme.hpp"

//linear solvers
#include "linear_solvers/linear_solver.h"


namespace Kratos
{

namespace Python
{

namespace py = pybind11;

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;

    typedef NewmarkQuasistaticUPwScheme< SparseSpaceType, LocalSpaceType >  NewmarkQuasistaticUPwSchemeType;
    typedef NewmarkQuasistaticDampedUPwScheme< SparseSpaceType, LocalSpaceType >  NewmarkQuasistaticDampedUPwSchemeType;
    typedef NewmarkDynamicUPwScheme< SparseSpaceType, LocalSpaceType >  NewmarkDynamicUPwSchemeType;

    typedef PoromechanicsNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > PoromechanicsNewtonRaphsonStrategyType;
    typedef PoromechanicsRammArcLengthStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > PoromechanicsRammArcLengthStrategyType;
    typedef PoromechanicsNewtonRaphsonNonlocalStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > PoromechanicsNewtonRaphsonNonlocalStrategyType;
    typedef PoromechanicsRammArcLengthNonlocalStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > PoromechanicsRammArcLengthNonlocalStrategyType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    py::class_< NewmarkQuasistaticUPwSchemeType, typename NewmarkQuasistaticUPwSchemeType::Pointer, BaseSchemeType >
    (m, "NewmarkQuasistaticUPwScheme")
    .def( py::init<  double, double, double >());
    py::class_< NewmarkQuasistaticDampedUPwSchemeType, typename NewmarkQuasistaticDampedUPwSchemeType::Pointer, BaseSchemeType >
    (m, "NewmarkQuasistaticDampedUPwScheme")
    .def( py::init<  double, double, double, double, double >());
    py::class_< NewmarkDynamicUPwSchemeType,typename NewmarkDynamicUPwSchemeType::Pointer, BaseSchemeType >
    (m, "NewmarkDynamicUPwScheme")
    .def( py::init<  double, double, double, double, double >());

    py::class_< PoromechanicsNewtonRaphsonStrategyType, typename PoromechanicsNewtonRaphsonStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "PoromechanicsNewtonRaphsonStrategy")
    .def( py::init< ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
        BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >());

    py::class_< PoromechanicsRammArcLengthStrategyType, typename PoromechanicsRammArcLengthStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "PoromechanicsRammArcLengthStrategy")
    .def( py::init< ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
        BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >())
    .def("UpdateLoads",&PoromechanicsRammArcLengthStrategyType::UpdateLoads);
    py::class_< PoromechanicsNewtonRaphsonNonlocalStrategyType, typename PoromechanicsNewtonRaphsonNonlocalStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "PoromechanicsNewtonRaphsonNonlocalStrategy")
    .def( py::init< ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
        BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >());
    py::class_< PoromechanicsRammArcLengthNonlocalStrategyType, typename PoromechanicsRammArcLengthNonlocalStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "PoromechanicsRammArcLengthNonlocalStrategy")
    .def( py::init< ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
        BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >());

}

}  // namespace Python.
} // Namespace Kratos
