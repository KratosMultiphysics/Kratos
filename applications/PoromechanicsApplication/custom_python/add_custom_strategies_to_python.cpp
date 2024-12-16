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
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "custom_strategies/strategies/poromechanics_newton_raphson_strategy.hpp"
#include "custom_strategies/strategies/poromechanics_ramm_arc_length_strategy.hpp"
#include "custom_strategies/strategies/poromechanics_newton_raphson_nonlocal_strategy.hpp"
#include "custom_strategies/strategies/poromechanics_ramm_arc_length_nonlocal_strategy.hpp"
#include "custom_strategies/strategies/poromechanics_explicit_strategy.hpp"
#include "custom_strategies/strategies/poromechanics_explicit_nonlocal_strategy.hpp"

//builders and solvers

//schemes
#include "custom_strategies/schemes/one-phase_flow/poro_explicit_cd_scheme.hpp"
#include "custom_strategies/schemes/one-phase_flow/poro_explicit_vv_scheme.hpp"

#include "custom_strategies/schemes/one-phase_flow/poro_newmark_quasistatic_U_Pl_scheme.hpp"
#include "custom_strategies/schemes/one-phase_flow/poro_newmark_quasistatic_damped_U_Pl_scheme.hpp"
#include "custom_strategies/schemes/one-phase_flow/poro_newmark_dynamic_U_Pl_scheme.hpp"

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
    typedef ImplicitSolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;

    typedef PoroExplicitCDScheme< SparseSpaceType, LocalSpaceType >  PoroExplicitCDSchemeType;
    typedef PoroExplicitVVScheme< SparseSpaceType, LocalSpaceType >  PoroExplicitVVSchemeType;

    typedef PoroNewmarkQuasistaticUPlScheme< SparseSpaceType, LocalSpaceType >  PoroNewmarkQuasistaticUPlSchemeType;
    typedef PoroNewmarkQuasistaticDampedUPlScheme< SparseSpaceType, LocalSpaceType >  PoroNewmarkQuasistaticDampedUPlSchemeType;
    typedef PoroNewmarkDynamicUPlScheme< SparseSpaceType, LocalSpaceType >  PoroNewmarkDynamicUPlSchemeType;

    typedef PoromechanicsNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > PoromechanicsNewtonRaphsonStrategyType;
    typedef PoromechanicsRammArcLengthStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > PoromechanicsRammArcLengthStrategyType;
    typedef PoromechanicsNewtonRaphsonNonlocalStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > PoromechanicsNewtonRaphsonNonlocalStrategyType;
    typedef PoromechanicsRammArcLengthNonlocalStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > PoromechanicsRammArcLengthNonlocalStrategyType;
    typedef PoromechanicsExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > PoromechanicsExplicitStrategyType;
    typedef PoromechanicsExplicitNonlocalStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > PoromechanicsExplicitNonlocalStrategyType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    py::class_< PoroExplicitCDSchemeType,typename PoroExplicitCDSchemeType::Pointer, BaseSchemeType >
    (m,"PoroExplicitCDScheme")
    .def(py::init< >());
    py::class_< PoroExplicitVVSchemeType,typename PoroExplicitVVSchemeType::Pointer, BaseSchemeType >
    (m,"PoroExplicitVVScheme")
    .def(py::init< >());
    
    py::class_< PoroNewmarkQuasistaticUPlSchemeType, typename PoroNewmarkQuasistaticUPlSchemeType::Pointer, BaseSchemeType >
    (m, "PoroNewmarkQuasistaticUPlScheme")
    .def( py::init<  double, double, double, double >());
    py::class_< PoroNewmarkQuasistaticDampedUPlSchemeType, typename PoroNewmarkQuasistaticDampedUPlSchemeType::Pointer, BaseSchemeType >
    (m, "PoroNewmarkQuasistaticDampedUPlScheme")
    .def( py::init<  double, double, double, double >());
    py::class_< PoroNewmarkDynamicUPlSchemeType,typename PoroNewmarkDynamicUPlSchemeType::Pointer, BaseSchemeType >
    (m, "PoroNewmarkDynamicUPlScheme")
    .def( py::init<  double, double, double, double >());

    py::class_< PoromechanicsNewtonRaphsonStrategyType, typename PoromechanicsNewtonRaphsonStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "PoromechanicsNewtonRaphsonStrategy")
    .def( py::init< ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer,
        BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >());
    py::class_< PoromechanicsRammArcLengthStrategyType, typename PoromechanicsRammArcLengthStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "PoromechanicsRammArcLengthStrategy")
    .def( py::init< ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer,
        BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >())
    .def("UpdateLoads",&PoromechanicsRammArcLengthStrategyType::UpdateLoads);
    py::class_< PoromechanicsNewtonRaphsonNonlocalStrategyType, typename PoromechanicsNewtonRaphsonNonlocalStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "PoromechanicsNewtonRaphsonNonlocalStrategy")
    .def( py::init< ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer,
        BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >());
    py::class_< PoromechanicsRammArcLengthNonlocalStrategyType, typename PoromechanicsRammArcLengthNonlocalStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "PoromechanicsRammArcLengthNonlocalStrategy")
    .def( py::init< ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer,
        BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >());
    py::class_< PoromechanicsExplicitStrategyType, typename PoromechanicsExplicitStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "PoromechanicsExplicitStrategy")
    .def( py::init< ModelPart&, BaseSchemeType::Pointer, Parameters&, bool, bool, bool >());
    py::class_< PoromechanicsExplicitNonlocalStrategyType, typename PoromechanicsExplicitNonlocalStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "PoromechanicsExplicitNonlocalStrategy")
    .def( py::init< ModelPart&, BaseSchemeType::Pointer, Parameters&, bool, bool, bool >());

}

} // namespace Python.
} // Namespace Kratos
