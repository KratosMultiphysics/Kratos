// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

// External includes
#include "spaces/ublas_space.h"

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "includes/kratos_parameters.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/geo_mechanics_newton_raphson_strategy.hpp"
#include "custom_strategies/strategies/geo_mechanics_ramm_arc_length_strategy.hpp"

//builders and solvers

//schemes
#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/newmark_quasistatic_damped_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/newmark_dynamic_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/newmark_quasistatic_Pw_scheme.hpp"
#include "custom_strategies/schemes/backward_euler_quasistatic_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/backward_euler_quasistatic_Pw_scheme.hpp"

//linear solvers
#include "linear_solvers/linear_solver.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

void AddCustomStrategiesToPython(pybind11::module& m)
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector>           LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType >                        LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >  BaseSolvingStrategyType;
    typedef Scheme< SparseSpaceType, LocalSpaceType >                             BaseSchemeType;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType >                ConvergenceCriteriaType;

    typedef NewmarkQuasistaticUPwScheme< SparseSpaceType, LocalSpaceType >        NewmarkQuasistaticUPwSchemeType;
    typedef NewmarkQuasistaticDampedUPwScheme< SparseSpaceType, LocalSpaceType >  NewmarkQuasistaticDampedUPwSchemeType;
    typedef NewmarkDynamicUPwScheme< SparseSpaceType, LocalSpaceType >            NewmarkDynamicUPwSchemeType;
    typedef NewmarkQuasistaticPwScheme< SparseSpaceType, LocalSpaceType >         NewmarkQuasistaticPwSchemeType;
    typedef BackwardEulerQuasistaticUPwScheme< SparseSpaceType, LocalSpaceType >  BackwardEulerQuasistaticUPwSchemeType;
    typedef BackwardEulerQuasistaticPwScheme< SparseSpaceType, LocalSpaceType >   BackwardEulerQuasistaticPwSchemeType;

    typedef GeoMechanicsNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > GeoMechanicsNewtonRaphsonStrategyType;
    typedef GeoMechanicsRammArcLengthStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > GeoMechanicsRammArcLengthStrategyType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    class_< NewmarkQuasistaticUPwSchemeType, typename NewmarkQuasistaticUPwSchemeType::Pointer, BaseSchemeType >
    (m, "NewmarkQuasistaticUPwScheme")
    .def(init<  double, double, double >());

    class_< NewmarkQuasistaticDampedUPwSchemeType, typename NewmarkQuasistaticDampedUPwSchemeType::Pointer, BaseSchemeType >
    (m, "NewmarkQuasistaticDampedUPwScheme")
    .def(init<  double, double, double >());

    class_< NewmarkDynamicUPwSchemeType,typename NewmarkDynamicUPwSchemeType::Pointer, BaseSchemeType >
    (m, "NewmarkDynamicUPwScheme")
    .def(init<  double, double, double >());

    class_< NewmarkQuasistaticPwSchemeType, typename NewmarkQuasistaticPwSchemeType::Pointer, BaseSchemeType >
    (m, "NewmarkQuasistaticPwScheme")
    .def(init<  double >());

    class_< BackwardEulerQuasistaticUPwSchemeType, typename BackwardEulerQuasistaticUPwSchemeType::Pointer, BaseSchemeType >
    (m, "BackwardEulerQuasistaticUPwScheme")
    .def(init< >());

    class_< BackwardEulerQuasistaticPwSchemeType, typename BackwardEulerQuasistaticPwSchemeType::Pointer, BaseSchemeType >
    (m, "BackwardEulerQuasistaticPwScheme")
    .def(init< >());

    class_< GeoMechanicsNewtonRaphsonStrategyType, typename GeoMechanicsNewtonRaphsonStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "GeoMechanicsNewtonRaphsonStrategy")
    .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
        BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >());

    class_< GeoMechanicsRammArcLengthStrategyType, typename GeoMechanicsRammArcLengthStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "GeoMechanicsRammArcLengthStrategy")
    .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
        BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >())
    .def("UpdateLoads",&GeoMechanicsRammArcLengthStrategyType::UpdateLoads);

}

}  // namespace Python.
} // Namespace Kratos
