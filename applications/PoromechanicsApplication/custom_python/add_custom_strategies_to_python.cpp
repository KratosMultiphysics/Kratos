//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "spaces/ublas_space.h"
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

using namespace boost::python;

void  AddCustomStrategiesToPython()
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
    
    class_< NewmarkQuasistaticUPwSchemeType,bases< BaseSchemeType >, boost::noncopyable >("NewmarkQuasistaticUPwScheme",
        init<  double, double, double >());
    class_< NewmarkQuasistaticDampedUPwSchemeType,bases< BaseSchemeType >, boost::noncopyable >("NewmarkQuasistaticDampedUPwScheme",
        init<  double, double, double, double, double >());
    class_< NewmarkDynamicUPwSchemeType,bases< BaseSchemeType >, boost::noncopyable >("NewmarkDynamicUPwScheme",
        init<  double, double, double, double, double >());

    class_< PoromechanicsNewtonRaphsonStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >("PoromechanicsNewtonRaphsonStrategy", 
        init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
                 BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >());
    class_< PoromechanicsRammArcLengthStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >("PoromechanicsRammArcLengthStrategy", 
        init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
                 BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >())
        .def("UpdateLoads",&PoromechanicsRammArcLengthStrategyType::UpdateLoads);
    class_< PoromechanicsNewtonRaphsonNonlocalStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >("PoromechanicsNewtonRaphsonNonlocalStrategy", 
        init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
                 BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >());
    class_< PoromechanicsRammArcLengthNonlocalStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >("PoromechanicsRammArcLengthNonlocalStrategy", 
        init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
                 BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >());

}

}  // namespace Python.
} // Namespace Kratos
