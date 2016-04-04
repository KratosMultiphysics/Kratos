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

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/newton_raphson_strategy.hpp"
#include "custom_strategies/ramm_arc_length_strategy.hpp"

//builders and solvers

//schemes
#include "custom_strategies/custom_schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "custom_strategies/custom_schemes/newmark_dynamic_U_Pw_scheme.hpp"

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

    typedef NewmarkQuasistaticUPwScheme< SparseSpaceType, LocalSpaceType >  NewmarkQuasistaticUPwSchemeType;
    typedef NewmarkDynamicUPwScheme< SparseSpaceType, LocalSpaceType >  NewmarkDynamicUPwSchemeType; 
    
    typedef NewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > NewtonRaphsonStrategyType;
    typedef RammArcLengthStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > RammArcLengthStrategyType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    class_< NewmarkQuasistaticUPwSchemeType,bases< BaseSchemeType >,  boost::noncopyable >("NewmarkQuasistaticUPwScheme",
        init<  ModelPart&, double, double, double >());

    class_< NewmarkDynamicUPwSchemeType,bases< BaseSchemeType >,  boost::noncopyable >("NewmarkDynamicUPwScheme",
        init<  ModelPart&, double, double, double, double, double >());
        

    class_< NewtonRaphsonStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >("NewtonRaphsonStrategy", 
        init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, double, double, int, bool, bool, bool >());

    class_< RammArcLengthStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >("RammArcLengthStrategy", 
        init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, double, double, int, int, double, double, bool, bool, bool >());
}

}  // namespace Python.
} // Namespace Kratos

