//   
//   Project Name:        
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: $
//

// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 

// Project includes
#include "includes/define.h"
#include "containers/flags.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "spaces/ublas_space.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/newton_raphson_strategy.hpp"

//builders and solvers

//schemes
#include "custom_strategies/custom_schemes/newmark_scheme.hpp"


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
    typedef NewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > NewtonRaphsonStrategyType;
    typedef NewmarkScheme< SparseSpaceType, LocalSpaceType >  NewmarkSchemeType;  

    // Newmark Scheme
    class_< NewmarkSchemeType,bases< BaseSchemeType >,  boost::noncopyable >("NewmarkScheme",init< bool >());
      
    // Newton-Raphson Strategy
    class_< NewtonRaphsonStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >("NewtonRaphsonStrategy", 
        init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer, double, double, int, bool, bool, bool >());
}

}  // namespace Python.
} // Namespace Kratos

