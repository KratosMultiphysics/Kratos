// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp>

// Project includes
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "includes/model_part.h"

// Application includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_strategies/custom_strategies/adjoint_fluid_strategy.h"

namespace Kratos
{

namespace Python
{
  using namespace boost::python;

  void AddCustomStrategiesToPython()
  {
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
    typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;

    class_< AdjointFluidStrategy<SparseSpaceType,
    				  LocalSpaceType,
    				  LinearSolverType>,
    	    bases<SolvingStrategyType>,
    	    boost::noncopyable >
        ( "AdjointFluidStrategy",
          init<ModelPart&,
          LinearSolverType::Pointer,int>() )
        .def( "SetDragForceDirection",
    	  &AdjointFluidStrategy<SparseSpaceType,
    				 LocalSpaceType,
    				 LinearSolverType>::SetDragForceDirection )
        .def( "ComputeSensitivity",
    	  &AdjointFluidStrategy<SparseSpaceType,
    				 LocalSpaceType,
    				 LinearSolverType>::ComputeSensitivity )
          ;

  }

} // namespace Python

} // namespace Kratos
