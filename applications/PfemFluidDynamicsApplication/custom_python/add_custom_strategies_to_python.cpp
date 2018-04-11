//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes 

// External includes 

// Project includes
#include "spaces/ublas_space.h"
#include "utilities/openmp_utils.h"
#include "custom_python/add_custom_strategies_to_python.h"

// Strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/two_step_v_p_strategy.h"
#include "custom_strategies/strategies/gauss_seidel_linear_strategy.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"

//schemes

namespace Kratos
{

namespace Python
{		
using namespace pybind11;

void  AddCustomStrategiesToPython(pybind11::module& m)
{      
  //base types
  typedef UblasSpace<double, CompressedMatrix, Vector>                                   SparseSpaceType;
  typedef UblasSpace<double, Matrix, Vector>                                              LocalSpaceType;
  typedef LinearSolver<SparseSpaceType, LocalSpaceType >                                LinearSolverType;
  typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >   BaseSolvingStrategyType;
  typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType >     BuilderAndSolverType;
  typedef Scheme< SparseSpaceType, LocalSpaceType >                                       BaseSchemeType;

  // Solution strategy types
  typedef TwoStepVPStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >           TwoStepVPStrategyType;
  typedef GaussSeidelLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > GaussSeidelStrategyType;
          
  //*************************STRATEGY CLASSES***************************
      
  class_<TwoStepVPStrategyType, typename TwoStepVPStrategyType::Pointer, BaseSolvingStrategyType>
      (m,"TwoStepVPStrategy")
      .def(init<ModelPart&,LinearSolverType::Pointer,LinearSolverType::Pointer,bool,double,double,int,unsigned int,unsigned int>())
      .def("CalculateAccelerations", &TwoStepVPStrategyType::CalculateAccelerations)
      .def("CalculateDisplacements", &TwoStepVPStrategyType::CalculateDisplacements)
      ;

  class_<GaussSeidelStrategyType, typename GaussSeidelStrategyType::Pointer, BaseSolvingStrategyType>
      (m,"GaussSeidelLinearStrategy")
      .def(init<ModelPart&,BaseSchemeType::Pointer,LinearSolverType::Pointer,BuilderAndSolverType::Pointer,bool,bool>())
      .def("GetResidualNorm", &GaussSeidelStrategyType::GetResidualNorm)
      .def("SetBuilderAndSolver", &GaussSeidelStrategyType::SetBuilderAndSolver)
      ;

}

}  // namespace Python.

} // Namespace Kratos

