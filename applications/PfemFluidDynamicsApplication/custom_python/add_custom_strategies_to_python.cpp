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
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"

#include "custom_strategies/strategies/two_step_v_p_strategy.h"
#include "custom_strategies/strategies/gauss_seidel_linear_strategy.h"

// builder_and_solvers

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"


//linear solvers
#include "linear_solvers/linear_solver.h"

//schemes

namespace Kratos
{

  namespace Python
  {
    using namespace pybind11;

    void  AddCustomStrategiesToPython(pybind11::module& m)
    {
      typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
      typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

      //base types
      typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
      typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
      typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
      typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
      //typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaBaseType;


      //********************************************************************
      //*************************SHCHEME CLASSES****************************
      //********************************************************************

      class_< TwoStepVPStrategy< SparseSpaceType,LocalSpaceType, LinearSolverType > >
	      (m, "TwoStepVPStrategy")
	      .def(init<ModelPart&,LinearSolverType::Pointer,LinearSolverType::Pointer,bool,double,double,int,unsigned int,unsigned int>())
        .def("CalculateAccelerations",&TwoStepVPStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::CalculateAccelerations)
        .def("CalculateDisplacements",&TwoStepVPStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::CalculateDisplacements)
      	// .def("InitializeStressStrain",&TwoStepVPStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::InitializeStressStrain)
	    ;

      class_< GaussSeidelLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > >
	      (m,"GaussSeidelLinearStrategy")
	      .def(init < ModelPart& ,  BaseSchemeType::Pointer, LinearSolverType::Pointer, BuilderAndSolverType::Pointer, bool, bool  >())
        .def("GetResidualNorm", &GaussSeidelLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetResidualNorm)
        .def("SetBuilderAndSolver", &GaussSeidelLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetBuilderAndSolver)
      ;

    }

  }  // namespace Python.

} // Namespace Kratos

