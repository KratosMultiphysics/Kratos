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
#include "custom_strategies/strategies/explicit_two_step_v_p_strategy.hpp"
#include "custom_strategies/strategies/nodal_two_step_v_p_strategy.h"

//schemes
#include "custom_strategies/schemes/first_order_forward_euler_scheme.hpp"

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
    namespace py = pybind11;

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
      typedef FirstOrderForwardEulerScheme< SparseSpaceType, LocalSpaceType >  FirstOrderForwardEulerSchemeType;

      //custom strategy types
      typedef TwoStepVPStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > TwoStepVPStrategyType;
      typedef ExplicitTwoStepVPStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ExplicitStrategyType;
      typedef NodalTwoStepVPStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > NodalTwoStepVPStrategyType;
      typedef GaussSeidelLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > GaussSeidelLinearStrategyType;
      //********************************************************************
      //*************************SHCHEME CLASSES****************************
      //********************************************************************


      py::class_< TwoStepVPStrategyType, TwoStepVPStrategyType::Pointer, BaseSolvingStrategyType >(m,"TwoStepVPStrategy")
      	.def(py::init<ModelPart&,LinearSolverType::Pointer,LinearSolverType::Pointer,bool,double,double,int,unsigned int,unsigned int>())
        .def("CalculateAccelerations",&TwoStepVPStrategyType::CalculateAccelerations)
        .def("CalculateDisplacements",&TwoStepVPStrategyType::CalculateDisplacementsAndPorosity)
      	// .def("InitializeStressStrain",&TwoStepVPStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::InitializeStressStrain)
      	    ;

      py::class_< ExplicitStrategyType, ExplicitStrategyType::Pointer, BaseSolvingStrategyType >(m,"ExplicitStrategyType")
	.def(py::init < ModelPart&, BaseSchemeType::Pointer,  LinearSolverType::Pointer, bool, bool, bool >())
	.def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer,  bool, bool, bool >())
	.def("SetInitializePerformedFlag", &ExplicitStrategyType::SetInitializePerformedFlag)
	.def("GetInitializePerformedFlag", &ExplicitStrategyType::GetInitializePerformedFlag)
	;

      py::class_< NodalTwoStepVPStrategyType, NodalTwoStepVPStrategyType::Pointer, BaseSolvingStrategyType >(m,"NodalTwoStepVPStrategy")
      	.def(py::init<ModelPart&,LinearSolverType::Pointer,LinearSolverType::Pointer,bool,double,double,int,unsigned int,unsigned int>())
        .def("CalculateAccelerations",&NodalTwoStepVPStrategyType::CalculateAccelerations)
        .def("CalculateDisplacements",&NodalTwoStepVPStrategyType::CalculateDisplacements)
	;

      py::class_< GaussSeidelLinearStrategyType,GaussSeidelLinearStrategyType::Pointer, BaseSolvingStrategyType >(m,"GaussSeidelLinearStrategy")
	      .def(py::init < ModelPart& ,  BaseSchemeType::Pointer, LinearSolverType::Pointer, BuilderAndSolverType::Pointer, bool, bool  >())
        .def("GetResidualNorm", &GaussSeidelLinearStrategyType::GetResidualNorm)
        .def("SetBuilderAndSolver", &GaussSeidelLinearStrategyType::SetBuilderAndSolver)
      ;

      // Explicit scheme: Central differences
      py::class_< FirstOrderForwardEulerSchemeType,FirstOrderForwardEulerSchemeType::Pointer, BaseSchemeType > (m,"FirstOrderForwardEulerSchemeType")
	.def(py::init< const double, const double, const double, const bool >() )
	.def("Initialize", &FirstOrderForwardEulerScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;
    }

  }  // namespace Python.

} // Namespace Kratos

