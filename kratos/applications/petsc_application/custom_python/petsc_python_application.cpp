//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: pooyan $
//   Date:                $Date: 2008-07-09 13:14:17 $
//   Revision:            $Revision: 1.1 $
//
//

// System includes 

#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>


// Project includes 
#include "includes/define.h"
#include "petsc_application.h"
#include "petsc_space.h"
#include "spaces/ublas_space.h"
#include "add_petsc_linear_solvers_to_python.h"
#include "custom_strategies/builder_and_solvers/petsc_residualbased_elimination_builder_and_solver.h"
#include "includes/model_part.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

//schemes
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/petsc_residualbased_incrementalupdate_static_scheme.h"

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
//#include "solving_strategies/convergencecriterias/residual_criteria.h"
//#include "solving_strategies/convergencecriterias/and_criteria.h"

//Builder And Solver
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"


//linear solvers
#include "linear_solvers/linear_solver.h"

//utilities
#include "python/pointer_vector_set_python_interface.h"


 
namespace Kratos
{

namespace Python
{

  using namespace boost::python;

  void EraseAll(std::string& ThisString, std::string ToBeRemoved)
  {
    int position;
    while( (position = ThisString.find_first_of(ToBeRemoved)) >= 0 ) 
      {
	ThisString.erase(position ,ToBeRemoved.size() );
      }

  }

  std::string ErrorCleaner(std::string const& Input)
  {
    std::string output(Input);

    EraseAll(output, "boost::numeric::");

    return output;
  }



  typedef PetscSpace<Mat, Vec> PetscSparseSpaceType;
  typedef UblasSpace<double, Matrix, Vector> PetscLocalSpaceType;

	//ADDED BY PAOLO (next two)
	double Dot(PetscSparseSpaceType& dummy, PetscSparseSpaceType::VectorType& rX, PetscSparseSpaceType::VectorType& rY)
	{   return dummy.Dot(rX,rY);} 

	void ScaleAndAdd(PetscSparseSpaceType& dummy, const double A, const PetscSparseSpaceType::VectorType& rX, const double B, PetscSparseSpaceType::VectorType& rY) 
	//(const double A,const  VectorType& rX, const double B, VectorType& rY) // rY = (A * rX) + (B * rY) 
	{	dummy.ScaleAndAdd(A,rX,B, rY);	}

	
	void Mult(PetscSparseSpaceType& dummy, PetscSparseSpaceType::MatrixType& rA, PetscSparseSpaceType::VectorType& rX, PetscSparseSpaceType::VectorType& rY)
	//rY=A*rX (the product is stored inside the rY)
	{   dummy.Mult(rA, rX, rY);}
	
// 	void TransposeMult(PetscSparseSpaceType& dummy, PetscSparseSpaceType::MatrixType& rA, PetscSparseSpaceType::VectorType& rX, PetscSparseSpaceType::VectorType& rY)
// 	//rY=A*rX (the product is stored inside the rY)
// 	{   dummy.TransposeMult(rA, rX, rY);}
	
	PetscSparseSpaceType::IndexType Size(PetscSparseSpaceType& dummy, PetscSparseSpaceType::VectorType const& rV)
	{return dummy.Size(rV);} 

    void ResizeMatrix(PetscSparseSpaceType& dummy, PetscSparseSpaceType::MatrixType& A, unsigned int i1, unsigned int i2)
	{	dummy.Resize(A,i1,i2);	}			

	void ResizeVector(PetscSparseSpaceType& dummy, PetscSparseSpaceType::VectorType& x, unsigned int i1)
	{	dummy.Resize(x,i1);	}

	void SetToZeroMatrix(PetscSparseSpaceType& dummy, PetscSparseSpaceType::MatrixType& A)
	{	dummy.SetToZero(A);	}

	void SetToZeroVector(PetscSparseSpaceType& dummy, PetscSparseSpaceType::VectorType& x)
	{	dummy.SetToZero(x);	}

	void ClearMatrix(PetscSparseSpaceType& dummy, PetscSparseSpaceType::MatrixType& A)
	{	dummy.Clear(A);	}

	void ClearVector(PetscSparseSpaceType& dummy, PetscSparseSpaceType::VectorType& x)
	{	dummy.Clear(x);	}

	double TwoNorm(PetscSparseSpaceType& dummy, PetscSparseSpaceType::VectorType& x)
	{	return dummy.TwoNorm(x);	}

	void UnaliasedAdd(PetscSparseSpaceType& dummy, PetscSparseSpaceType::VectorType& x, const double A, const PetscSparseSpaceType::VectorType& rY) // x+= a*Y
	{	dummy.UnaliasedAdd(x,A,rY);	}

  void MoveMesh( Scheme< PetscSparseSpaceType, PetscLocalSpaceType >& dummy, ModelPart::NodesContainerType& rNodes)
  {	
    for(ModelPart::NodeIterator i = rNodes.begin() ; i != rNodes.end() ; ++i)
      {
	const array_1d<double,3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);
	(i)->X() = (i)->X0() + disp[0];
	(i)->Y() = (i)->Y0() + disp[1];
	(i)->Z() = (i)->Z0() + disp[2];
      }
  }	


  
  BOOST_PYTHON_MODULE(KratosPetscApplication)
  {

	  class_<KratosPetscApplication, 
			  KratosPetscApplication::Pointer, 
			  bases<KratosApplication>, boost::noncopyable >("KratosPetscApplication")
			;

	  typedef LinearSolver<PetscSparseSpaceType, PetscLocalSpaceType > PetscLinearSolverType;
	bool (PetscLinearSolverType::*pointer_to_solve)(PetscLinearSolverType::SparseMatrixType& rA, PetscLinearSolverType::VectorType& rX, PetscLinearSolverType::VectorType& rB) = &PetscLinearSolverType::Solve;
	  class_<PetscLinearSolverType, PetscLinearSolverType::Pointer>("PetscLinearSolver")
	    .def("Initialize",&PetscLinearSolverType::Initialize)
	    .def("Solve",pointer_to_solve)
		  //.def("",&LinearSolverType::)
		  .def(self_ns::str(self))
		  ;



	  typedef SolvingStrategy< PetscSparseSpaceType, PetscLocalSpaceType, PetscLinearSolverType > PetscBaseSolvingStrategyType;
	  typedef Scheme< PetscSparseSpaceType, PetscLocalSpaceType > PetscBaseSchemeType;
	  typedef BuilderAndSolver< PetscSparseSpaceType, PetscLocalSpaceType, PetscLinearSolverType > BuilderAndSolverType;

	  //********************************************************************
	  //********************************************************************
	  //strategy base class
	  class_< PetscBaseSolvingStrategyType, boost::noncopyable >("PetscSolvingStrategy", init< ModelPart&, bool >() )
	    .def("Predict", &PetscBaseSolvingStrategyType::Predict )
	    .def("Solve", &PetscBaseSolvingStrategyType::Solve )
				.def("IsConverged", &PetscBaseSolvingStrategyType::IsConverged )
				.def("CalculateOutputData", &PetscBaseSolvingStrategyType::CalculateOutputData )
				.def("SetEchoLevel", &PetscBaseSolvingStrategyType::SetEchoLevel )
				.def("GetEchoLevel", &PetscBaseSolvingStrategyType::GetEchoLevel )
				.def("SetRebuildLevel", &PetscBaseSolvingStrategyType::SetRebuildLevel )
				.def("GetRebuildLevel", &PetscBaseSolvingStrategyType::GetRebuildLevel )
				.def("SetMoveMeshFlag", &PetscBaseSolvingStrategyType::SetMoveMeshFlag )
				.def("MoveMeshFlag", &PetscBaseSolvingStrategyType::MoveMeshFlag )
				.def("MoveMesh", &PetscBaseSolvingStrategyType::MoveMesh )
				.def("Clear", &PetscBaseSolvingStrategyType::Clear )
				//.def("GetModelPart", &PetscBaseSolvingStrategyType::GetModelPart )
				; 

			
			class_< ResidualBasedLinearStrategy< PetscSparseSpaceType, PetscLocalSpaceType, PetscLinearSolverType >,bases< PetscBaseSolvingStrategyType >,  boost::noncopyable >
				("PetscResidualBasedLinearStrategy", 
				 init<ModelPart&,PetscBaseSchemeType::Pointer, PetscLinearSolverType::Pointer, BuilderAndSolverType::Pointer, bool, bool, bool, bool	>() )
				;

			
 			typedef ConvergenceCriteria< PetscSparseSpaceType, PetscLocalSpaceType > TConvergenceCriteriaType;


 			class_< ResidualBasedNewtonRaphsonStrategy< PetscSparseSpaceType, PetscLocalSpaceType, PetscLinearSolverType >,bases< PetscBaseSolvingStrategyType >,  boost::noncopyable >
 				("PetscResidualBasedNewtonRaphsonStrategy", 
 				init<ModelPart&, PetscBaseSchemeType::Pointer, PetscLinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool
 				>() )
 				;


			//********************************************************************
			//********************************************************************
			class_< PetscBaseSchemeType, boost::noncopyable >
			("PetscScheme", init< >() )
			.def("Initialize", &PetscBaseSchemeType::Initialize )
			.def("SchemeIsInitialized", &PetscBaseSchemeType::SchemeIsInitialized )
			.def("ElementsAreInitialized", &PetscBaseSchemeType::ElementsAreInitialized )
			.def("InitializeElements", &PetscBaseSchemeType::InitializeElements )
			.def("InitializeSolutionStep", &PetscBaseSchemeType::InitializeSolutionStep )
			.def("FinalizeSolutionStep", &PetscBaseSchemeType::FinalizeSolutionStep )
			.def("InitializeNonLinIteration", &PetscBaseSchemeType::InitializeNonLinIteration )
			.def("FinalizeNonLinIteration", &PetscBaseSchemeType::FinalizeNonLinIteration )
			.def("Predict", &PetscBaseSchemeType::Predict )
			.def("Update", &PetscBaseSchemeType::Update )
			.def("CalculateOutputData", &PetscBaseSchemeType::CalculateOutputData )
			.def("Clean", &PetscBaseSchemeType::Clean )
			.def("MoveMesh", MoveMesh )
			;

			class_< PetscResidualBasedIncrementalUpdateStaticScheme< PetscSparseSpaceType, PetscLocalSpaceType>,	
					bases< PetscBaseSchemeType >,  boost::noncopyable >
				(
					"PetscResidualBasedIncrementalUpdateStaticScheme", init< >() 
				);


			//********************************************************************
			//********************************************************************
			//convergence criteria base class
			 class_< ConvergenceCriteria< PetscSparseSpaceType, PetscLocalSpaceType >, boost::noncopyable >("PetscConvergenceCriteria", init<>() )
 				.def("SetActualizeRHSFlag", &ConvergenceCriteria<PetscSparseSpaceType, PetscLocalSpaceType >::SetActualizeRHSFlag )
 			 .def("GetActualizeRHSflag", &ConvergenceCriteria<PetscSparseSpaceType, PetscLocalSpaceType >::GetActualizeRHSflag )
// 			 .def("PreCriteria", &ConvergenceCriteria<PetscSparseSpaceType, PetscLocalSpaceType >::PreCriteria )
// 			 .def("PostCriteria", &ConvergenceCriteria<PetscSparseSpaceType, PetscLocalSpaceType >::PostCriteria )
// 			 .def("Initialize", &ConvergenceCriteria<PetscSparseSpaceType, PetscLocalSpaceType >::Initialize )
// 			 .def("InitializeSolutionStep", &ConvergenceCriteria<PetscSparseSpaceType, PetscLocalSpaceType >::InitializeSolutionStep )
// 			 .def("FinalizeSolutionStep", &ConvergenceCriteria<PetscSparseSpaceType, PetscLocalSpaceType >::FinalizeSolutionStep )
			 ;                  
			
			class_< DisplacementCriteria<PetscSparseSpaceType, PetscLocalSpaceType >,
			         bases<ConvergenceCriteria< PetscSparseSpaceType, PetscLocalSpaceType > >,  
			         boost::noncopyable >
			        ("PetscDisplacementCriteria", init< double, double>() );
			
/*			class_< ResidualCriteria< PetscSparseSpaceType >,
			         bases<ConvergenceCriteria< PetscSparseSpaceType > >,  
			         boost::noncopyable >
			        ("PetscResidualCriteria", init<Model::Pointer, double >() );
			
			class_< AndCriteria< PetscSparseSpaceType >,
			         bases<ConvergenceCriteria< PetscSparseSpaceType > >,  
			         boost::noncopyable >
			        ("PetscAndCriteria", init<Model::Pointer, ConvergenceCriteria< PetscSparseSpaceType >::Pointer, ConvergenceCriteria< PetscSparseSpaceType >::Pointer >()*/
			//);

			//********************************************************************
			//********************************************************************
			//
			
			//Builder and Solver

			class_< BuilderAndSolverType::DofsArrayType, boost::noncopyable >("PetscDofsArrayType",	init<>() );
			 
			class_< BuilderAndSolverType, boost::noncopyable >("PetscBuilderAndSolver",	init<PetscLinearSolverType::Pointer>() )
			 .def("SetCalculateReactionsFlag", &BuilderAndSolverType::SetCalculateReactionsFlag )
			 .def("GetCalculateReactionsFlag", &BuilderAndSolverType::GetCalculateReactionsFlag )
			 .def("SetDofSetIsInitializedFlag", &BuilderAndSolverType::SetDofSetIsInitializedFlag )
			 .def("GetDofSetIsInitializedFlag", &BuilderAndSolverType::GetDofSetIsInitializedFlag )
			 .def("SetReshapeMatrixFlag", &BuilderAndSolverType::SetReshapeMatrixFlag )
			 .def("GetReshapeMatrixFlag", &BuilderAndSolverType::GetReshapeMatrixFlag )
			 .def("GetEquationSystemSize", &BuilderAndSolverType::GetEquationSystemSize )
			 .def("BuildLHS", &BuilderAndSolverType::BuildLHS )
			 .def("BuildRHS", &BuilderAndSolverType::BuildRHS )
			 .def("Build", &BuilderAndSolverType::Build )
			 .def("SystemSolve", &BuilderAndSolverType::SystemSolve )
			 .def("BuildAndSolve", &BuilderAndSolverType::BuildAndSolve )
			 .def("BuildRHSAndSolve", &BuilderAndSolverType::BuildRHSAndSolve )
			 .def("ApplyDirichletConditions", &BuilderAndSolverType::ApplyDirichletConditions )
			 .def("SetUpDofSet", &BuilderAndSolverType::SetUpDofSet )
			 .def("GetDofSet", &BuilderAndSolverType::GetDofSet, return_internal_reference<>() )
			 .def("SetUpSystem", &BuilderAndSolverType::SetUpSystem )
			 .def("ResizeAndInitializeVectors", &BuilderAndSolverType::ResizeAndInitializeVectors )
			 .def("InitializeSolutionStep", &BuilderAndSolverType::InitializeSolutionStep )
			 .def("FinalizeSolutionStep", &BuilderAndSolverType::FinalizeSolutionStep )
			 .def("CalculateReactions", &BuilderAndSolverType::CalculateReactions )
			 .def("Clear", &BuilderAndSolverType::Clear )
			 .def("SetEchoLevel", &BuilderAndSolverType::SetEchoLevel )
			 .def("GetEchoLevel", &BuilderAndSolverType::GetEchoLevel )
			 ;


			//********************************************************************
			//********************************************************************

			class_< PetscSparseSpaceType, boost::noncopyable >("PetscSparseSpace",	init<>() )
			 .def("ClearMatrix", ClearMatrix )
			 .def("ClearVector", ClearVector )
			 .def("ResizeMatrix", ResizeMatrix )
			 .def("ResizeVector", ResizeVector )
			 .def("SetToZeroMatrix", SetToZeroMatrix )
			 .def("SetToZeroVector", SetToZeroVector )
			 .def("TwoNorm", TwoNorm )
			 //the dot product of two vectors
			 .def("Dot", Dot)
			 //the matrix-vector multiplication
			 .def("Mult", Mult)
// 			 .def("TransposeMult", TransposeMult)
			 .def("Size", Size)
			 .def("UnaliasedAdd",UnaliasedAdd)
			 .def("ScaleAndAdd",ScaleAndAdd)
			;
	  typedef PetscResidualBasedEliminationBuilderAndSolver< PetscSparseSpaceType, PetscLocalSpaceType> PetscResidualBasedEliminationBuilderAndSolverType;
	  
	  class_< PetscResidualBasedEliminationBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable> ("PetscResidualBasedEliminationBuilderAndSolver", init<>() );

	AddPetscLinearSolversToPython();

	def("ErrorCleaner", ErrorCleaner);

  }
  
  
}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
