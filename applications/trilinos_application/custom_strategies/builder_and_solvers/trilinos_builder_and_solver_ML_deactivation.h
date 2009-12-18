/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
/* *********************************************************   
*          
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2009-01-13 15:39:56 $
*   Revision:            $Revision: 1.5 $
*
* ***********************************************************/


#if !defined(KRATOS_TRILINOS_BUILDER_AND_SOLVER_ML_DEACTIVATION )
#define  KRATOS_TRILINOS_BUILDER_AND_SOLVER_ML_DEACTIVATION


/* System includes */
#include <set>
/* #include <omp.h> */

/* External includes */
#include "boost/smart_ptr.hpp"
#include "boost/timer.hpp" 


/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "Epetra_MpiComm.h"

//trilinos includes
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
// #include "epetra_test_err.h"



//aztec solver includes
#include "AztecOO.h"

#include "Amesos.h"
// #include "AmesosClassType.h"
#include "Epetra_LinearProblem.h"
#include "ml_MultiLevelPreconditioner.h"


namespace Kratos
{

	/**@name Kratos Globals */
	/*@{ */


	/*@} */
	/**@name Type Definitions */       
	/*@{ */

	/*@} */


	/**@name  Enum's */       
	/*@{ */


	/*@} */
	/**@name  Functions */       
	/*@{ */



	/*@} */
	/**@name Kratos Classes */
	/*@{ */

	/** Short class definition.

	Detail class definition.

	Current class provides an implementation for standard builder and solving operations.

	the RHS is constituted by the unbalanced loads (residual)

	Degrees of freedom are reordered putting the restrained degrees of freedom at 
	the end of the system ordered in reverse order with respect to the DofSet.

	Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
	this information.

	Calculation of the reactions involves a cost very similiar to the calculation of the total residual

	\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

	\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

	\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

	\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


	\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

	\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

	\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

	\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


	*/
	template<class TSparseSpace, 
	class TDenseSpace,
	class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
	>
	class TrilinosBuilderAndSolverMLDeactivation 
		: public BuilderAndSolver< TSparseSpace,TDenseSpace, LinearSolver<TSparseSpace,TDenseSpace> >
	{
	public:
		/**@name Type Definitions */       
		/*@{ */
		//typedef boost::shared_ptr< TrilinosBuilderAndSolverMLDeactivation<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;		
		KRATOS_CLASS_POINTER_DEFINITION( TrilinosBuilderAndSolverMLDeactivation );


		typedef BuilderAndSolver<TSparseSpace,TDenseSpace, LinearSolver<TSparseSpace,TDenseSpace> > BaseType;

		typedef TSparseSpace SparseSpaceType;

		typedef typename BaseType::TSchemeType TSchemeType;

		typedef typename BaseType::TDataType TDataType;

		typedef typename BaseType::DofsArrayType DofsArrayType;

		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

		typedef typename BaseType::TSystemVectorType TSystemVectorType;

		typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

		typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
		typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
		typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
 
		typedef typename BaseType::NodesArrayType NodesArrayType;
		typedef typename BaseType::ElementsArrayType ElementsArrayType;
		typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

		typedef typename BaseType::ElementsContainerType ElementsContainerType;

		/*@} */
		/**@name Life Cycle 
		*/    
		/*@{ */

		/** Constructor.
		*/
		TrilinosBuilderAndSolverMLDeactivation(
			Epetra_MpiComm& Comm,
			int guess_row_size,
			typename TLinearSolver::Pointer pNewLinearSystemSolver)
			: BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >(pNewLinearSystemSolver)
			, mrComm(Comm),mguess_row_size(guess_row_size)
		{
			

			/* 			std::cout << "using the standard builder and solver " << std::endl; */

		}
// 		TrilinosBuilderAndSolverMLDeactivation(
// 			)
// 		  : BaseType(typename LinearSolver<TSparseSpace,TDenseSpace>::Pointer(new LinearSolver<TSparseSpace,TDenseSpace>))
// 		{
// 
// 			/* 			std::cout << "using the standard builder and solver " << std::endl; */
// 
// 		}


		/** Destructor.
		*/
		virtual ~TrilinosBuilderAndSolverMLDeactivation(){}


		/*@} */
		/**@name Operators 
		*/  
		/*@{ */

		//**************************************************************************
		//**************************************************************************
		void Build(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& b)
		{
			KRATOS_TRY
				if(!pScheme)
					KRATOS_ERROR(std::runtime_error, "No scheme provided!", "");

			//getting the elements from the model
			ElementsArrayType& pElements = r_model_part.Elements();

			//getting the array of the conditions
			ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

			//resetting to zero the vector of reactions
 			TSparseSpace::SetToZero(*BaseType::mpReactionsVector);

			//contributions to the system
			LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
			LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

			//vector containing the localization in the system of the different
			//terms
			Element::EquationIdVectorType EquationId;

			//			int rank = A.Comm().MyPID(); //getting the processor Id

			ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
			// assemble all elements
			for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
			{
                                if( ! (*it)->GetValue( IS_INACTIVE ) )
                                {
				    //calculate elemental contribution
                                        pScheme->CalculateSystemContributions(*it,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);

				    //assemble the elemental contribution
				    TSparseSpace::AssembleLHS(A,LHS_Contribution,EquationId);
    				    TSparseSpace::AssembleRHS(b,RHS_Contribution,EquationId);

				    // clean local elemental memory
				    pScheme->CleanMemory(*it);
                                }
			}

			LHS_Contribution.resize(0,0,false);
			RHS_Contribution.resize(0,false);

			// assemble all conditions
			for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
			{
                                if( ! (*it)->GetValue( IS_INACTIVE ) )
                                {
				    //calculate elemental contribution
				    pScheme->Condition_CalculateSystemContributions(*it,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);

				    //assemble the elemental contribution
				    TSparseSpace::AssembleLHS(A,LHS_Contribution,EquationId);
				    TSparseSpace::AssembleRHS(b,RHS_Contribution,EquationId);
                                }
			}

			//finalizing the assembly
			A.GlobalAssemble();
			b.GlobalAssemble();

			KRATOS_CATCH("")

		}

		//**************************************************************************
		//**************************************************************************
		void BuildLHS(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemMatrixType& A)
		{
			KRATOS_TRY

			//getting the elements from the model
			ElementsArrayType& pElements = r_model_part.Elements();

			//getting the array of the conditions
			ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

			//resetting to zero the vector of reactions
// 			TSparseSpace::SetToZero(BaseType::mReactionsVector);

			//contributions to the system
			LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);

			//vector containing the localization in the system of the different
			//terms
			Element::EquationIdVectorType EquationId;

			ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

			// assemble all elements
			for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
			{
                                if( ! (*it)->GetValue( IS_INACTIVE ) )
                                {

				//calculate elemental contribution
				pScheme->Calculate_LHS_Contribution(*it,LHS_Contribution,EquationId,CurrentProcessInfo);

				//assemble the elemental contribution
				TSparseSpace::AssembleLHS(A,LHS_Contribution,EquationId);

				// clean local elemental memory
				pScheme->CleanMemory(*it);
                                }
			}

			LHS_Contribution.resize(0,0,false);

			// assemble all conditions
			for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
			{
                                if( ! (*it)->GetValue( IS_INACTIVE ) )
                                {

				//calculate elemental contribution
				pScheme->Condition_Calculate_LHS_Contribution(*it,LHS_Contribution,EquationId,CurrentProcessInfo);

				//assemble the elemental contribution
				TSparseSpace::AssembleLHS(A,LHS_Contribution,EquationId);
                                }
			}

			//finalizing the assembly
			A.GlobalAssemble();
			KRATOS_CATCH("")

		}
		//**************************************************************************
		//**************************************************************************
		void SystemSolve(
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b
			)
		{
			KRATOS_TRY


			double norm_b;
			if(TSparseSpace::Size(b) != 0)
				norm_b = TSparseSpace::TwoNorm(b);
			else
				norm_b = 0.00;

			if(norm_b != 0.00)
			  {
				Epetra_LinearProblem AztecProblem(&A,&Dx,&b);

				Teuchos::ParameterList MLList;
				ML_Epetra::SetDefaults("SA",MLList);
				MLList.set("print unused",1);
				MLList.set("output",10);
				MLList.set("increasing or decreasing","increasing");				
				MLList.set("max levels",3);

				MLList.set("aggregation: threshold",0.0);
				MLList.set("aggregation: damping factor",1.33);
				MLList.set("aggregation: type","Uncoupled");
				
				MLList.set("smoother: pre or post","both");
				//MLList.set("smoother: type (level 0)","MLS");
				//MLList.set("smoother: MLS polynomial order (level 0)",3);

				//MLList.set("smoother: type (level 0)","Gauss-Seidel");
				MLList.set("smoother: type (level 0)","symmetric Gauss-Seidel");
				MLList.set("smoother: sweeps (level 0)",2);
				MLList.set("smoother: damping factor (level 0)",0.89);

				MLList.set("coarse: max size",8000);
				MLList.set("coarse: type","Amesos-Superludist");
				//MLList.set("coarse: type","Amesos-KLU");

				// computation of the nullspace
				int numdf = 2; 
				int dimns = 3; //
//				int lrows = A.NumMyCols (); //number of rows for calling processor
				MLList.set("PDE equations",numdf);
				MLList.set("null space: dimension",dimns);	
				//MLList.set("null space: type","pre-computed");
				//MLList.set("null space: add default vectors",false);
				
				// create the preconditioner
				ML_Epetra::MultiLevelPreconditioner* MLPrec = new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);
		
				// create an AztecOO solver
				AztecOO Solver(AztecProblem); 
				// set preconditioner and solve
				Solver.SetPrecOperator(MLPrec);
				Solver.SetAztecOption(AZ_solver, AZ_gmres);
				Solver.SetAztecOption(AZ_kspace, 200);
				//Solver.SetAztecOption(AZ_output,32);
				
				int mmax_iter=300;
				Solver.Iterate(mmax_iter, 1e-12);
				delete MLPrec;

//				BaseType::mpLinearSystemSolver->Solve(A,Dx,b);
	
			  }
			else
			{
				TSparseSpace::SetToZero(Dx);
			}

			//prints informations about the current time
			if (this->GetEchoLevel()>1)
			{
				std::cout << *(BaseType::mpLinearSystemSolver) << std::endl;
			}

			KRATOS_CATCH("")

		}


		//**************************************************************************
		//**************************************************************************
		void SystemSolveML(
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b,
			ModelPart& r_model_part
			)
		{
			KRATOS_TRY

			double norm_b;
			if(TSparseSpace::Size(b) != 0)
				norm_b = TSparseSpace::TwoNorm(b);
			else
				norm_b = 0.00;

			if(norm_b != 0.00)
			  {
				KRATOS_WATCH("entering in the solver");

				int rank;
				MPI_Comm_rank(MPI_COMM_WORLD,&rank);				

				Epetra_LinearProblem AztecProblem(&A,&Dx,&b);

				Teuchos::ParameterList MLList;
				ML_Epetra::SetDefaults("SA",MLList);
				MLList.set("print unused",1);
				MLList.set("output",10);
				MLList.set("increasing or decreasing","increasing");				
				MLList.set("max levels",4);

				MLList.set("aggregation: threshold",0.0);
				MLList.set("aggregation: damping factor",1.33);
				MLList.set("aggregation: type","Uncoupled");
				
				MLList.set("smoother: pre or post","both");
				//MLList.set("smoother: type (level 0)","MLS");
				//MLList.set("smoother: MLS polynomial order (level 0)",3);

				//MLList.set("smoother: type (level 0)","Aztec");
                //MLList.set("smoother: type (level 0)","Gauss-Seidel");
				MLList.set("smoother: type (level 0)","symmetric Gauss-Seidel");
				MLList.set("smoother: sweeps (level 0)",2);
				MLList.set("smoother: damping factor (level 0)",0.89);

				MLList.set("coarse: max size",8000);
				//MLList.set("coarse: type","Amesos-Superludist");
				MLList.set("coarse: type","Amesos-KLU");

				// computation of the nullspace
				int numdf; // dofs per node
				int dimns; // dimension of the null space
//				int lrows =  A.NumMyRows(); //number of rows for calling processor

				//Teuchos::RCP<vector<double> >  ns;
				boost::shared_ptr<vector<double> > ns;
				double* nullsp;				

				GenerateNullSpace(A, r_model_part, nullsp, ns, numdf, dimns);							

				nullsp= &((*ns)[0]);				
				
				MLList.set("PDE equations",numdf);
				MLList.set("null space: dimension",dimns);	
				MLList.set("null space: type","pre-computed");
				MLList.set("null space: add default vectors",false);
				MLList.set("null space: vectors",nullsp);
				
				// create the preconditioner
				ML_Epetra::MultiLevelPreconditioner* MLPrec = new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);
		
				// create an AztecOO solver
				AztecOO Solver(AztecProblem); 

				// set preconditioner and solve
				Solver.SetPrecOperator(MLPrec);
				Solver.SetAztecOption(AZ_solver, AZ_gmres);
				Solver.SetAztecOption(AZ_kspace, 5000);
				Solver.SetAztecOption(AZ_output,50); //SetAztecOption(AZ_output, AZ_none);
								
				int mmax_iter=5000;
				Solver.Iterate(mmax_iter, 1e-9);
				delete MLPrec;

		//		BaseType::mpLinearSystemSolver->Solve(A,Dx,b);
	
			  }
			else
			{
				TSparseSpace::SetToZero(Dx);
			}

			//prints informations about the current time
			if (this->GetEchoLevel()>1)
			{
				std::cout << *(BaseType::mpLinearSystemSolver) << std::endl;
			}

			KRATOS_CATCH("")

		}

		//**************************************************************************
		//**************************************************************************


		void BuildAndSolve(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
			KRATOS_TRY

			boost::timer building_time;

			Build(pScheme,r_model_part,A,b);

			if(BaseType::GetEchoLevel()>0)
			{
				std::cout << "Building Time : " << building_time.elapsed() << std::endl;
			}
			
			//apply dirichlet conditions
			ApplyDirichletConditions(pScheme,r_model_part,A,Dx,b);

			if (BaseType::GetEchoLevel()== 3)
			{
				std::cout << "before the solution of the system" << std::endl;
				std::cout << "System Matrix = " << A << std::endl;
				std::cout << "unknowns vector = " << Dx << std::endl;
				std::cout << "RHS vector = " << b << std::endl;
			}

			boost::timer solve_time;

			SystemSolveML(A,Dx,b,r_model_part);

			if(BaseType::GetEchoLevel()>0)
			{
				std::cout << "System Solve Time : " << solve_time.elapsed() << std::endl;
			}
			if (BaseType::GetEchoLevel()== 3)
			{
				std::cout << "after the solution of the system" << std::endl;
				std::cout << "System Matrix = " << A << std::endl;
				std::cout << "unknowns vector = " << Dx << std::endl;
				std::cout << "RHS vector = " << b << std::endl;
			}
			
			KRATOS_CATCH("")
		}

		//**************************************************************************
		//**************************************************************************
		void BuildRHSAndSolve(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
			KRATOS_TRY

			BuildRHS(pScheme,r_model_part,b);
			SystemSolve(A,Dx,b);

			KRATOS_CATCH("")
		}

		//**************************************************************************
		//**************************************************************************
		void BuildRHS(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemVectorType& b)
		{
			KRATOS_TRY		

			  //Getting the Elements
			  ElementsArrayType& pElements = r_model_part.Elements();

			  //getting the array of the conditions
			  ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

			  ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
			  
			  //resetting to zero the vector of reactions
			  // 			TSparseSpace::SetToZero(BaseType::mReactionsVector);
			  
			  //contributions to the system
			  LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
			  
			  //vector containing the localization in the system of the different
			  //terms
			  Element::EquationIdVectorType EquationId;
			  
			  // assemble all elements
			  for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
			    {
                                if( ! (*it)->GetValue( IS_INACTIVE ) )
                                {

			      //calculate elemental Right Hand Side Contribution
			      pScheme->Calculate_RHS_Contribution(*it,RHS_Contribution,EquationId,CurrentProcessInfo);
			      
			      //assemble the elemental contribution
			     TSparseSpace::AssembleRHS(b,RHS_Contribution,EquationId);
                                }
			    }
			  
			  RHS_Contribution.resize(0,false);
			  
			  // assemble all conditions
			  for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
			    {
                                if( ! (*it)->GetValue( IS_INACTIVE ) )
                                {

			      //calculate elemental contribution
			      pScheme->Condition_Calculate_RHS_Contribution(*it,RHS_Contribution,EquationId,CurrentProcessInfo);
			      
			      //assemble the elemental contribution
			      TSparseSpace::AssembleRHS(b,RHS_Contribution,EquationId);
                                }
			    }
			  
			//finalizing the assembly
			b.GlobalAssemble();

			  KRATOS_CATCH("")
			    
		}				
		//**************************************************************************
		//**************************************************************************
		void SetUpDofSet(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part
			) 
		{
			KRATOS_TRY

			  //Gets the array of elements from the modeler
 			  ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();
/*  			  ElementsArrayType& pElements = r_model_part.Elements(ModelPart::Kratos_Local); */

			  Element::DofsVectorType ElementalDofList;

			  ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

			  DofsArrayType Doftemp;
			  BaseType::mDofSet = DofsArrayType();

			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);

			  for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
			    {
                                if( ! (*it)->GetValue( IS_INACTIVE ) )
                                {

			      // gets list of Dof involved on every element
				pScheme->GetElementalDofList(*it,ElementalDofList,CurrentProcessInfo);
				
				for(typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ; i != ElementalDofList.end() ; ++i)
				{
					Doftemp.push_back(*i);
				}
                                }
			}
	
			//taking in account conditions
			ConditionsArrayType& pConditions = r_model_part.Conditions();
			for (typename ConditionsArrayType::ptr_iterator it=pConditions.ptr_begin(); it!=pConditions.ptr_end(); ++it)
			{
                                if( ! (*it)->GetValue( IS_INACTIVE ) )
                                {

				// gets list of Dof involved on every element
				pScheme->GetConditionDofList(*it,ElementalDofList,CurrentProcessInfo);

				for(typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ; i != ElementalDofList.end() ; ++i)
				{
					//mDofSet.push_back(*i);
					Doftemp.push_back(*i);
				}
                                }
			}


			Doftemp.Unique();

			BaseType::mDofSet = Doftemp;

			//throws an execption if there are no Degrees of freedom involved in the analysis
			if (BaseType::mDofSet.size()==0)
				KRATOS_ERROR(std::logic_error, "No degrees of freedom!", "");

			BaseType::mDofSetIsInitialized = true;

			KRATOS_CATCH("")
		}

		//**************************************************************************
		//**************************************************************************
		void SetUpSystem(
			ModelPart& r_model_part
			)
		{ 

			// Set equation id for degrees of freedom

			int free_size = 0;
			//int fixed_size = 0;

			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);

			// Calculating number of fixed and free dofs
			for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
			  if(dof_iterator->GetSolutionStepValue(PARTITION_INDEX) == rank)
			    {
				  free_size++;
			    }

			// Calculating the total size and required offset
			//int fixed_offset;
			int free_offset;
			int global_size;

			// The correspounding offset by the sum of the sizes in thread with inferior rank
			MPI_Scan(&free_size, &free_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

			// The total size by the sum of all size in all threads
			MPI_Allreduce(&free_size, &global_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

			// finding the offset for the begining of the partition
			free_offset -= free_size;

			// Now setting the equation id with .
			for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
			  if(dof_iterator->GetSolutionStepValue(PARTITION_INDEX) == rank)
			    {
				dof_iterator->SetEquationId(free_offset++);
//  				std::cout << rank << " : set eq. id for dof " << dof_iterator->Id() << " to " << dof_iterator->EquationId() << std::endl;
			    }


			BaseType::mEquationSystemSize = global_size;
			mLocalSystemSize = free_size;
			std::cout << rank << " : BaseType::mEquationSystemSize = " << BaseType::mEquationSystemSize << std::endl; 
			std::cout << rank << " : mLocalSystemSize = " << mLocalSystemSize << std::endl;
			std::cout << rank << " : free_offset = " << free_offset << std::endl; 
			//std::cout << rank << " : fixed_offset = " << fixed_offset << std::endl;

			//by Riccardo ... it may be wrong!	
			mFirstMyId = free_offset-mLocalSystemSize;
			mLastMyId = mFirstMyId+mLocalSystemSize;

			r_model_part.GetCommunicator().SynchronizeDofs();	

		}
 
	  void UpdateGhostDofs(ModelPart& rThisModelPart)
		{
		  int rank;
		  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

// 		  std::cout << rank << " : Strarting UpdateGhostDofs...." << std::endl;

		  //int source=rank;
		  int destination=0;

// 		  vector<int>& neighbours_indices = rThisModelPart[NEIGHBOURS_INDICES];
		  vector<int>& neighbours_indices = rThisModelPart.GetCommunicator().NeighbourIndices();

// 		  std::cout << rank << " starting domain loop " << std::endl;
		  for(unsigned int i_domain = 0 ; i_domain <  neighbours_indices.size() ; i_domain++)
		    if((destination = neighbours_indices[i_domain]) >= 0)
		      {
// 			std::cout << rank << " domian #" << i_domain << std::endl;
			unsigned int send_buffer_size = 0;
			unsigned int receive_buffer_size = 0;

// 			std::cout << rank;
// 			KRATOS_WATCH(destination);
			// Calculating send and received buffer size
			// The interface meshes are stored after all, local and ghost meshes
			NodesArrayType& r_interface_nodes = rThisModelPart.GetCommunicator().LocalMesh(i_domain).Nodes();
			NodesArrayType& r_ghost_nodes = rThisModelPart.GetCommunicator().GhostMesh().Nodes();

// 			std::cout << rank << " : 2...." << std::endl;
			for(typename NodesArrayType::iterator i_node = r_interface_nodes.begin(); i_node != r_interface_nodes.end(); ++i_node)
			  send_buffer_size += i_node->GetDofs().size();

// 			std::cout << rank << " : 3...." << std::endl;
			for(typename NodesArrayType::iterator i_node = r_ghost_nodes.begin(); i_node != r_ghost_nodes.end(); ++i_node)
			  if(i_node->GetSolutionStepValue(PARTITION_INDEX) == destination)
			    {
			    receive_buffer_size += i_node->GetDofs().size();

			    }
			unsigned int position = 0;
			int* send_buffer = new int[send_buffer_size];
			int* receive_buffer = new int[receive_buffer_size];


			// Filling the buffer
			std::cout << rank << " :  Filling the buffer...." << std::endl;
			for(ModelPart::NodeIterator i_node = r_interface_nodes.begin(); i_node != r_interface_nodes.end(); ++i_node)
 			  for(ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++)
			    {
			      send_buffer[position++] = i_dof->EquationId();

			    }


			MPI_Status status;


			if(position > send_buffer_size)
			  std::cout << rank << " Error in estimating send buffer size...." << std::endl;
		  

			int send_tag = 1;//i_domain;
			int receive_tag = 1;//i_domain;


			MPI_Sendrecv (send_buffer, send_buffer_size, MPI_INT, destination, send_tag, receive_buffer, receive_buffer_size, MPI_INT, destination, receive_tag,
 				 MPI_COMM_WORLD, &status);

// 			std::cout << rank << " : Send and receive Finished" << std::endl;
			
			// Updating nodes
			position = 0;
			for(ModelPart::NodeIterator i_node = rThisModelPart.GetCommunicator().GhostMesh().NodesBegin() ; 
			    i_node != rThisModelPart.GetCommunicator().GhostMesh().NodesEnd() ; i_node++)
// 			for(ModelPart::NodeIterator i_node = rThisModelPart.NodesBegin(ModelPart::Kratos_Ghost) ; 
// 			    i_node != rThisModelPart.NodesEnd(ModelPart::Kratos_Ghost) ; i_node++)
			  if(i_node->GetSolutionStepValue(PARTITION_INDEX) == destination)
			    for(ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++)
			  {
			    i_dof->SetEquationId(receive_buffer[position++]);
			    
			  }

			
			if(position > receive_buffer_size)
			  std::cout << rank << " Error in estimating receive buffer size...." << std::endl;
			
			delete [] send_buffer;
			delete [] receive_buffer;
		      }
		  
		}

		//**************************************************************************
		//**************************************************************************
		void ResizeAndInitializeVectors(
			TSystemMatrixPointerType& pA,
			TSystemVectorPointerType& pDx,
			TSystemVectorPointerType& pb,
			ElementsArrayType& rElements,
			ConditionsArrayType& rConditions,
			ProcessInfo& CurrentProcessInfo
			)
		{
			KRATOS_TRY

				std::cout << "entering ResizeAndInitializeVectors" << std::endl;

				//resizing the system vectors and matrix
				if ( pA == NULL || TSparseSpace::Size1(*pA) == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
				{
					//creating a work array
					unsigned int number_of_local_dofs = mLastMyId - mFirstMyId;

					int temp_size = number_of_local_dofs;
					if(temp_size <1000) temp_size = 1000;
					int* temp = new int[temp_size]; //

				
					//generate map - use the "temp" array here 
					for(unsigned int i=0; i!=number_of_local_dofs; i++)
						temp[i] = mFirstMyId+i;
					Epetra_Map my_map(-1, number_of_local_dofs, temp, 0, mrComm);							
					
					//create and fill the graph of the matrix --> the temp array is reused here with a different meaning
					Epetra_FECrsGraph Agraph(Copy, my_map, mguess_row_size);
					//int ierr;
					Element::EquationIdVectorType EquationId;
					// assemble all elements
					for (typename ElementsArrayType::ptr_iterator it=rElements.ptr_begin(); it!=rElements.ptr_end(); it++)
					{
                                            if( ! (*it)->GetValue( IS_INACTIVE ) )
                                            {

						//TODO! this should go through the scheme!!!!!!!!!1
						(*it)->EquationIdVector(EquationId,CurrentProcessInfo);
						
						//filling the list of active global indices (non fixed)
						unsigned int num_active_indices = 0;
						for(unsigned int i=0; i<EquationId.size(); i++)
							//if ( EquationId[i] < BaseType::mEquationSystemSize ) //check!!!
							{
								temp[num_active_indices] =  EquationId[i];
								num_active_indices += 1;
// KRATOS_WATCH(temp[i]);
							}

// KRATOS_WATCH(" ");
						if(num_active_indices != 0)
						{
							int ierr = Agraph.InsertGlobalIndices(num_active_indices,temp,num_active_indices, temp);
//							if(ierr != 0) KRATOS_ERROR(std::logic_error,"Epetra failure found in Agraph.InsertGlobalIndices --> ln 942","");
						}						
                                            }
					}
// KRATOS_WATCH("assemble conditions");		
					// assemble all conditions
					for (typename ConditionsArrayType::ptr_iterator it=rConditions.ptr_begin(); it!=rConditions.ptr_end(); it++)
					{
                                            if( ! (*it)->GetValue( IS_INACTIVE ) )
                                            {

						//TODO! this should go through the scheme!!!!!!!!!1
						(*it)->EquationIdVector(EquationId,CurrentProcessInfo);
						
						//filling the list of active global indices (non fixed)
						unsigned int num_active_indices = 0;
						for(unsigned int i=0; i<EquationId.size(); i++)
							//if ( EquationId[i] < BaseType::mEquationSystemSize ) //check!!!
							{
								temp[num_active_indices] =  EquationId[i];
								num_active_indices += 1;

// KRATOS_WATCH(temp[i]);
							}

						if(num_active_indices != 0)
						{
							int ierr = Agraph.InsertGlobalIndices(num_active_indices,temp,num_active_indices, temp);
//							if(ierr != 0) KRATOS_ERROR(std::logic_error,"Epetra failure found in Agraph.InsertGlobalIndices --> ln 966","");
						}
                                            }
					}

					//finalizing graph construction
					int graph_assemble_ierr = Agraph.GlobalAssemble();
					if(graph_assemble_ierr != 0) KRATOS_ERROR(std::logic_error,"Epetra failure found in Agraph.GlobalAssemble()","");

// KRATOS_WATCH(Agraph);
				
					//generate a new matrix pointer according to this graph
					TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(Copy,Agraph) );	
					pA.swap(pNewA);
// KRATOS_WATCH(*pA);


					//generate new vector pointers according to the given map
					if( pb == NULL || TSparseSpace::Size(*pb) != BaseType::mEquationSystemSize)
					{
						TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(my_map) );
						pb.swap(pNewb);	
					}
					if( pDx == NULL || TSparseSpace::Size(*pDx) != BaseType::mEquationSystemSize)
					{
						TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(my_map) );
						pDx.swap(pNewDx);	
					}
					if( BaseType::mpReactionsVector == NULL) //if the pointer is not initialized initialize it to an empty matrix
					{
						TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(my_map) );
						BaseType::mpReactionsVector.swap(pNewReactionsVector); 
					}

KRATOS_WATCH( TSparseSpace::Size1(*pA) );
KRATOS_WATCH( TSparseSpace::Size(*pb) );
					delete [] temp;



					
				}
				else
				{
					if(TSparseSpace::Size1(*pA) == 0 || TSparseSpace::Size1(*pA) != BaseType::mEquationSystemSize || TSparseSpace::Size2(*pA) != BaseType::mEquationSystemSize)
					{
						KRATOS_ERROR(std::logic_error,"it should not come here resizing is not allowed this way!!!!!!!! ... ","");
					}
				}

				//


				//if needed resize the vector for the calculation of reactions
				if(BaseType::mCalculateReactionsFlag == true)
				{

					KRATOS_ERROR(std::logic_error,"calculation of reactions not yet implemented with Trilinos","");
				}

std::cout << "finished ResizeAndInitializeVectors" << std::endl;

				KRATOS_CATCH("")

		}


		//**************************************************************************
		//**************************************************************************
		void InitializeSolutionStep(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
			KRATOS_TRY
			KRATOS_CATCH("")
		}

		//**************************************************************************
		//**************************************************************************
		void FinalizeSolutionStep(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
		}


		//**************************************************************************
		//**************************************************************************
		void CalculateReactions(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{

		}

		void BuildLHS_CompleteOnFreeRows(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemMatrixType& A)
		{
			KRATOS_ERROR(std::logic_error,"method BuildLHS_CompleteOnFreeRows not implemented in Trilinos Builder And Solver ","");
		}

		//**************************************************************************
		//**************************************************************************
		void ApplyDirichletConditions(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{

			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);	
			
			int NumEntries;    // number of nonzero entries extracted	

			
			for (typename DofsArrayType::iterator dof_it = BaseType::mDofSet.begin(); dof_it != BaseType::mDofSet.end(); ++dof_it)
			{ 
			    if(dof_it->GetSolutionStepValue(PARTITION_INDEX) == rank)	
			    {
				if (dof_it->IsFixed()) 
				{

					int GlobalRow = dof_it->EquationId();  // row to extract
					int Length = A.NumGlobalEntries(dof_it->EquationId());  // length of Values and Indices
					
   					double* Values = new double[Length];     // extracted values for this row     
					int* Indices = new int[Length];          // extracted global column indices for the corresponding values

					A.ExtractGlobalRowCopy(GlobalRow, Length, NumEntries, Values, Indices);
		
					// put 0.0 in each row A[ii] and 1.0 on the diagonal
					for (int ii=0; ii<Length; ii++)
					{
						if (Indices[ii] == GlobalRow) 
							Values[ii]=1.0; 
						else
							Values[ii]=0.0;
					}

					A.ReplaceGlobalValues(GlobalRow, Length, Values, Indices);
					
					// redo better !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					int* gb= new int[1]; 
					gb[0]=GlobalRow;
					A.ReplaceGlobalValues(Length, Indices, 1, gb, Values);

					double* bb=new double[1];
					bb[0]=0.0;
					
					// put 0.0 in vector b[GlobalRow] if GlobalRow is a fixed dof
					b.ReplaceGlobalValues(1,gb,bb);

					delete [] Values;	
					delete [] Indices;
					delete [] gb;
					delete [] bb;
					
			        }
			    }	
			
			} 
			

		}

		//**************************************************************************
		//**************************************************************************
		void ApplyPointLoads(
			typename TSchemeType::Pointer pScheme,
			ModelPart& r_model_part,
			TSystemVectorType& b)
		{}

		/** 
		this function is intended to be called at the end of the solution step to clean up memory
		storage not needed 
		*/
		void Clear() 
		{
			this->mDofSet = DofsArrayType();
			//this->mReactionsVector = TSystemVectorType();

			if (this->GetEchoLevel()>0)
			{

				KRATOS_WATCH("TrilinosBuilderAndSolverMLDeactivation Clear Function called");
			}
		}


		/*@} */
		/**@name Operations */
		/*@{ */


		/*@} */  
		/**@name Access */
		/*@{ */


		/*@} */
		/**@name Inquiry */
		/*@{ */


		/*@} */      
		/**@name Friends */
		/*@{ */


		/*@} */

	protected:
		/**@name Protected static Member Variables */
		/*@{ */


		/*@} */
		/**@name Protected member Variables */
		/*@{ */
		Epetra_MpiComm& mrComm;
		int mguess_row_size;
		int mFirstMyId;
		int mLastMyId;


		/*@} */
		/**@name Protected Operators*/
		/*@{ */

		/*@} */
		/**@name Protected Operations*/
		/*@{ */


		/*@} */
		/**@name Protected  Access */
		/*@{ */


		/*@} */     
		/**@name Protected Inquiry */
		/*@{ */


		/*@} */   
		/**@name Protected LifeCycle */  
		/*@{ */



		/*@} */    
		

		//virtual void GenerateNullSpace(TSystemMatrixType& A,ModelPart& r_model_part , double* nullsp , Teuchos::RCP<vector<double> >&  ns )
		virtual void GenerateNullSpace(TSystemMatrixType& A,ModelPart& r_model_part , double* nullsp , boost::shared_ptr<vector<double> >&  ns,
						int& numdf, int& dimns )
			{};


	private:
		/**@name Static Member Variables */
		/*@{ */


		/*@} */
		/**@name Member Variables */
		/*@{ */

		unsigned int mLocalSystemSize;

		/*@} */
		/**@name Private Operators*/
		/*@{ */


		/*@} */
		/**@name Private Operations*/
		/*@{ */


		//**************************************************************************
		void AssembleLHS_CompleteOnFreeRows(
			TSystemMatrixType& A,
			LocalSystemMatrixType& LHS_Contribution,
			Element::EquationIdVectorType& EquationId
			)
		{
		  KRATOS_ERROR(std::logic_error, "This method is not implemented for Trilinos", "");
		}		

		/*@} */
		/**@name Private  Access */
		/*@{ */


		/*@} */     
		/**@name Private Inquiry */
		/*@{ */


		/*@} */   
		/**@name Un accessible methods */
		/*@{ */


		/*@} */   

	}; /* Class TrilinosBuilderAndSolverMLDeactivation */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_BUILDER_AND_SOLVER_ML_DEACTIVATION  defined */

