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
*   Date:                $Date: 2008-11-10 14:23:33 $
*   Revision:            $Revision: 1.12 $
*
* ***********************************************************/


#if !defined(KRATOS_EXPLICIT_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY )
#define  KRATOS_EXPLICIT_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY

/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/matrix_market_interface.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//default builder and solver
#include "custom_strategies/builder_and_solvers/explicit_residualbased_builder.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"


namespace Kratos
{


	template<class TSparseSpace,
	class TDenseSpace, // = DenseSpace<double>,
	class TLinearSolver> //= LinearSolver<TSparseSpace,TDenseSpace>

	class ExplicitResidualBasedNewtonRaphsonStrategy : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
	{
	public:
		/**@name Type Definitions */       
		/*@{ */
		typedef ConvergenceCriteria<TSparseSpace,TDenseSpace> TConvergenceCriteriaType;

		/** Counted pointer of ClassName */
		//typedef boost::shared_ptr< ResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;
		KRATOS_CLASS_POINTER_DEFINITION( ExplicitResidualBasedNewtonRaphsonStrategy );

		typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;
                typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver> NewtonBaseType;
		typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

		typedef typename BaseType::TDataType TDataType;
		
		typedef TSparseSpace SparseSpaceType;
		
		typedef typename BaseType::TSchemeType TSchemeType;

		//typedef typename BaseType::DofSetType DofSetType;

		typedef typename BaseType::DofsArrayType DofsArrayType;

		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

		typedef typename BaseType::TSystemVectorType TSystemVectorType;

		typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

		typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

		typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
		typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;


		/*@} */
		/**@name Life Cycle 
		*/    
		/*@{ */

		/** Constructor.
		*/
		ExplicitResidualBasedNewtonRaphsonStrategy(
			ModelPart& model_part, 
			typename TSchemeType::Pointer pScheme,
			typename TLinearSolver::Pointer pNewLinearSolver,
			typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
			int MaxIterations = 30,
			bool CalculateReactions = false,
			bool ReformDofSetAtEachStep = false,
			bool MoveMeshFlag = false
			)
			: ResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, 
				pScheme,
				pNewLinearSolver,
				pNewConvergenceCriteria,
				MaxIterations,
				CalculateReactions,
				ReformDofSetAtEachStep,
				MoveMeshFlag)
		{
			KRATOS_TRY
				//set flags to default values
				NewtonBaseType::SetMaxIterationNumber(MaxIterations);
			NewtonBaseType::mCalculateReactionsFlag = CalculateReactions;
			
			
			NewtonBaseType::mReformDofSetAtEachStep = ReformDofSetAtEachStep;

			//saving the convergence criteria to be used 
			NewtonBaseType::mpConvergenceCriteria = pNewConvergenceCriteria;

			//saving the scheme
			NewtonBaseType::mpScheme = pScheme;

			//saving the linear solver
			NewtonBaseType::mpLinearSolver = pNewLinearSolver;

			//setting up the default builder and solver
			NewtonBaseType::mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer
				(
				new ExplicitResidualBasedBuilder<TSparseSpace,TDenseSpace,TLinearSolver>(NewtonBaseType::mpLinearSolver)
				);

			//set flags to start correcty the calculations
			NewtonBaseType::mSolutionStepIsInitialized = false;
			NewtonBaseType::mInitializeWasPerformed = false;

			//tells to the builder and solver if the reactions have to be Calculated or not
			NewtonBaseType::GetBuilderAndSolver()->SetCalculateReactionsFlag(NewtonBaseType::mCalculateReactionsFlag);

			//tells to the Builder And Solver if the system matrix and vectors need to
			//be reshaped at each step or not
			NewtonBaseType::GetBuilderAndSolver()->SetReshapeMatrixFlag(NewtonBaseType::mReformDofSetAtEachStep);

			//set EchoLevel to the default value (only time is displayed)
			NewtonBaseType::SetEchoLevel(1);

			//by default the matrices are rebuilt at each iteration
			NewtonBaseType::SetRebuildLevel(2);

			KRATOS_CATCH("")			
		}


		/** Destructor.
		*/
		virtual ~ExplicitResidualBasedNewtonRaphsonStrategy() {}

		/** Destructor.
		*/
		void Predict()
		{
			KRATOS_TRY
				//OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
				//if the operations needed were already performed this does nothing
				if(NewtonBaseType::mInitializeWasPerformed == false)
					NewtonBaseType::Initialize();
KRATOS_WATCH("0000000000000000000000");
			//initialize solution step
			if (NewtonBaseType::mSolutionStepIsInitialized == false)
				NewtonBaseType::InitializeSolutionStep();
KRATOS_WATCH("111111111111111111111111111111");
			DofsArrayType& rDofSet = NewtonBaseType::GetBuilderAndSolver()->GetDofSet();

			TSystemMatrixType& mA = *NewtonBaseType::mpA;
			TSystemVectorType& mDx = *NewtonBaseType::mpDx;
			TSystemVectorType& mb = *NewtonBaseType::mpb;


			NewtonBaseType::GetScheme()->Predict(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
KRATOS_WATCH("222222222222222222222222222222222222");
			//move the mesh if needed
			if(this->MoveMeshFlag() == true) BaseType::MoveMesh();

			KRATOS_CATCH("")
		}


		//*********************************************************************************
		/** 
		the problem of interest is solved
		*/
		//**********************************************************************
		double Solve()
		{
			KRATOS_TRY


				//pointers needed in the solution
				typename TSchemeType::Pointer pScheme = NewtonBaseType::GetScheme();
			typename TBuilderAndSolverType::Pointer pBuilderAndSolver = NewtonBaseType::GetBuilderAndSolver();

			DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();

			//OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
			//if the operations needed were already performed this does nothing
			if(NewtonBaseType::mInitializeWasPerformed == false)
				NewtonBaseType::Initialize();

			//prints informations about the current time
			if (this->GetEchoLevel()!=0)
			{
				std::cout << " " << std::endl;
				std::cout << "CurrentTime = " << BaseType::GetModelPart().GetProcessInfo()[TIME] << std::endl;
			}

			//updates the database with a prediction of the solution
KRATOS_WATCH("BEFORE PREDICt");
			Predict();
KRATOS_WATCH("AFTER PREDICt");
			//initialize solution step
			if (NewtonBaseType::mSolutionStepIsInitialized == false)
				NewtonBaseType::InitializeSolutionStep();
// 
			TSystemMatrixType& mA = *NewtonBaseType::mpA;
			TSystemVectorType& mDx = *NewtonBaseType::mpDx;
			TSystemVectorType& mb = *NewtonBaseType::mpb;

// KRATOS_WATCH(mA);
// KRATOS_WATCH(mDx);
// KRATOS_WATCH(mb);


// 			//initializing the parameters of the Newton-Raphson cicle
			unsigned int iteration_number=1; 
			BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;	
// //			BaseType::GetModelPart().GetProcessInfo().SetNonLinearIterationNumber(iteration_number);
 			bool is_converged = false;
			//bool ResidualIsUpdated = false;
			pScheme->InitializeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);
			is_converged = NewtonBaseType::mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

			//function to perform the building 


KRATOS_WATCH("BEFORE BuildRHS");
			pBuilderAndSolver->BuildRHS(pScheme,BaseType::GetModelPart(),mb);
KRATOS_WATCH("AfTER BuildRHS");



			//update results
			pScheme->Update(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
			//move the mesh if needed
			if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();

			pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);
			is_converged = NewtonBaseType::mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
 


			//Iteration Cicle... performed only for NonLinearProblems
			while(	is_converged == false &&
				iteration_number++<NewtonBaseType::mMaxIterationNumber) 
			{
				//setting the number of iteration
				BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;	

				pScheme->InitializeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);

				is_converged = NewtonBaseType::mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

				//call the linear system solver to find the correction mDx for the 
				//it is not called if there is no system to solve

				KRATOS_WATCH(iteration_number);
				pBuilderAndSolver->BuildRHS(pScheme,BaseType::GetModelPart(),mb);			


				//Updating the results stored in the database
// 				rDofSet = pBuilderAndSolver->GetDofSet();
				pScheme->Update(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

				//move the mesh if needed
				if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();

				pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);
				is_converged = NewtonBaseType::mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

// 				ResidualIsUpdated = false;

			}



			//Finalisation of the solution step, 
			//operations to be done after achieving convergence, for example the 
			//Final Residual Vector (mb) has to be saved in there 
			//to avoid error accumulation
			pScheme->FinalizeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
			pBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);

			//Cleaning memory after the solution
			pScheme->Clean();

			//reset flags for next step
			NewtonBaseType::mSolutionStepIsInitialized = false;



			return 0.00;

			KRATOS_CATCH("")

		}


		/*@} */
		/**@name Operators 
		*/  
		/*@{ */

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

	private:
		/**@name Protected static Member Variables */
		/*@{ */


		/*@} */
		/**@name Protected member Variables */
		/*@{ */


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

	protected:
		/**@name Static Member Variables */
		/*@{ */


		/*@} */
		/**@name Member Variables */
		/*@{ */

		//**********************************************************************
		//**********************************************************************
// 		void InitializeSolutionStep()
// 		{
// 			KRATOS_TRY
// 
// 				typename TBuilderAndSolverType::Pointer pBuilderAndSolver = NewtonBaseType::GetBuilderAndSolver();
// 			typename TSchemeType::Pointer pScheme = NewtonBaseType::GetScheme();
// 
// KRATOS_WATCH("EXPLICITTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
// KRATOS_WATCH("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
// 			//setting up the Vectors involved to the correct size 
// 			pBuilderAndSolver->ResizeAndInitializeVectors(mpA,mpDx,mpb,BaseType::GetModelPart().Elements(),BaseType::GetModelPart().Conditions(),BaseType::GetModelPart().GetProcessInfo());
// KRATOS_WATCH("QQQQQQQQQQQQQQQQQQQQQQQQQ");
// 			TSystemMatrixType& mA = *NewtonBaseType::mpA;
// 			TSystemVectorType& mDx = *NewtonBaseType::mpDx;
// 			TSystemVectorType& mb = *NewtonBaseType::mpb;
// 
// 
// 			//initial operations ... things that are constant over the Solution Step
// 			pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
// KRATOS_WATCH("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
// 			//initial operations ... things that are constant over the Solution Step
// 			pScheme->InitializeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
// KRATOS_WATCH("SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS");
// 			NewtonBaseType::mSolutionStepIsInitialized = true;			
// 
// 			KRATOS_CATCH("")
// 		}	


		//**********************************************************************

		/*@} */
		/**@name Private Operations*/
		/*@{ */


		/*@} */
		/**@name Private  Access */
		/*@{ */


		/*@} */     
		/**@name Private Inquiry */
		/*@{ */


		/*@} */   
		/**@name Un accessible methods */
		/*@{ */

		/** Copy constructor.
		*/
		ExplicitResidualBasedNewtonRaphsonStrategy(const ExplicitResidualBasedNewtonRaphsonStrategy& Other);


		/*@} */   

	}; /* Class ResidualBasedNewtonRaphsonStrategy */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY  defined */

