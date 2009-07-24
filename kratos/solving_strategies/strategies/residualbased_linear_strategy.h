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
*   Revision:            $Revision: 1.7 $
*
* ***********************************************************/


#if !defined(KRATOS_NEW_RESIDUALBASED_LINEAR_STRATEGY )
#define  KRATOS_NEW_RESIDUALBASED_LINEAR_STRATEGY


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "solving_strategies/strategies/solving_strategy.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

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
	
	/// Short class definition.
	/**   Detail class definition.
	
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
			 class TDenseSpace, //= DenseSpace<double>,
			 class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
			>
		class ResidualBasedLinearStrategy 
		: public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
    {
    public:
		/**@name Type Definitions */       
		/*@{ */

		/** Counted pointer of ClassName */
		//typedef boost::shared_ptr< ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;
		KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedLinearStrategy );
		
    typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

	typedef typename BaseType::TSchemeType TSchemeType;
	typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
		
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
		ResidualBasedLinearStrategy(
			ModelPart& model_part, 
			typename TSchemeType::Pointer pScheme,
			typename TLinearSolver::Pointer pNewLinearSolver,
			bool CalculateReactionFlag = false,
			bool ReformDofSetAtEachStep = false,
			bool CalculateNormDxFlag = false,
			bool MoveMeshFlag = false
			)
			: SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag)
		{
			KRATOS_TRY

			mCalculateReactionsFlag = CalculateReactionFlag;
			mReformDofSetAtEachStep = ReformDofSetAtEachStep;
			mCalculateNormDxFlag = CalculateNormDxFlag;


			//saving the scheme
			mpScheme = pScheme;

			//saving the linear solver
			mpLinearSolver = pNewLinearSolver;

			//setting up the default builder and solver
			mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer
				(
				new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(mpLinearSolver)
				);

			//set flag to start correcty the calculations
			mSolutionStepIsInitialized = false;
			mInitializeWasPerformed = false;

			//tells to the builder and solver if the reactions have to be Calculated or not
			GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);
			
			//tells to the Builder And Solver if the system matrix and vectors need to
			//be reshaped at each step or not
			GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

			//set EchoLevel to the default value (only time is displayed)
			this->SetEchoLevel(1);

			//by default the matrices are rebuilt at each solution step
			BaseType::SetRebuildLevel(1);

			KRATOS_CATCH("")
		}

		//constructor specifying the builder and solver
		ResidualBasedLinearStrategy(
			ModelPart& model_part, 
			typename TSchemeType::Pointer pScheme,
			typename TLinearSolver::Pointer pNewLinearSolver,
			typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
			bool CalculateReactionFlag = false,
			bool ReformDofSetAtEachStep = false,
			bool CalculateNormDxFlag = false,
			bool MoveMeshFlag = false
			)
			: SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part,MoveMeshFlag)
		{
			KRATOS_TRY

			mCalculateReactionsFlag = CalculateReactionFlag;
			mReformDofSetAtEachStep = ReformDofSetAtEachStep;
			mCalculateNormDxFlag = CalculateNormDxFlag;

			//saving the scheme
			mpScheme = pScheme;

			//saving the linear solver
			mpLinearSolver = pNewLinearSolver;

			//setting up the  builder and solver
			mpBuilderAndSolver = pNewBuilderAndSolver;

			//set flag to start correcty the calculations
			mSolutionStepIsInitialized = false;
			mInitializeWasPerformed = false;

			//tells to the builder and solver if the reactions have to be Calculated or not
			GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);
			
			//tells to the Builder And Solver if the system matrix and vectors need to
			//be reshaped at each step or not
			GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

			//set EchoLevel to the default value (only time is displayed)
			this->SetEchoLevel(1);

			//by default the matrices are rebuilt at each solution step
			BaseType::SetRebuildLevel(1);

			KRATOS_CATCH("")
		}

		/** Destructor.
		*/
		virtual ~ResidualBasedLinearStrategy() {}
		
		/** Destructor.
		*/
		
		//Set and Get Scheme ... containing Builder, Update and other
		void SetScheme(typename TSchemeType::Pointer pScheme ) {mpScheme = pScheme;};
		typename TSchemeType::Pointer GetScheme() {return mpScheme;};
		
		//Set and Get the BuilderAndSolver
		void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver ) {mpBuilderAndSolver = pNewBuilderAndSolver;};
		typename TBuilderAndSolverType::Pointer GetBuilderAndSolver() {return mpBuilderAndSolver;};

		
		void SetCalculateReactionsFlag(bool CalculateReactionsFlag) 
			{
				mCalculateReactionsFlag = CalculateReactionsFlag;
				GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);
			}
		bool GetCalculateReactionsFlag() {return mCalculateReactionsFlag;}
		
		void SetReformDofSetAtEachStepFlag(bool flag) 
		{
			mReformDofSetAtEachStep = flag;
			GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
		}
		bool GetReformDofSetAtEachStepFlag() {return mReformDofSetAtEachStep;}
		
		//level of echo for the solving strategy
		// 0 -> mute... no echo at all
		// 1 -> printing time and basic informations
		// 2 -> printing linear solver data
		// 3 -> Print of debug informations:
		//		Echo of stiffness matrix, Dx, b...
		void SetEchoLevel(int Level) 
		{
			BaseType::SetEchoLevel( Level );
			GetBuilderAndSolver()->SetEchoLevel(Level);	
		}
		
		
		//*********************************************************************************
		/**OPERATIONS ACCESSIBLE FROM THE INPUT:*/
		
		/**
		operation to predict the solution ... if it is not called a trivial predictor is used in which the 
		values of the solution step of interest are assumed equal to the old values
		*/
		void Predict()
		{
			KRATOS_TRY
				//OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
				//if the operations needed were already performed this does nothing
			//if(mInitializeWasPerformed == false)
			//{
			//	Initialize();
			//	mInitializeWasPerformed = true;
			//}

			////initialize solution step
			//if (mSolutionStepIsInitialized == false)
			//	InitializeSolutionStep();


			TSystemMatrixType& mA = *mpA;
			TSystemVectorType& mDx = *mpDx;
			TSystemVectorType& mb = *mpb;

			DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();

			this->GetScheme()->Predict(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

			KRATOS_CATCH("")
		}
		
		//*********************************************************************************
		/** 
		the problem of interest is solved
		a double containing norm(Dx) is returned if CalculateNormDxFlag == true, else 0 is returned
		*/
		//**********************************************************************
		double Solve()
		{
			KRATOS_TRY

			//pointers needed in the solution
			typename TSchemeType::Pointer pScheme = GetScheme();
			typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
			int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

			ProcessInfo& pCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

			//OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
			//if the operations needed were already performed this does nothing
			if(mInitializeWasPerformed == false)
			{
				Initialize();
				mInitializeWasPerformed = true;
			}

			//prints informations about the current time
			if (BaseType::GetEchoLevel()!=0 && rank == 0)
			{
				std::cout << " " << std::endl;
				std::cout << "CurrentTime = " << pCurrentProcessInfo[TIME] << std::endl;
			}		
			
			//initialize solution step
			if (mSolutionStepIsInitialized == false)
			{
				InitializeSolutionStep();
				mSolutionStepIsInitialized = true;
			}

			//updates the database with a prediction of the solution
			Predict();

			TSystemMatrixType& mA = *mpA;
			TSystemVectorType& mDx = *mpDx;
			TSystemVectorType& mb = *mpb;

				
			if(BaseType::mRebuildLevel >0 || BaseType::mStiffnessMatrixIsBuilt == false)
			{
				TSparseSpace::SetToZero(mA);
				TSparseSpace::SetToZero(mDx);
				TSparseSpace::SetToZero(mb);
				pBuilderAndSolver->BuildAndSolve(pScheme,BaseType::GetModelPart(),mA,mDx,mb);
				BaseType::mStiffnessMatrixIsBuilt = true;
			}
			else
			{
				TSparseSpace::SetToZero(mDx); 
				TSparseSpace::SetToZero(mb); 
				pBuilderAndSolver->BuildRHSAndSolve(pScheme,BaseType::GetModelPart(),mA,mDx,mb);
			}

			if (BaseType::GetEchoLevel()==3) //if it is needed to print the debug info
			{
				std::cout << "SystemMatrix = " << mA << std::endl;
				std::cout << "solution obtained = " << mDx << std::endl;
				std::cout << "RHS  = " << mb << std::endl;
			}

			//update results
			DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();
			pScheme->Update(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

			//move the mesh if needed
			if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();

			//calculate if needed the norm of Dx
			double normDx = 0.00;
			if(mCalculateNormDxFlag == true)
			{
				normDx = TSparseSpace::TwoNorm(mDx);
			}

			//calculate reactions if required
			if (mCalculateReactionsFlag ==true)
			{
				pBuilderAndSolver->CalculateReactions(pScheme,BaseType::GetModelPart(),mA,mDx,mb);
			}

			//Finalisation of the solution step, 
			//operations to be done after achieving convergence, for example the 
			//Final Residual Vector (mb) has to be saved in there 
			//to avoid error accumulation
			pScheme->FinalizeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
			pBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
			
			//deallocate the systemvectors if needed
			if (mReformDofSetAtEachStep == true) 
			{
				if(rank == 0) std::cout << "Clearing System" << std::endl;
				this->Clear();
				//std::cout << "Clearing System" << std::endl;
				//TSparseSpace::ClearData(mA);
				//TSparseSpace::ClearData(mDx);
				//TSparseSpace::ClearData(mb);
			}

			//Cleaning memory after the solution
			pScheme->Clean();
			
			//reset flags for next step
			mSolutionStepIsInitialized = false;

			return normDx;
			
			KRATOS_CATCH("")
				
	}
	
		TSystemMatrixType& GetSystemMatrix()
		    {
			TSystemMatrixType& mA = *mpA;

			return mA;
		    }
	
	//*********************************************************************************
	/** 
	this operations should be called before printing the results when non trivial results (e.g. stresses)
	need to be calculated given the solution of the step
	
	  This operations should be called only when needed, before printing as it can involve a non negligible cost
	*/
	void CalculateOutputData()
	{
		TSystemMatrixType& mA = *mpA;
		TSystemVectorType& mDx = *mpDx;
		TSystemVectorType& mb = *mpb;

		DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();
		GetScheme()->CalculateOutputData(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);		
	}
	
	//*********************************************************************************
	
	
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
	
    protected:
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
		
    private:
        /**@name Static Member Variables */
        /*@{ */
        
        
        /*@} */
        /**@name Member Variables */
        /*@{ */

		typename TLinearSolver::Pointer mpLinearSolver;

		typename TSchemeType::Pointer mpScheme;     
  
		typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;
	
		TSystemVectorPointerType mpDx;
		TSystemVectorPointerType mpb;
		TSystemMatrixPointerType mpA;
		
		/** 
		Flag telling if it is needed to reform the DofSet at each
		solution step or if it is possible to form it just once
		- true  => reforme at each time step
		- false => form just one (more efficient)
		
		  Default = false
		*/
		bool mReformDofSetAtEachStep;

		//calculates if required the norm of the correction term Dx
		bool mCalculateNormDxFlag; 
		
		/** 
		Flag telling if it is needed or not to compute the reactions
		
		  default = true
		*/
		bool mCalculateReactionsFlag;
		
		bool mSolutionStepIsInitialized;

		bool mInitializeWasPerformed;
		
		///default = 30
		unsigned int mMaxIterationNumber;
		
		
		
		/*@} */
		/**@name Private Operators*/
		/*@{ */
		//**********************************************************************
		//**********************************************************************
		void Initialize()
		{
			KRATOS_TRY

			if(BaseType::GetEchoLevel()>2)
				std::cout << "entering in the  Initialize of the ResidualBasedLinearStrategy" << std::endl;
			
			//pointers needed in the solution
			typename TSchemeType::Pointer pScheme = GetScheme();

			//Initialize The Scheme - OPERATIONS TO BE DONE ONCE
			if (pScheme->SchemeIsInitialized() == false) 
				pScheme->Initialize(BaseType::GetModelPart());

                        //Initialize The Elements - OPERATIONS TO BE DONE ONCE
			if (pScheme->ElementsAreInitialized() == false) 
				pScheme->InitializeElements(BaseType::GetModelPart());

                        if(BaseType::GetEchoLevel()>2)
				std::cout << "exiting the  Initialize of the ResidualBasedLinearStrategy" << std::endl;

						
			KRATOS_CATCH("")
		}
		
		
		//**********************************************************************
		//**********************************************************************
		void InitializeSolutionStep()
		{
			KRATOS_TRY
			
			typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
			typename TSchemeType::Pointer pScheme = GetScheme();
                        
                        int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

			if(BaseType::GetEchoLevel()>2 && rank==0)
				std::cout << "entering in the  InitializeSolutionStep of the ResidualBasedLinearStrategy" << std::endl;

				
			//loop to reform the dofset
			boost::timer system_construction_time;
			if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false || 
				mReformDofSetAtEachStep == true )
			{
				boost::timer setup_dofs_time;
				//setting up the list of the DOFs to be solved
				pBuilderAndSolver->SetUpDofSet(pScheme,BaseType::GetModelPart());
				if(BaseType::GetEchoLevel()>0 && rank == 0)
					std::cout << "setup_dofs_time : " << setup_dofs_time.elapsed() << std::endl;
	  			
				//shaping correctly the system 
				boost::timer setup_system_time;
				pBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());	
				if(BaseType::GetEchoLevel()>0 && rank == 0)
					std::cout << "setup_system_time : " << setup_system_time.elapsed() << std::endl;

				//setting up the Vectors involved to the correct size 
				boost::timer system_matrix_resize_time;
				pBuilderAndSolver->ResizeAndInitializeVectors(mpA,mpDx,mpb,BaseType::GetModelPart().Elements(),BaseType::GetModelPart().Conditions(),BaseType::GetModelPart().GetProcessInfo());
				if(BaseType::GetEchoLevel()>0 && rank == 0)
					std::cout << "system_matrix_resize_time : " << system_matrix_resize_time.elapsed() << std::endl;
			}
			if(BaseType::GetEchoLevel()>0 && rank == 0)
				std::cout << "System Construction Time : " << system_construction_time.elapsed() << std::endl;


			TSystemMatrixType& mA = *mpA;
			TSystemVectorType& mDx = *mpDx;
			TSystemVectorType& mb = *mpb;

			//initial operations ... things that are constant over the Solution Step
			pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
		
			//initial operations ... things that are constant over the Solution Step
			pScheme->InitializeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
		
			
			KRATOS_CATCH("")
		}	
		
		//**********************************************************************
		//**********************************************************************
		void MaxIterationsExceeded()
		{
			std::cout << "***************************************************" << std::endl;
			std::cout << "******* ATTENTION: max iterations exceeded ********" << std::endl;
			std::cout << "***************************************************" << std::endl;
		}

		//**********************************************************************
		//**********************************************************************
		void Clear()
		{
			KRATOS_TRY
				std::cout << "strategy Clear function used" << std::endl;

		
		
			
			
			if(mpA != NULL) {
			  TSystemMatrixType& mA = *mpA;
			  SparseSpaceType::Clear(mpA);
			  SparseSpaceType::Resize(mA,0,0);
			}			  
			if(mpDx != NULL) {
			  TSystemVectorType& mDx = *mpDx;
			  SparseSpaceType::Clear(mpDx);
			  SparseSpaceType::Resize(mDx,0);
			}
			if(mpb != NULL) {
			  TSystemVectorType& mb = *mpb;
			  SparseSpaceType::Clear(mpb);
			  SparseSpaceType::Resize(mb,0);}
			
			//setting to zero the internal flag to ensure that the dof sets are recalculated
			GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);

			GetBuilderAndSolver()->Clear();
				
				
			KRATOS_CATCH("");
		}		
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
		ResidualBasedLinearStrategy(const ResidualBasedLinearStrategy& Other);
		
        
        /*@} */   
        
    }; /* Class ResidualBasedLinearStrategy */
	
	/*@} */
	
	/**@name Type Definitions */       
	/*@{ */
	
	
	/*@} */
	
}  /* namespace Kratos.*/

#endif /* KRATOS_NEW_RESIDUALBASED_LINEAR_STRATEGY  defined */

