/* *********************************************************   
*          
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2008-11-10 14:23:32 $
*   Revision:            $Revision: 1.12 $
*
* ***********************************************************/


#if !defined(KRATOS_LAP_MODIFIED_LINEAR_STRATEGY)
#define  KRATOS_LAP_MODIFIED_LINEAR_STRATEGY


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"

#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "ULF_application.h" 


#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"


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
	template<unsigned int TDim, class TSparseSpace,
			 class TDenseSpace, //= DenseSpace<double>,
			 class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
			>
		class LapModifiedLinearStrategy 
		: public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
    {
    public:
		/**@name Type Definitions */       
		/*@{ */

		/** Counted pointer of ClassName */
		typedef boost::shared_ptr< LapModifiedLinearStrategy<TDim, TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;
		//KRATOS_CLASS_POINTER_DEFINITION( LapModifiedLinearStrategy );
		
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
		
		typedef ModelPart::NodesContainerType NodesArrayType;
		
		/*@} */
		/**@name Life Cycle 
		*/    
		/*@{ */
		
		/** Constructor.
		*/
		LapModifiedLinearStrategy(
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
		LapModifiedLinearStrategy(
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
		virtual ~LapModifiedLinearStrategy() {}
		
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

			TSystemMatrixType& mD = mpD;
				
			if(BaseType::mRebuildLevel >0 || BaseType::mStiffnessMatrixIsBuilt == false)
			{
				TSparseSpace::SetToZero(mA);
				TSparseSpace::SetToZero(mDx);
				TSparseSpace::SetToZero(mb);

				TSparseSpace::SetToZero(mD);

				//pBuilderAndSolver->BuildAndSolve(pScheme,BaseType::GetModelPart(),mA,mDx,mb);
				
				pBuilderAndSolver->Build(pScheme,BaseType::GetModelPart(),mA,mb);
				//ADD HERE THE INTERFACE LAPLACIAN
				
				ConstructMatrixStructure_FluidDivergenceMatrixD( mD,  BaseType::GetModelPart());
				/*
				double density_str=1000.0;
				BuildAuxiliariesFSI(mD, density_str, BaseType::GetModelPart());
				mA+=mD;
				*/
				pBuilderAndSolver->ApplyDirichletConditions(pScheme,BaseType::GetModelPart(),mA,mDx,mb);
				//solve!
				pBuilderAndSolver->SystemSolve(mA,mDx,mb);
				
				
				//pBuilderAndSolver->BuildAndSolve(pScheme,BaseType::GetModelPart(),mA,mDx,mb);
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
	double GetResidualNorm()
	{
		if(TSparseSpace::Size(*mpb) !=0)
			return TSparseSpace::TwoNorm(*mpb);
		else
			return 0.0;		
		
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
	
		TSystemMatrixType mpD;

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
				std::cout << "entering in the  Initialize of the LapModifiedLinearStrategy" << std::endl;
			
			//pointers needed in the solution
			typename TSchemeType::Pointer pScheme = GetScheme();

			//Initialize The Scheme - OPERATIONS TO BE DONE ONCE
			if (pScheme->SchemeIsInitialized() == false) 
				pScheme->Initialize(BaseType::GetModelPart());

                        //Initialize The Elements - OPERATIONS TO BE DONE ONCE
			if (pScheme->ElementsAreInitialized() == false) 
				pScheme->InitializeElements(BaseType::GetModelPart());

                        if(BaseType::GetEchoLevel()>2)
				std::cout << "exiting the  Initialize of the LapModifiedLinearStrategy" << std::endl;

						
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
				std::cout << "entering in the  InitializeSolutionStep of the LapModifiedLinearStrategy" << std::endl;

				
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

			//TSystemMatrixType& mD = mpD;

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
			/*
			if(mpD.size1() != 0 && mpD.size2()!=0) {
			  TSystemMatrixType& mD = mpD;
			  SparseSpaceType::Resize(mD,0,0);}	
			*/
			
			//setting to zero the internal flag to ensure that the dof sets are recalculated
			GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);

			GetBuilderAndSolver()->Clear();
				
				
			KRATOS_CATCH("");
		}	

//////////////////////////////////////////////////			
		void AddFSILaplacian()
		{
		std::cout<<"Now empty";
		}
		///////////////////////////////////////////////////
		void ConstructMatrixStructure_FluidDivergenceMatrixD(
					TSystemMatrixType& mD,  ModelPart& r_model_part)
		{
			KRATOS_TRY
			KRATOS_WATCH("Constructing fluid div matrix")
			unsigned int reduced_dim = this->mpb->size() / TDim;

			mD.resize(reduced_dim,this->mpb->size(),false);

			KRATOS_WATCH(mD)
			
			std::vector<int>  indices;
			indices.reserve(1000);

			//count non zeros
			int total_nnz = 0;
			for (typename NodesArrayType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); ++it)
			{
				//not to add lonely nodes to the system
				if(it->FastGetSolutionStepValue(IS_FLUID)!=0 && (it->GetValue(NEIGHBOUR_NODES)).size() != 0 )
				{
				int count_fluid_neighb=0;
					WeakPointerVector< Node<3> >& neighb_nodes = it->GetValue(NEIGHBOUR_NODES);
					for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin();
						i != neighb_nodes.end(); i++)
					{
					if (i->FastGetSolutionStepValue(IS_FLUID)!=0)
						count_fluid_neighb++;
					}				
				//we add one because we have to account for the contribution of the node itself
				total_nnz += 1 + count_fluid_neighb;//(it->GetValue(NEIGHBOUR_NODES)).size();
				}
			}
			
			mD.reserve(total_nnz* TDim,false);
			KRATOS_WATCH("Constructing fluid div matrix 1")
			
			unsigned int row_index;
			//fill the matrix row by row
			//unsigned int dof_position = r_model_part.NodesBegin()->GetDofPosition(DISPLACEMENT_X);
			unsigned int dof_position = r_model_part.NodesBegin()->GetDofPosition(PRESSURE);
			KRATOS_WATCH("Pressure DOF")
			KRATOS_WATCH(dof_position)
			
			for (typename NodesArrayType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); ++it)
			{
				WeakPointerVector< Node<3> >& neighb_nodes = it->GetValue(NEIGHBOUR_NODES);
				if( neighb_nodes.size() != 0 && it->FastGetSolutionStepValue(IS_FLUID)!=0)
				{
					//first row in the block 
					//row_index = it->GetDof(DISPLACEMENT_X,dof_position).EquationId();
					row_index = it->GetDof(PRESSURE, dof_position).EquationId();
					
					//add id of the current node
					//NOTE: here and in the following we ASSUME that the ids of DISPLACEMENT_X _Y and _Z are sequential
					for(unsigned int kk = 0; kk<TDim; kk++)
					{
						indices.push_back(TDim*row_index + kk);
					}

					//filling and order the first neighbours list
					for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin();
						i != neighb_nodes.end(); i++)
					{
						if ( i->FastGetSolutionStepValue(IS_FLUID)!=0)
						{
							//unsigned int tmp = (i->GetDof(DISPLACEMENT_X,dof_position)).EquationId();
							unsigned int tmp = (i->GetDof(PRESSURE,dof_position)).EquationId();
							KRATOS_WATCH(tmp)
							for(unsigned int kk = 0; kk<TDim; kk++)
							{
								indices.push_back(TDim*tmp + kk);
							}
						}
					}
					std::sort(indices.begin(),indices.end());


					//fill D (the divergence matrix) - note that the "pressure index" is assumed to concide the	DISPLACEMENT_X	index divided by 3
					/*					
					for(unsigned int j=0; j<indices.size(); j++)
					{
						//mD.push_back(row_index/TDim, indices[j] , 0.00);
						mD.push_back(row_index, indices[j] , 0.00);
					}
					*/
		
					//clean the indices (it is a work array)
					indices.erase(indices.begin(),indices.end());
				}
			}
			
			KRATOS_WATCH("FSI D")
			KRATOS_WATCH(mD)
			KRATOS_CATCH("")
		
		}
		void BuildAuxiliariesFSI(
			TSystemMatrixType& mD, double density_str,
			ModelPart& r_model_part) //TSystemMatrixType& WorkMatrix, double density_str,	ModelPart& r_model_part)
		{
			KRATOS_TRY			
			KRATOS_WATCH("BUILDING MATRIX DM-1G for the FSI interface nodes")
			typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();

			//mb-right hand side vector (global)
			unsigned int reduced_dim = this->mpb->size() / TDim;
//mb.size()/
//this->mEquationSystemSize / TDim;
				
				//WorkMatrix.resize(reduced_dim, reduced_dim, false);

			TSystemMatrixType WorkMatrix(reduced_dim,reduced_dim);
			KRATOS_WATCH(WorkMatrix)

			boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX;
			array_1d<double,TDim+1> N;
			array_1d<unsigned int ,TDim+1> local_indices;
			//array_1d<double,TDim+1> rhs_contribution;
			double Volume;
			double temp;
			
			TSystemVectorType mMdiagInv(this->mpb->size());
//this->mEquationSystemSize);
			KRATOS_WATCH(mMdiagInv)


			//getting the dof position
			unsigned int dof_position = (r_model_part.NodesBegin())->GetDofPosition(DISPLACEMENT_X);
			//unsigned int dof_position = TDim*((r_model_part.NodesBegin())->GetDofPosition(PRESSURE));

			double aaa = 1.0/(TDim+1.0);
			//if the element is not having all the nodes IS_STRUCTURE, assemble it, otherwise do nothing
			for(ModelPart::ElementsContainerType::iterator i = r_model_part.ElementsBegin(); 
				i!=r_model_part.ElementsEnd(); i++)
			{	

				
				Geometry< Node<3> >& geom = i->GetGeometry();
				//counting the n-r of structure nodes  
				unsigned int str_nr=0;
				
				//for (int k = 0;k<TDim+1;k++)
				for (unsigned int k = 0;k<geom.size();k++)
				{
				str_nr+=(unsigned int)(i->GetGeometry()[k].FastGetSolutionStepValue(IS_STRUCTURE));
				}
				///////////////////////////////////////////////////////////////////////////////////////////////
				//if the element is not having all the nodes IS_STRUCTURE, assemble it, otherwise do nothing
				// that means, that the entries corresponding to the structural elements are zero
				///////////////////////////////////////////////////////////////////////////////////////////////
				if (geom.size()!=str_nr)
				{				
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
							
					//finiding local indices
					//for(int ii = 0; ii<TDim+1; ii++)
					for(unsigned int ii = 0; ii<geom.size(); ii++)
					{
						local_indices[ii] = geom[ii].GetDof(DISPLACEMENT_X,dof_position).EquationId();
					}
					//building matrix D (transpose of the gradient integrated by parts)
							
					//now we check if this is an interface element (i.e. FSI element)

					int n_interface=0;
					for(unsigned int ii = 0; ii<geom.size(); ii++)
					{
						n_interface = geom[ii].FastGetSolutionStepValue(IS_INTERFACE);
					}
					//if the element does not contain interface node, we set the matrix to be zero
					if (n_interface!=0)
						temp = Volume*aaa;
					else 
						temp=0.0;


					for(unsigned int row = 0; row<TDim+1; row++)
					{
						unsigned int row_index = local_indices[row] / (TDim); 
						//it is a structural mass matrix
						mMdiagInv[2*row_index] += (density_str*temp);
						mMdiagInv[(2*row_index)+1] +=(density_str*temp);

						//KRATOS_WATCH(row_index)
						for(unsigned int col = 0; col<TDim+1; col++)
						{											
							for(unsigned int kkk = 0; kkk<TDim; kkk++)
							{
								//check if the below is correct (copied it from Mass matrix)
								unsigned int col_index = local_indices[col]+kkk;
								//unsigned int col_index = col + kkk;
								mD(row_index,col_index) += temp * DN_DX(col,kkk);
							}
						}

					
					}
				}
				
			}
			//inverting the mass matrix
			KRATOS_WATCH(mMdiagInv.size())
			KRATOS_WATCH(mD.size1())
			KRATOS_WATCH(mD.size2())
			for(unsigned int i = 0; i<TSparseSpace::Size(mMdiagInv); i++)
			{
				if (mMdiagInv[i]>1e-26)
					mMdiagInv[i] = 1.0/mMdiagInv[i];
				else{ //if (mMdiagInv[i]==0.0)
					//KRATOS_WATCH(mMdiagInv[i])
					//KRATOS_ERROR(std::logic_error,"something is wrong with the mass matrix entry - ZERO!!!","")					
					mMdiagInv[i] = 1000000000000.0;					
					}
			}
			//finally we will compute the matrix DM-1G (interface Lapalcian)
			//first we multiply G with Minv
			TSystemMatrixType matG;
			
			TSparseSpace::Transpose(mD, matG);
			KRATOS_WATCH(matG)
			for (int i=0;i<mMdiagInv.size();i++)
			{
			for (int j=0;j<matG.size2();j++)
				{
				matG(i,j)*=mMdiagInv[i];
				}
			}

			TSparseSpace::SetToZero(WorkMatrix);
			//TSparseSpace::SetToZero(destination);

			//WorkMatrix = D * G
			//TSparseSpace::Mult(mD, matG, WorkMatrix);
			WorkMatrix=boost::numeric::ublas::prod(mD, matG);
			KRATOS_WATCH(WorkMatrix)

			KRATOS_CATCH (" ")
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
		LapModifiedLinearStrategy(const LapModifiedLinearStrategy& Other);
		
        
        /*@} */   
        
    }; /* Class LapModifiedLinearStrategy */
	
	/*@} */
	
	/**@name Type Definitions */       
	/*@{ */
	
	
	/*@} */
	
}  /* namespace Kratos.*/
#endif /* KRATOS_LAP_MODIFIED_LINEAR_STRATEGY  defined */

