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
*   Last Modified by:    $Author: Nelson $
*   Date:                $Date: 2009-02-5 
*   Revision:            $Revision: 1.01 $
*
* ***********************************************************/


#if !defined(KRATOS_RESIDUALBASED_ARC_LENGHT_STRATEGY )
#define  KRATOS_RESIDUALBASED_ARC_LENGHT_STRATEGY


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"



/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "custom_utilities/line_searches_utility.h"

#include <cmath>

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
	class TDenseSpace, // = DenseSpace<double>,
	class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
	>
	class ResidualBasedArcLengthStrategy 
		: public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
		 
		  
	{
	public:
		/**@name Type Definitions */       
		/*@{ */
		typedef ConvergenceCriteria<TSparseSpace,TDenseSpace> TConvergenceCriteriaType;

		/** Counted pointer of ClassName */
		//typedef boost::shared_ptr< ResidualBasedArcLenghtStrategy<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;
		KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedArcLengthStrategy );

		typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

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
		ResidualBasedArcLengthStrategy(
			ModelPart& model_part, 
			typename TSchemeType::Pointer pScheme,
			typename TLinearSolver::Pointer pNewLinearSolver,
			typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
			unsigned  int Ide,
			unsigned  int MaxIterations,
			double factor_delta_lmax,
			bool CalculateReactions     = true,
			bool ReformDofSetAtEachStep = true,
			bool MoveMeshFlag           = true,	
			bool ApplyBodyForce         = true // if it is true step should begin in 0
			)
			: SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag)

		 {
			KRATOS_TRY
			//set flags to default values
			SetMaxIterationNumber(MaxIterations);
			mCalculateReactionsFlag = CalculateReactions;
			//KRATOS_WATCH(Parameters_Line_Searches)
			//mstep = step;
			// creating models part for analysis
			mAuxElementModelPart.SetBufferSize(model_part.GetBufferSize());
			mAuxConditionModelPart.SetBufferSize(model_part.GetBufferSize());
			InitializeAuxiliaryModelParts(model_part);
			
		        mfactor_delta_lmax      = factor_delta_lmax;                            		
		        mIde	                = Ide;
			mApplyBodyForce         = ApplyBodyForce;
                        
                        
			
			mReformDofSetAtEachStep = ReformDofSetAtEachStep;

			//saving the convergence criteria to be used 
			mpConvergenceCriteria = pNewConvergenceCriteria;

			//saving the scheme
			mpScheme = pScheme;

			//saving the linear solver
			mpLinearSolver = pNewLinearSolver;

			//setting up the default builder and solver
			mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer
				(
				new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(mpLinearSolver)
				);

			//set flags to start correcty the calculations
			mSolutionStepIsInitialized = false;
			mInitializeWasPerformed = false;

			//tells to the builder and solver if the reactions have to be Calculated or not
			GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

			//tells to the Builder And Solver if the system matrix and vectors need to
			//be reshaped at each step or not
			GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

			//set EchoLevel to the default value (only time is displayed)
			SetEchoLevel(1);

			//by default the matrices are rebuilt at each iteration
			this->SetRebuildLevel(2);

			KRATOS_CATCH("")			
		}
                             

		/** Destructor.
		*/
		virtual ~ResidualBasedArcLengthStrategy() {}

		/** Destructor.
		*/

		//Set and Get Scheme ... containing Builder, Update and other
		void SetScheme(typename TSchemeType::Pointer pScheme ) {mpScheme = pScheme;};
		typename TSchemeType::Pointer GetScheme() {return mpScheme;};

		//Set and Get the BuilderAndSolver
		void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver ) {mpBuilderAndSolver = pNewBuilderAndSolver;};
		typename TBuilderAndSolverType::Pointer GetBuilderAndSolver() {return mpBuilderAndSolver;};

		void SetCalculateReactionsFlag(bool CalculateReactionsFlag) {mCalculateReactionsFlag = CalculateReactionsFlag;}
		bool GetCalculateReactionsFlag() {return mCalculateReactionsFlag;}

		void SetReformDofSetAtEachStepFlag(bool flag) 
                {
					mReformDofSetAtEachStep = flag;
                    GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
                }
		bool GetReformDofSetAtEachStepFlag() {return mReformDofSetAtEachStep;}

		void SetMaxIterationNumber(unsigned int  MaxIterationNumber) {mMaxIterationNumber = MaxIterationNumber;}
		unsigned int GetMaxIterationNumber() {return mMaxIterationNumber;}

		//level of echo for the solving strategy
		// 0 -> mute... no echo at all
		// 1 -> printing time and basic informations
		// 2 -> printing linear solver data
		// 3 -> Print of debug informations:
		//		Echo of stiffness matrix, Dx, b...
		void SetEchoLevel(int Level) 
		{
			BaseType::mEchoLevel = Level;
			GetBuilderAndSolver()->SetEchoLevel(Level);		
		}

		void VariablesArcLength()
		{
		  KRATOS_TRY
                  mlamda = BaseType::GetModelPart().GetProcessInfo()[LAMNDA];     
		  mdelta_lamda      = 1.00;
		  meta              = 1.00;
                  KRATOS_CATCH("")	  
  	
		}

		void Solve_Polinomial_Equation( const double& a, const double& b, const double& c, Vector& x_sol)
		{

		    KRATOS_TRY
		    double disc = b*b - 4.00*a*c;
                    //if (fabs(disc)<1E-40) {disc=0.00;}
                    double q = 0.00;
                    x_sol.resize(2,false);
		    if (b > 0) 
		      q = -0.5 * (b + sqrt(disc));
		    else       
		      q = -0.5 * (b - sqrt(disc));
		    
		    x_sol(0) = q / a;
		    x_sol(1) = c / q;

		    KRATOS_CATCH("")
		    
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
				if(mInitializeWasPerformed == false)
					Initialize();
			//initialize solution step
			if (mSolutionStepIsInitialized == false)
				InitializeSolutionStep();
			
			DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();

			TSystemMatrixType& mA = *mpA;
			TSystemVectorType& mDx = *mpDx;
			//TSystemVectorType& mDx_old = *mpDx_old;
			TSystemVectorType& mb = *mpb;
			//TSystemVectorType& mRHS_cond =  *mpRHS_cond;
			//TSystemVectorType& mDelta_p    = *mpDelta_p;
			//TSystemVectorType& mDelta_pold = *mpDelta_pold;



			GetScheme()->Predict(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

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
			
		       TSystemVectorPointerType pSigma_q_cond; // desplazamientos debido a  las cond
		       TSystemVectorPointerType pSigma_q_body; // desplazamientos debido al peso propio
		       TSystemVectorPointerType pSigma_q;      // suma de los dos anteriores
		       TSystemVectorPointerType pSigma_h;      // desplazamiento debido al desequilibrado
		       TSystemVectorPointerType ph;	       // Ortogonal component of h	
		       TSystemVectorPointerType pAux_Vector;
	               TSystemVectorPointerType pe;            // out of balance load
		       //TSystemVectorPointerType pX_old;        // desplazamioento de la ultima iteracion convergida.

		       unsigned int iteration_number=0;	
		       double Ao   = 0.00;
		       double A    = 1.00;
		       double aux  = 0.00;
		       double miu  = 0.00;
                       double g    = 0.00;
                       bool reduce_arc_lenght = true;
		       bool is_converged      = false;
		       bool ResidualIsUpdated = false;
		       
			
		       typename TSchemeType::Pointer pScheme = GetScheme();
		       typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
			
		       std::cout<<"***************************************"<<std::endl;
		       std::cout<<"Solving System With Arc Lenght Strategy "<< std::endl;
		       std::cout<<"***************************************"<<std::endl;

		
			DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();

			// Creating models part for analysis
			InitializeAuxiliaryModelParts(BaseType::GetModelPart());
			mstep = BaseType::GetModelPart().GetProcessInfo()[TIME_STEPS];                
			KRATOS_WATCH(mstep)
                        KRATOS_WATCH(mIde)
                      
			//OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
			//if the operations needed were already performed this does nothing
			//Pointers needed in the solution
			//Initialize The Scheme - OPERATIONS TO BE DONE ONCE
			//Initialize The Elements - OPERATIONS TO BE DONE ONCE
			//Initialisation of the convergence criteria and variables of arc lenght
			if(mInitializeWasPerformed == false)
				Initialize(); 
			//set up the system, operation performed just once unless it is required 
			//to reform the dof set at each iteration
			if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false || 
				mReformDofSetAtEachStep == true )
			{
				//setting up the list of the DOFs to be solved
				pBuilderAndSolver->SetUpDofSet(pScheme,BaseType::GetModelPart());
				
				//shaping correctly the system 				
				pBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());
				
			}

			//prints informations about the current time
			if (this->GetEchoLevel()!=0)
			{
				std::cout << " " << std::endl;
				// model part original
				std::cout << "CurrentTime = " << BaseType::GetModelPart().GetProcessInfo()[TIME] << std::endl;
			}

			//updates the database with a prediction of the solution
			Predict();
                        
			//initialize solution step
			if (mSolutionStepIsInitialized == false)
				InitializeSolutionStep();
			
			// Asignando tamaño
			pBuilderAndSolver->ResizeAndInitializeVectors(mpA,pSigma_q_cond,pSigma_q_body,BaseType::GetModelPart().Elements(),BaseType::GetModelPart().Conditions(),BaseType::GetModelPart().GetProcessInfo());

			pBuilderAndSolver->ResizeAndInitializeVectors(mpA,pSigma_q,pSigma_h,BaseType::GetModelPart().Elements(),BaseType::GetModelPart().Conditions(),BaseType::GetModelPart().GetProcessInfo());

			pBuilderAndSolver->ResizeAndInitializeVectors(mpA,mpX_old,ph,BaseType::GetModelPart().Elements(),BaseType::GetModelPart().Conditions(),BaseType::GetModelPart().GetProcessInfo());

			pBuilderAndSolver->ResizeAndInitializeVectors(mpA,pe,pAux_Vector,BaseType::GetModelPart().Elements(),BaseType::GetModelPart().Conditions(),BaseType::GetModelPart().GetProcessInfo());


			TSystemMatrixType& mA            = *mpA;
			TSystemVectorType& mDx           = *mpDx;
			TSystemVectorType& mb            = *mpb;
			TSystemVectorType& mRHS_cond     = *mpRHS_cond;
			TSystemVectorType& Sigma_q_cond  = *pSigma_q_cond;
		        TSystemVectorType& Sigma_q_body  = *pSigma_q_body;
			TSystemVectorType& mDelta_p      = *mpDelta_p;
			TSystemVectorType& mDelta_pold   = *mpDelta_pold;
			TSystemVectorType& Sigma_q       = *pSigma_q;
			TSystemVectorType& h             = *ph;
			TSystemVectorType& Sigma_h       = *pSigma_h;
			TSystemVectorType& Aux_Vector    = *pAux_Vector;
			TSystemVectorType& e             = *pe;
			TSystemVectorType& mX_old        = *mpX_old;

			
			
			//function to perform the building and the solving phase.
			if(BaseType::mRebuildLevel >1 || BaseType::mStiffnessMatrixIsBuilt == false)
			{
				TSparseSpace::SetToZero(mA);
				TSparseSpace::SetToZero(mDx);
				TSparseSpace::SetToZero(mb); // fuera interna
				TSparseSpace::SetToZero(mRHS_cond);
				TSparseSpace::SetToZero(Sigma_q_cond);   
				TSparseSpace::SetToZero(Sigma_q_body);
				TSparseSpace::SetToZero(Sigma_q);
				TSparseSpace::SetToZero(mDelta_p); 
				TSparseSpace::SetToZero(h); 
				TSparseSpace::SetToZero(Sigma_h);
				TSparseSpace::SetToZero(Aux_Vector);
				TSparseSpace::SetToZero(e);
				TSparseSpace::SetToZero(mX_old);
	    
				// mb        = Fuerza Interna
				// RHS_cond  = Fuerza  Externa
				pBuilderAndSolver->Build(pScheme,mAuxElementModelPart,mA,mb); 
				pBuilderAndSolver->BuildRHS(pScheme,mAuxConditionModelPart,mRHS_cond);	
			}
			else
			{
				TSparseSpace::SetToZero(mA);
				TSparseSpace::SetToZero(mDx); //mDx=0.00;
				TSparseSpace::SetToZero(mb); // Fuerza interna
				TSparseSpace::SetToZero(mRHS_cond);
				TSparseSpace::SetToZero(Sigma_q_cond);   
				TSparseSpace::SetToZero(Sigma_q_body); 
				TSparseSpace::SetToZero(mDelta_p); 
				TSparseSpace::SetToZero(h); 
				TSparseSpace::SetToZero(Sigma_h);
				TSparseSpace::SetToZero(Aux_Vector);
				TSparseSpace::SetToZero(e);
				TSparseSpace::SetToZero(Sigma_q);
				TSparseSpace::SetToZero(mX_old); 

				
				pBuilderAndSolver->Build(pScheme,mAuxElementModelPart,mA,mb); 
				pBuilderAndSolver->BuildRHS(pScheme,mAuxConditionModelPart,mRHS_cond);
			}
			
		    
			// Guardando los desplazamientos antiguos en X_old. En 1era iteracion X_old = 0
			this->BackupDatabase(rDofSet,mX_old); 
 
			// Siempre resolveremos primero bajo peso propio.
                        // en caso de no considrerar peso propio colocar en body_force un valor cercano a cero.(1E-10)
			if (mstep==0 && mApplyBodyForce==true) 
			{
			    pBuilderAndSolver->SystemSolve(mA,Sigma_q_body,mb);
			    TSparseSpace::Copy(Sigma_q_body,mDelta_p); 
			}
                        else

			{
			    pBuilderAndSolver->SystemSolve(mA,Sigma_q_cond,mRHS_cond);
			    TSparseSpace::Copy(Sigma_q_cond,Sigma_q);
			}			
			//Iteration Cicle... performed only for NonLinearProblems
		        std::cout<<"************************************************************************"<<std::endl;
		        std::cout<<"Begininning Arc Lenght Method.A Pseudo-Line Searches Included. Please Wait...."<<std::endl;
			std::cout<<"************************************************************************"<<std::endl;
			
			while(reduce_arc_lenght==true)
                        {
			while(	is_converged == false &&
				iteration_number++<mMaxIterationNumber) 
			{
				//setting the number of iteration
                                std::cout<<std::endl;				
				KRATOS_WATCH(iteration_number);
				KRATOS_WATCH(mstep);
				KRATOS_WATCH(mMaxIterationNumber);
				is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,mb); // retorna true
			        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;	
			        pScheme->InitializeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);
				
			   if (iteration_number==1)
				{
				   
				  if (mstep==0 && mApplyBodyForce==true)
				      {
					 std::cout<<"Resolving Structure For Dead Load "<< std::endl;
					 mlamda_old       = 0.00; 

				      }

				   else if (mstep==1)
					{  
					    mlamda_old       = 0.00;  
					    Ao               = TSparseSpace::Dot(Sigma_q,Sigma_q);     //Ao = inner_prod (Sigma_q, Sigma_q);
					    mdelta_l         = sqrt(2*Ao*mdelta_lamda*mdelta_lamda);
                                            mdelta_lold      = mdelta_l;
                                            mdelta_lmax      = mdelta_l*mfactor_delta_lmax;              //  Tamaño maximo de arc-length
                                            TSparseSpace::Assign(mDelta_p, mdelta_lamda,Sigma_q);       // mDelta_p = mdelta_lamda*Sigma_q; 
					    mdelta_lamda_old  = mdelta_lamda;
                                            TSparseSpace::Copy(mDelta_p ,mDelta_pold);                  // mDelta_pold      = mDelta_p;
					    KRATOS_WATCH(mdelta_l)
                                            KRATOS_WATCH(mdelta_lamda)
					    KRATOS_WATCH(Ao)
                             
					}
				    else
					{
					    
					    
					    KRATOS_WATCH(mlamda_old)
					    KRATOS_WATCH(mdelta_lamda_old)
					      
					    aux          = TSparseSpace::Dot(mDelta_pold,mDelta_pold);   //inner_prod(mDelta_pold,mDelta_pold);
                                            Ao           = aux/(mlamda_old*mlamda_old); 
					    miu          = mdelta_l/sqrt(aux + Ao*mdelta_lamda_old*mdelta_lamda_old);
					    TSparseSpace::Assign(mDelta_p, miu, mDelta_pold);            //mDelta_p     = miu*mDelta_pold;
				            mdelta_lamda = miu*mdelta_lamda_old;

					    KRATOS_WATCH(mdelta_l)
					    KRATOS_WATCH(Ao)
					    KRATOS_WATCH(miu)
					    KRATOS_WATCH(mdelta_lamda)
					    
                                           
					}
				  }
			    else
			        {
				    if (mstep==0)
				      {
					pBuilderAndSolver->SystemSolve(mA,Sigma_q_body,mb);
					TSparseSpace::Copy(Sigma_q_body,mDelta_p); 
				      }
				    else
				    { 
				    TSystemVectorType& Sigma_h = *pSigma_h;
				    TSparseSpace::SetToZero(Sigma_h);
				    //TSparseSpace::SetToZero(mRHS_cond);
				    TSparseSpace::SetToZero(mDx);
				    TSparseSpace::SetToZero(mDelta_p);
				    TSparseSpace::SetToZero(h);
				    
				    //pBuilderAndSolver->BuildRHS(pScheme,mAuxConditionModelPart,mRHS_cond);
				    KRATOS_WATCH(mlamda)
                                    KRATOS_WATCH(meta)
                                    //KRATOS_WATCH(mRHS_cond)
                                    TSparseSpace::ScaleAndAdd(mlamda,mRHS_cond,A,mb,e);    //e = mlamda*mRHS_cond + mb; // desequilibrado=> e ver ecuacion (7)
				    g = TSparseSpace::Dot(e,mRHS_cond)/TSparseSpace::Dot(mRHS_cond,mRHS_cond);
				    TSparseSpace::ScaleAndAdd(A,e,-g,mRHS_cond,h);            //noalias(h) = (e - g*mRHS_cond); 
				    
					
				    // precaucion si los componenetes de h = 0.00
                                    if (norm_2(h) < 1.0E-10)
					{
					    TSparseSpace::Copy(e,h);
                                            g = 0.00;
					} 
				    
				    pBuilderAndSolver->SystemSolve(mA,Sigma_h,h);
				    				    
				    //Funcion Recursiva
				    Recursive_Function_Arc_Length(pSigma_q, pAux_Vector,pSigma_h, g, Ao);
				    }
				    
				    }


			 // Computing Convergence
                        TSparseSpace::SetToZero(mDx); //mDx=0.00;
                        TSparseSpace::SetToZero(mb); //mDx=0.00; 
			mlamda  = mlamda_old + mdelta_lamda;
                        TSparseSpace::Copy(mDelta_p,mDx); //mDx= mDelta_p;
		        KRATOS_WATCH(mlamda_old);
			//KRATOS_WATCH(mb);
			//KRATOS_WATCH(mDelta_p);

                         // hay que resolver simultaneamente las ecuaciones para lamda y meta.
			 // leer articulo An Al icluding line searches de crisfield, para proxima version 1.02
			 // las contantes para resolver sistema de ecuacion cambia, ya que se multiplican por meta.

			 /*
                         if ( this->mApplyLineSearches==true)
				    {
					   
					  Satisfactory_Line_Search = this->LineSearches(BaseType::GetModelPart(),
											pScheme,
											pBuilderAndSolver,
											rDofSet,
											X_old,mDelta_p,mDx,mb,mA);
					  if ( Satisfactory_Line_Search== true)		
					      {
						    std::cout<<"***************************************************"<<std::endl; 
						    std::cout<<"******Line Searches Has Succecfuly Finished*******"<<std::endl; 
						    std::cout<<"***************************************************"<<std::endl;
					      }
					   
				      
					      
					      rDofSet = pBuilderAndSolver->GetDofSet();
					      this->SetDatabaseToValue(rDofSet, X_old);
		 			      TSparseSpace::Assign(mDx,LineSearchesUtility<TSparseSpace,TDenseSpace,TLinearSolver>::meta,mDx);
					      pScheme->Update(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
				

					 }

				 else
				      {
					    //Updating the results stored in the database
					    rDofSet = pBuilderAndSolver->GetDofSet();
					    pScheme->Update(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
				      }

				  */
			//update results
			rDofSet = pBuilderAndSolver->GetDofSet();
			pScheme->Update(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
			//move the mesh if needed
			if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
			

			//KRATOS_WATCH(mA);	   
		
			// No hace nada  
			//pScheme->InitializeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);

			// return true
			//is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
				
			// bucle sobre los elementos y condiciones
			pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);
			
                        //KRATOS_WATCH(mA)
                        //KRATOS_WATCH(mb)

			
			ResidualIsUpdated = false;

				if (is_converged==true) 
				{

					if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
					{
						TSparseSpace::SetToZero(mb); 
						pBuilderAndSolver->BuildRHS(pScheme,BaseType::GetModelPart(),mb);
						ResidualIsUpdated = true;
						//std::cout << "mb is calculated" << std::endl;
					}
					 std::cout<<"Comprobando Convergencia"<<std::endl;
					 is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,h);
				}

                        // recalculando mA y sigma_q
		     
		       if (is_converged == false)
			   {

				TSparseSpace::SetToZero(Sigma_q_cond);   
				TSparseSpace::SetToZero(Sigma_q_body);
				TSparseSpace::SetToZero(Sigma_q);
                                TSparseSpace::SetToZero(mA);
				TSparseSpace::SetToZero(mb);
                                TSparseSpace::SetToZero(mRHS_cond);
                                TSparseSpace::SetToZero(mDx);
				TSparseSpace::SetToZero(mDelta_p);
				TSparseSpace::SetToZero(Aux_Vector);
				TSparseSpace::SetToZero(e);	  
  
				pBuilderAndSolver->Build(pScheme,mAuxElementModelPart,mA,mb); 
				if (mstep != 0)
				{
				pBuilderAndSolver->BuildRHS(pScheme,mAuxConditionModelPart,mRHS_cond);
				pBuilderAndSolver->SystemSolve(mA,Sigma_q_cond,mRHS_cond);
				TSparseSpace::Copy(Sigma_q_cond,Sigma_q);
				}
			   }
 
				      
			if (is_converged == false && iteration_number>=mMaxIterationNumber)
				 {
				    
                                    reduce_arc_lenght = true;
                                    mdelta_l = sqrt(mIde/iteration_number)*mdelta_l;
				    mdelta_lamda = mdelta_lamda_old;
				    MaxIterationsExceeded();
                                    std::cout<<"Longitud de Arco Modificada"<<std::endl;
                                    KRATOS_WATCH(mdelta_l);
				 }



                        // Si hemos alcanzado convergencia

			if (is_converged==true) 
			{	  
			//recalculate residual if needed 
			// (note that some convergence criteria need it to be recalculated)
			if (ResidualIsUpdated==false)
			{
				TSparseSpace::SetToZero(mb);

				pBuilderAndSolver->BuildRHS(pScheme,BaseType::GetModelPart(),mb);

				//std::cout << "mb is calculated" << std::endl;
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
                        this->FinalizeSolutionStep(iteration_number,reduce_arc_lenght);
			BaseType::GetModelPart().GetProcessInfo()[LAMNDA] = mlamda;
                        KRATOS_WATCH(mdelta_l)

 
			//Cleaning memory after the solution
			pScheme->Clean();

			
			mSolutionStepIsInitialized = false;

			if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
			{
				//TSparseSpace::Clear(mA);
				//TSparseSpace::Clear(mDx);
				//TSparseSpace::Clear(mb);
				
				SparseSpaceType::Clear(mpA);
				SparseSpaceType::Clear(mpA);
				SparseSpaceType::Clear(mpDx);
				SparseSpaceType::Clear(mpb);
				SparseSpaceType::Clear(mpDelta_p);
			}

                        }
                        }
                        }

			return 0.00;

			KRATOS_CATCH("")

		}

		/** 
		this should be considered as a "post solution" convergence check which is useful for coupled analysis
		- the convergence criteria used is the one used inside the "solve" step
		*/
		//**********************************************************************
		bool IsConverged()
		{
			KRATOS_TRY

			TSystemMatrixType& mA = *mpA;
			TSystemVectorType& mDx = *mpDx;
			TSystemVectorType& mb = *mpb;


			if (mpConvergenceCriteria->mActualizeRHSIsNeeded == true)
			{
				GetBuilderAndSolver()->BuildRHS(GetScheme(),BaseType::GetModelPart(),mb);
			}

			DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();

			return mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
			KRATOS_CATCH("")

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

		//**********************************************************************
		//**********************************************************************
		void Clear()
		{
			KRATOS_TRY
			std::cout << "Arc Length Strategy  Clear function used" << std::endl;

			TSystemMatrixType& mA   = *mpA;
			TSystemVectorType& mDx  = *mpDx;
			TSystemVectorType& mb   = *mpb;
			TSystemVectorType& mRHS_cond = *mpRHS_cond;
			TSystemVectorType& mDelta_p = *mpDelta_p;

			SparseSpaceType::Clear(mpA);
			SparseSpaceType::Resize(mA,0,0);

			SparseSpaceType::Clear(mpDx);
			SparseSpaceType::Resize(mDx,0);

			SparseSpaceType::Clear(mpb);
			SparseSpaceType::Resize(mb,0);
	  
			SparseSpaceType::Clear(mpRHS_cond);
			SparseSpaceType::Resize(mRHS_cond,0);

			SparseSpaceType::Clear(mpDelta_p);
			SparseSpaceType::Resize(mDelta_p,0);	


			//setting to zero the internal flag to ensure that the dof sets are recalculated
			GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
			GetBuilderAndSolver()->Clear();
				
				
			KRATOS_CATCH("");
		}

		//**********************************************************************
		//**********************************************************************

		 void Recursive_Function_Arc_Length(TSystemVectorPointerType& pSigma_q,TSystemVectorPointerType& pAux_Vector,TSystemVectorPointerType& pSigma_h, double& g, double& Ao)
		      {
			KRATOS_TRY
			  
		         // coeficients to resolve 2nd equation.
                         double a   = 0.00;
		         double b   = 0.00;
                         double c   = 0.00;
		         double disc = 0.00;           // discriminante de la ecuacion cuadratica.
		         double x    = 0.00;
		         Vector x_sol(2);              // Solucion de la ecuacion de segundo grado
                         x_sol.resize(2,false);

			 double lamda_cr       = 0.00;
			 double delta_lamda_cr = 0.00;
			 double delta_lcr      = 0.00;
                         double miu            = 0.00;
			 bool  imag = false;

                          // Variables vectoriales y matriciles
			TSystemMatrixType& mA            = *mpA;
			TSystemVectorType& mDx           = *mpDx;
			TSystemVectorType& mb            = *mpb;
			TSystemVectorType& mRHS_cond     = *mpRHS_cond;
			TSystemVectorType& mDelta_p      = *mpDelta_p;
			TSystemVectorType& mDelta_pold   = *mpDelta_pold;
			TSystemVectorType& Sigma_q       = *pSigma_q;
			//TSystemVectorType& h             = *ph;
			TSystemVectorType& Sigma_h       = *pSigma_h;
			TSystemVectorType& Aux_Vector    = *pAux_Vector;
			TSystemVectorType& mX_old        = *mpX_old;
			//TSystemVectorType& e             = *pe;

	                 // Constantes necesarias para realizar la operacion de Ublas  
		         double A = 1.00;
                         double B = 1.00;

			 typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
                         DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();
			 rDofSet = pBuilderAndSolver->GetDofSet();
			 typename TSchemeType::Pointer pScheme = GetScheme();
			 
			    
		         // necesario para encontrar raices. mDelta_p se convierte aDx de la iteracion anterior
		         /*for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
			        if(i_dof->IsFree())
				 mDelta_p[i_dof->EquationId()] = i_dof->GetSolutionStepValue(0)-i_dof->GetSolutionStepValue(1);
			  */	    
			  Calculate_Delta_pold(rDofSet,mDelta_p);
                          a = Ao + TSparseSpace::Dot(Sigma_q,Sigma_q);     
                          TSparseSpace::ScaleAndAdd(A,mDelta_p,meta,Sigma_h,Aux_Vector);   
		          b = 2.00*(Ao*(mdelta_lamda-g) + TSparseSpace::Dot(Sigma_q,Aux_Vector)); 
		          c = Ao*(mdelta_lamda-g)*(mdelta_lamda-g) - mdelta_l*mdelta_l + TSparseSpace::Dot(Aux_Vector,Aux_Vector); 
				    
			  KRATOS_WATCH(a);
		          KRATOS_WATCH(b);
		          KRATOS_WATCH(c);

			  TSystemVectorPointerType pDelta_p1;
		          TSystemVectorPointerType pDelta_p2;
				   
		          // asignando tamaño a delta_p1 y delta_p2   
			  pBuilderAndSolver->ResizeAndInitializeVectors(mpA,pDelta_p1,pDelta_p2,BaseType::GetModelPart().Elements(),BaseType::GetModelPart().Conditions(),BaseType::GetModelPart().GetProcessInfo());
				   
			  TSystemVectorType& Delta_p1  = *pDelta_p1;
		          TSystemVectorType& Delta_p2  = *pDelta_p2;
				    
		          // setting to zero elements 
		          TSparseSpace::SetToZero(Delta_p1);
		          TSparseSpace::SetToZero(Delta_p2);
                          TSparseSpace::SetToZero(mDelta_p);
                                   
			  disc = b*b - 4.00*a*c;
		          KRATOS_WATCH(disc);

		          if (disc>=0.00)
			      {
			         Solve_Polinomial_Equation(a,b,c,x_sol);
                                 TSparseSpace::ScaleAndAdd(x_sol(0),Sigma_q,meta,Sigma_h,Delta_p1); //Delta_p1 = x_sol(0)*Sigma_q + meta*Sigma_h
			         TSparseSpace::ScaleAndAdd(x_sol(1),Sigma_q,meta,Sigma_h,Delta_p2); //Delta_p2 = x_sol(1)*Sigma_q + meta*Sigma_h
			         
				 KRATOS_WATCH(x_sol)		 
				 x = x_sol(1);
                                 TSparseSpace::Copy(Delta_p2,mDelta_p); //mDelta_p = Delta_p2;
                                 

                                 // Eligiendo valor de x
			          if (TSparseSpace::Dot(Delta_p1,mDelta_pold) > TSparseSpace::Dot(Delta_p2,mDelta_pold))
				    {
					 x= x_sol(0);
					 TSparseSpace::Copy(Delta_p1,mDelta_p); //mDelta_p = Delta_p1;
					 
				     }
					  
				  mdelta_lamda= mdelta_lamda - g + x;
                                  KRATOS_WATCH(mdelta_lamda)
				  KRATOS_WATCH(x)  
				    
				}
					
					  
			  else
				{	
					
				   std::cout<<"Discriminante Negativo"<<std::endl;
			           std::cout<<"Calculating eta"<<std::endl;
				   Calculate_eta(Ao,pSigma_q,pSigma_h,g,imag);
				   KRATOS_WATCH(meta)
				   KRATOS_WATCH(imag)
				   KRATOS_WATCH(g)
				   KRATOS_WATCH(Ao)

				   if (imag==false && disc >= 0.00)
					{
					   Recursive_Function_Arc_Length(pSigma_q,pAux_Vector,pSigma_h,g,Ao);
					}

				    else
					{   /*	
					    for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
						if(i_dof->IsFree())
						mDelta_p[i_dof->EquationId()] = i_dof->GetSolutionStepValue(0)-i_dof->GetSolutionStepValue(1);
					     */

					    TSparseSpace::SetToZero(mX_old);
					    this->BackupDatabase(rDofSet,mX_old); 
					    Calculate_Delta_pold(rDofSet,mDelta_p);
					    TSparseSpace::ScaleAndAdd(A,mDelta_p,B,Sigma_h,mDelta_p);  //mDelta_p = mDelta_p + Sigma_h; // delta p_cr
					    TSparseSpace::Copy(mDelta_p,mDx);    //mDx      = mDelta_p;
                                                
                                             //update results
					    rDofSet = pBuilderAndSolver->GetDofSet();
					    pScheme->Update(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
					    if(this->MoveMeshFlag() == true) BaseType::MoveMesh();	
					   
					    TSparseSpace::SetToZero(mb);
					    TSparseSpace::SetToZero(mRHS_cond);
					    pBuilderAndSolver->BuildRHS(pScheme,mAuxElementModelPart,mb);
					    pBuilderAndSolver->BuildRHS(pScheme,mAuxConditionModelPart,mRHS_cond);
					    lamda_cr = TSparseSpace::Dot(mb,mRHS_cond)/TSparseSpace::Dot(mRHS_cond,mRHS_cond);
	                                    delta_lamda_cr = lamda_cr - mlamda_old;
                                            delta_lcr = sqrt(inner_prod(mDelta_p,mDelta_p) + Ao*(delta_lamda_cr)*(delta_lamda_cr));
 				            miu = mdelta_l/delta_lcr;
					    TSparseSpace::ScaleAndAdd((miu-1.00),mDelta_p,miu,Sigma_h,mDelta_p); 
					    //TSparseSpace::Assign(mDelta_p,(miu-1),mDelta_p);
                                            //mDelta_p = (1.00-miu)*mDelta_p;
                                            mdelta_lamda = miu*delta_lamda_cr;
                                            KRATOS_WATCH(lamda_cr)
				            KRATOS_WATCH(delta_lamda_cr)
				            KRATOS_WATCH(miu)
				            KRATOS_WATCH(mdelta_lamda)
				  
					    // dejandolo como estaba antes
					    // mDx = -mDx;
					     this->SetDatabaseToValue(rDofSet, mX_old);
					
					    //update results
					    //rDofSet = pBuilderAndSolver->GetDofSet();
				            //pScheme->Update(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
					    //pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);

					    TSparseSpace::SetToZero(mb);
					    TSparseSpace::SetToZero(mDx);
					    TSparseSpace::SetToZero(mA);
                                                		        
					}
				           
				 } 

			KRATOS_CATCH("");
		      }

		//**********************************************************************
		//**********************************************************************
		
		TSystemMatrixType& GetSystemMatrix()
		    {
			TSystemMatrixType& mA = *mpA;

			return mA;
		    }
	

	protected:

	private:

		typename TSchemeType::Pointer mpScheme;     

		typename TLinearSolver::Pointer mpLinearSolver;

		typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;

		typename TConvergenceCriteriaType::Pointer mpConvergenceCriteria;

		TSystemVectorPointerType mpDx;
		TSystemVectorPointerType mpb;
		TSystemMatrixPointerType mpA;
		TSystemVectorPointerType mpRHS_cond;
		TSystemVectorPointerType mpX_old;  // Solucion en paso anterior convergido.
		TSystemVectorPointerType mpDelta_p;
		TSystemVectorPointerType mpDelta_pold;
		
  		/** 
		Flag telling if it is needed to reform the DofSet at each
		solution step or if it is possible to form it just once
		- true  => reforme at each time step
		- false => form just one (more efficient)

		Default = false
		*/
		bool mReformDofSetAtEachStep;	

		/** 
		Flag telling if it is needed or not to compute the reactions

		default = true
		*/
		bool mCalculateReactionsFlag;
		bool mInitializeWasPerformed;
                bool mApplyBodyForce; 

		bool mSolutionStepIsInitialized;
		unsigned int mMaxIterationNumber;
		unsigned int mstep; 
		

		// old = incremento anterior convergido

		double mdelta_l;         // Longitud de Arco
                double mdelta_lold;      // Longitud de Arco del incremeto anterior convergido.
		double mdelta_lmax;      // Longitud de Arco Maxima Permitida
                double mfactor_delta_lmax;
                double meta;			
                double mIde;             
		double mlamda;           
		double mlamda_old;       
                double mdelta_lamda;     
		double mdelta_lamda_old;
                ModelPart mAuxElementModelPart;
		ModelPart mAuxConditionModelPart; 
		

		/*@} */
		/**@name Private Operators*/
		/*@{ */
		//**********************************************************************
		//**********************************************************************
		void Initialize()
		{
			KRATOS_TRY

			//pointers needed in the solution
			typename TSchemeType::Pointer pScheme = GetScheme();
			typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;

			//Initialize The Scheme - OPERATIONS TO BE DONE ONCE
			if (pScheme->SchemeIsInitialized() == false) 
				pScheme->Initialize(BaseType::GetModelPart());

			//Initialize The Elements - OPERATIONS TO BE DONE ONCE
			if (pScheme->ElementsAreInitialized() == false) 
				pScheme->InitializeElements(BaseType::GetModelPart());

			//initialisation of the convergence criteria
			if (mpConvergenceCriteria->mConvergenceCriteriaIsInitialized == false)
				mpConvergenceCriteria->Initialize(BaseType::GetModelPart());


			mInitializeWasPerformed = true;
			VariablesArcLength(); // iniciando variables de arc lenght
		    
			KRATOS_CATCH("")
		}


		//**********************************************************************
		//**********************************************************************
		void InitializeSolutionStep()
		{
			KRATOS_TRY

			typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
			typename TSchemeType::Pointer pScheme = GetScheme();	

			//setting up the Vectors involved to the correct size with value cero
			pBuilderAndSolver->ResizeAndInitializeVectors(mpA,mpDx,mpb,BaseType::GetModelPart().Elements(),BaseType::GetModelPart().Conditions(),BaseType::GetModelPart().GetProcessInfo());
		    
			 pBuilderAndSolver->ResizeAndInitializeVectors(mpA,mpX_old,mpRHS_cond,BaseType::GetModelPart().Elements(),BaseType::GetModelPart().Conditions(),BaseType::GetModelPart().GetProcessInfo()); 


			pBuilderAndSolver->ResizeAndInitializeVectors(mpA,mpDelta_p,mpDelta_pold,BaseType::GetModelPart().Elements(),BaseType::GetModelPart().Conditions(),BaseType::GetModelPart().GetProcessInfo()); 
			
			
			TSystemMatrixType& mA          = *mpA;
			TSystemVectorType& mDx         = *mpDx;
			//TSystemVectorType& mDx_old     = *mpDx_old;
			TSystemVectorType& mb          = *mpb;
			//TSystemVectorType& mRHS_cond   = *mpRHS_cond;
			//TSystemVectorType& mDelta_p    = *mpDelta_p;
			//TSystemVectorType& mDelta_pold = *mpDelta_pold;
			

			//initial operations ... things that are constant over the Solution Step
			pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
			//pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),mA,mDx,mRHS_cond);	
			
			//initial operations ... things that are constant over the Solution Step
			//pScheme->InitializeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
		
			pScheme->InitializeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
			
			mSolutionStepIsInitialized = true;	

			KRATOS_CATCH("")
		}	



	    void FinalizeSolutionStep( unsigned int& iteration_number, bool &reduce_arc_lenght)
		{
		  KRATOS_TRY

		      typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();

		      pBuilderAndSolver->ResizeAndInitializeVectors(mpA,mpDelta_p,mpDelta_pold,BaseType::GetModelPart().Elements(),BaseType::GetModelPart().Conditions(),BaseType::GetModelPart().GetProcessInfo());
			
		      DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();

	              TSystemVectorType& mDelta_pold = *mpDelta_pold;
		      TSystemVectorType& mX_old = *mpX_old;
		      double factor = 0.00;			
		      mdelta_lamda_old = mdelta_lamda;	
                      //mdelta_l = mdelta_lold; 
		      factor = sqrt(mIde/iteration_number);   
		      if (factor > 1.5)
			  {
			    factor = 1.50;
			  } 
		      if (factor < 0.25)  
			  {
                            factor = 0.25;
			  } 
		      mdelta_l = factor*mdelta_lold;                                            
		      mlamda_old = mlamda;
                      //KRATOS_WATCH(iteration_number)
		      TSparseSpace::SetToZero(mDelta_pold);
		      /*
		      for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
					if(i_dof->IsFree())
						mDelta_pold[i_dof->EquationId()] = i_dof->GetSolutionStepValue(0) - i_dof->GetSolutionStepValue(1);
		      */
			
		     Calculate_Delta_pold(rDofSet,mDelta_pold);
		     this->BackupDatabase(rDofSet,mX_old);  

		      //KRATOS_WATCH(mDelta_pold);
		      // Controlo la longitud del arco.       
		      if (mdelta_lold > mdelta_lmax)
			 {
				mdelta_lold = mdelta_lmax;
				mdelta_l    = mdelta_lmax;
                                KRATOS_WATCH(mdelta_lmax);
			 }     

		    //BaseType::GetModelPart().GetProcessInfo()[ARC_LENGTH_REDUCED] = 0;
                      reduce_arc_lenght= false;  
		    
		     // setting conditions
		    KRATOS_CATCH("")
		}	



		//**********************************************************************
		//**********************************************************************
		void MaxIterationsExceeded()
		{
			std::cout << "***************************************************" << std::endl;
			std::cout << "******* ATTENTION: Max Iterations Exceeded ********" << std::endl;
			std::cout << "***************************************************" << std::endl;
		}

		//**********************************************************************
		//**********************************************************************

		void InitializeAuxiliaryModelParts(ModelPart& ThisModelPart)
		{
		    //mAuxElementModelPart.SetBufferSize(ThisModelPart.GetBufferSize());
		    //mAuxConditionModelPart.SetBufferSize(ThisModelPart.GetBufferSize());

		    mAuxElementModelPart.Nodes() = ThisModelPart.Nodes();
		    mAuxConditionModelPart.Nodes() = ThisModelPart.Nodes();

/* 		    mAuxElementModelPart.Properties() = ThisModelPart.Properties(); */
/* 		    mAuxConditionModelPart.Properties() = ThisModelPart.Properties(); */

		    mAuxElementModelPart.Elements() = ThisModelPart.Elements();
		    mAuxConditionModelPart.Conditions() = ThisModelPart.Conditions();

		}

		//**********************************************************************
		//**********************************************************************

		    void Calculate_eta(double& Ao,TSystemVectorPointerType& pSigma_q, TSystemVectorPointerType& pSigma_h,double& g, bool& imag)
			{
			    double a_prima= 0.00; 
			    double b_prima= 0.00;
			    double c_prima= 0.00;
			    double disc   = 0.00;
                            Vector SOL;
			    SOL.resize(2,false);


			    TSystemVectorType& Sigma_q     = *pSigma_q;
			    TSystemVectorType& Sigma_h     = *pSigma_h;
			    TSystemVectorType& mDelta_p    = *mpDelta_p;
			    //KRATOS_WATCH(Sigma_q)
                            //KRATOS_WATCH(Sigma_h)
                             
                           
			    typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
                            DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();
			    rDofSet = pBuilderAndSolver->GetDofSet();

			    // necesario para encontrar raices. mDelta_p se convierte aDx de la iteracion anterior
			    /* for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
				if(i_dof->IsFree())
				    mDelta_p[i_dof->EquationId()] = i_dof->GetSolutionStepValue(0)-i_dof->GetSolutionStepValue(1);
			    */
			    //KRATOS_WATCH(mDelta_p)
	  
			    Calculate_Delta_pold(rDofSet,mDelta_p);  
			    a_prima = (Ao + inner_prod(Sigma_q,Sigma_q))*(inner_prod(Sigma_h,Sigma_h))-(inner_prod(Sigma_q,Sigma_h))*(inner_prod(Sigma_q,Sigma_h)); 

			    b_prima = 2.00*((Ao + inner_prod(Sigma_q,Sigma_q))*(inner_prod(mDelta_p,Sigma_h))-((Ao*(mdelta_lamda-g)+ inner_prod(Sigma_q,mDelta_p)))*(inner_prod(Sigma_q,Sigma_h)));

			    c_prima = (Ao + inner_prod(Sigma_q,Sigma_q))*((inner_prod(mDelta_p,mDelta_p)-mdelta_l*mdelta_l))-(2.00*Ao*(mdelta_lamda-g)+inner_prod(Sigma_q,mDelta_p))*inner_prod(Sigma_q,mDelta_p) + inner_prod(Sigma_q,Sigma_q)*Ao*(mdelta_lamda-g)*(mdelta_lamda-g);

			    disc    = b_prima*b_prima*-4.00*a_prima*c_prima;
			    KRATOS_WATCH(a_prima)
		            KRATOS_WATCH(b_prima)
		            KRATOS_WATCH(c_prima)

                            KRATOS_WATCH(disc)
			    if (disc >= 0.00)
				 {
				    imag = false;
				    Solve_Polinomial_Equation(a_prima,b_prima,c_prima,SOL);
                                    KRATOS_WATCH(SOL)
				    KRATOS_WATCH(a_prima)
				    KRATOS_WATCH(b_prima)
				    KRATOS_WATCH(c_prima)

				    if (SOL(0)>SOL(1))
				    {
				      double a =   SOL(1);
				      SOL(1)   =   SOL(0);
				      SOL(0)   =    a;
				    }

				    double shi  = 0.05*fabs((SOL(1)-SOL(0)));
				    if (SOL(1) < 1.00)
				      {
					meta  = SOL(1)-shi;
				      }
				    else if ((1 > (-b_prima/a_prima)) && (SOL(1)>1.00))
				      {
					meta  = SOL(1) + shi;
				      }
				    else if ((1 < (-b_prima/a_prima)) && (SOL(0)<1.00))
				      {
					meta  = SOL(0) -shi;
				      }
				    else
				      {
				      meta  = SOL(0) + shi;
				      }
				  }
			    
			    else
				{
				  std::cout<<"eta imaginario"<<std::endl;;
				  meta   = 1.00;
				  imag  = true; 
				 }
		          }


		// Calcula el incremento total producido desde la iteracion i hasta la i+1  
		void Calculate_Delta_pold(
			DofsArrayType const & rDofSet,
			TSystemVectorType& Delta_pold
			) 
			{ 
			KRATOS_TRY
			      for(typename DofsArrayType::const_iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
				  if(i_dof->IsFree())
				      {Delta_pold[i_dof->EquationId()] = i_dof->GetSolutionStepValue(0) - i_dof->GetSolutionStepValue(1);}

			KRATOS_CATCH("")
			 }
		void SetDatabaseToValue(
			DofsArrayType& rDofSet,
			const TSystemVectorType& X_old
			) 
			{ 
			KRATOS_TRY

				for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
				{
					if(i_dof->IsFree())
					{i_dof->GetSolutionStepValue() = X_old[i_dof->EquationId()];}
				}
			KRATOS_CATCH("")
			 }
 
                  
		     // Permite escribir los desplazamientos antiguos en el X_old
		     void BackupDatabase(
			DofsArrayType const& rDofSet,
			TSystemVectorType& X_old
			) 
			{ 
			KRATOS_TRY

				 for(typename DofsArrayType::const_iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
				{
					    if(i_dof->IsFree())
					     { X_old[i_dof->EquationId()] = i_dof->GetSolutionStepValue();}
				}
			KRATOS_CATCH("")
			 }


		     

		ResidualBasedArcLengthStrategy(const ResidualBasedArcLengthStrategy& Other);


		/*@} */   

	}; /* Class ResidualBasedArcLenghtStrategy */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_ARC_LENGTH_STRATEGY  defined */

