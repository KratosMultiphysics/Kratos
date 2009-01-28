/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2007-03-13 15:16:46 $
*   Revision:            $Revision: 1.7 $
*
* ***********************************************************/

#if !defined(KRATOS_RESIDUALBASED_UZAWA_NEWTON_RAPHSON_STRATEGY )
#define  KRATOS_RESIDUALBASED_UZAWA_NEWTON_RAPHSON_STRATEGY

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

#include "custom_conditions/contact_link_3D.h"
#include "structural_application.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

namespace Kratos
{
    /**
     * Uzawa strategy for the solution of structural contact problems.
     * This strategy implements an uzawa algorithm for the solution
     * of contact problems using the augmented lagrange multiplier
     * method.
     */
    template< class TSparseSpace, class TDenseSpace, class TLinearSolver > 
            class ResidualBasedUzawaNewtonRaphsonStrategy 
    : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
    {
        public:
            
            /**
             * Type Definitions 
             */
            
            typedef ConvergenceCriteria<TSparseSpace,TDenseSpace> TConvergenceCriteriaType;
            
            /** 
             * Counted pointer of ResidualBasedUzawaNewtonRaphsonStrategy 
             */
	    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedUzawaNewtonRaphsonStrategy );
//             typedef boost::shared_ptr< ResidualBasedUzawaNewtonRaphsonStrategy< TSparseSpace, 
//             TDenseSpace,TLinearSolver> > Pointer;
            
            typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;
			
			typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> ResidualBasedNewtonRaphsonStrategyType;
            
            typedef std::size_t IndexType;
            
            typedef ModelPart::ConditionsContainerType ConditionsArrayType;
            
            typedef Geometry<Node<3> > GeometryType;
            
            typedef Properties PropertiesType;
            
            typedef typename BaseType::TDataType TDataType;

	    typedef typename BaseType::TSchemeType TSchemeType;

	    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
            
            //typedef typename BaseType::DofSetType DofSetType;
            
            typedef typename BaseType::DofsArrayType DofsArrayType;

            typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

            typedef typename BaseType::TSystemVectorType TSystemVectorType;

            typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

            typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
            
            typedef ResidualBasedEliminationBuilderAndSolver< TSparseSpace, 
            TDenseSpace, TLinearSolver > DefaultBuilderAndSolverType;
            
            
            /**
             * Life Cycle 
             */
            
            /**
             * Constructor.
             */
            ResidualBasedUzawaNewtonRaphsonStrategy( ModelPart& model_part, 
                    typename TSchemeType::Pointer pNewScheme,
                    typename TLinearSolver::Pointer pNewLinearSolver,
                    typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
					int MaxIterations = 30,
					bool CalculateReactions = false,
					bool ReformDofSetAtEachStep = false,
					bool MoveMeshFlag = false
                                                   )
			: SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part,MoveMeshFlag)
            {
                KRATOS_TRY
                //set flags to default values
				SetMaxIterationNumber(MaxIterations);
				mCalculateReactionsFlag = CalculateReactions;
				mReformDofSetAtEachStep = ReformDofSetAtEachStep;

				
				BaseType::GetModelPart().ElementsBegin()->GetProperties()[MAX_UZAWA_ITERATIONS] = MaxIterations;
				
// 				std::cout << "Constructor of ResidualBasedUzawaNewtonRaphsonStrategy" << std::endl;
// 				std::cout << model_part.GetProperties(1).GetValue(MAX_UZAWA_ITERATIONS) << std::endl;
                
                //setting up standard solution strategy for inner
                //iteration loop
                mpNewtonRaphsonStrategy = typename ResidualBasedNewtonRaphsonStrategyType::Pointer(
                        new ResidualBasedNewtonRaphsonStrategy< TSparseSpace,  
                TDenseSpace, TLinearSolver > ( model_part, pNewScheme, pNewLinearSolver,
                                               pNewConvergenceCriteria )  );
                //mesh is always moved!
                mpNewtonRaphsonStrategy->SetMoveMeshFlag( true );
                mpNewtonRaphsonStrategy->SetReformDofSetAtEachStepFlag( true );
                
                //saving the convergence criteria to be used 
                //this criteria is used only for the convergence check of uzawa loop
                  
//                 mpConvergenceCriteria = pNewConvergenceCriteria;
                
                //saving the scheme
//                 mpScheme = pNewScheme;
                
                //saving the linear solver
//                 mpLinearSolver = pNewLinearSolver;
                
                //setting up the default builder and solver
//                 mpBuilderAndSolver = TBuilderAndSolverType::Pointer( new
//                         ResidualBasedEliminationBuilderAndSolver< TSparseSpace, 
//                 TDenseSpace, TLinearSolver >(mpLinearSolver)
//                                                                    );
                
                //set flags to start correcty the calculations
                mSolutionStepIsInitialized = false;
                mInitializeWasPerformed = false;
                
                //tells to the builder and solver if the reactions have to be Calculated or not
//                 GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);
                
                //tells to the Builder And Solver if the system matrix and vectors need to
                //be reshaped at each step or not
//                 GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
                
                //set EchoLevel to the default value (only time is displayed)
                SetEchoLevel(1);
                
                //by default the matrices are rebuilt at each iteration
                this->SetRebuildLevel(2);
                
                KRATOS_CATCH("")
            }
            
            /**
             * Destructor.
             */
            virtual ~ResidualBasedUzawaNewtonRaphsonStrategy() {}
            
			//Set and Get the BuilderAndSolver
			void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver ) {mpNewtonRaphsonStrategy->SetBuilderAndSolver(pNewBuilderAndSolver);};
			typename TBuilderAndSolverType::Pointer GetBuilderAndSolver() {return mpNewtonRaphsonStrategy->GetBuilderAndSolver();};

            
            /**
             * Operations
             */
            
            void SetScheme(typename TSchemeType::Pointer pNewScheme ) {this->mpScheme = pNewScheme;};
            typename TSchemeType::Pointer GetScheme() {return this->mpScheme;};
            
            //Set and Get the BuilderAndSolver
//             void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver ) {mpBuilderAndSolver = pNewBuilderAndSolver;};
//             typename TBuilderAndSolverType::Pointer GetBuilderAndSolver() {return mpBuilderAndSolver;};

            void SetCalculateReactionsFlag(bool CalculateReactionsFlag) {mCalculateReactionsFlag = CalculateReactionsFlag;}
            bool GetCalculateReactionsFlag() {return mCalculateReactionsFlag;}

			void SetReformDofSetAtEachStepFlag(bool flag) {mpNewtonRaphsonStrategy->SetReformDofSetAtEachStepFlag( flag);}
			bool GetReformDofSetAtEachStepFlag() {return mpNewtonRaphsonStrategy->GetReformDofSetAtEachStepFlag();}
			void SetMaxIterationNumber(unsigned int  MaxIterationNumber) {mMaxIterationNumber = MaxIterationNumber;}
			unsigned int GetMaxIterationNumber() {return mMaxIterationNumber;}

//             void SetMaxIterationNumber(unsigned int  MaxIterationNumber) {mMaxIterationNumber = MaxIterationNumber;}
//             unsigned int GetMaxIterationNumber() {return mMaxIterationNumber;}
            
            /**
             * Setting the echo level
             * @param Level the requested echo level:
             * 0 -> mute... no echo at all
             * 1 -> printing time and basic informations
             * 2 -> printing linear solver data
             * 3 -> Print of debug informations:
             */
            void SetEchoLevel(int Level)
            {
                this->mpNewtonRaphsonStrategy->SetEchoLevel( Level );
                this->mEchoLevel = Level;
            }
            
            void SetMoveMeshFlag( bool MoveMeshFlag )
            {
                this->mpNewtonRaphsonStrategy->SetMoveMeshFlag( MoveMeshFlag );
            }
            
            //*********************************************************************************
            /**
             * OPERATIONS ACCESSIBLE FROM THE INPUT:
             */
            
            /**
             * operation to predict the solution.
             * if it is not called a trivial predictor is used in which the
             * values of the solution step of interest are assumed 
             * equal to the old values
             */
            void Predict()
            {
//                 KRATOS_TRY
                //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetition
                //if the operations needed were already performed this does nothing
//                 if(mInitializeWasPerformed == false)
                Initialize();
                
                //initialize solution step
//                 if (mSolutionStepIsInitialized == false)
//                     InitializeSolutionStep();
                
//                 DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();
                
                //ProcessInfo& pCurrentProcessInfo = GetCurrentProcessInfo();
//                 ProcessInfo& pCurrentProcessInfo = GetModelPart().GetProcessInfo();
                
//                 GetScheme()->Predict(GetModelPart(),rDofSet,mA,mDx,mb);
//                 
                //move the mesh if needed
//                 if(MoveMeshFlag() == true) MoveMesh();
                
//                 KRATOS_CATCH("")
            }
            
            
            //*********************************************************************************
            /**
             * the problem of interest is solved
             */
            //**********************************************************************
            double Solve()
            {
// 				std::cout << "UZAWA algorithm started" << std::endl;
// 				std::cout << BaseType::GetModelPart().ElementsBegin()->GetProperties().Data() << std::endl;
								
                KRATOS_TRY
                        
                /**
                 * - storing original condition size before adding virtual
                 *   conditions.
                 * - performing contact search
                 * - creating virtual link conditions for the assembling
                 */
// 				std::cout << "setting up contact conditions: " << std::endl;
				int originalPosition =  SetUpContactConditions();
                    
                bool uzawaConverged = false;
                        
                        
                //KRATOS_WATCH( GetModelPart().GetProcessInfo() );
                //KRATOS_WATCH( GetModelPart() );
				/**
				 * First step: reform DOF set and check if uzawa iteration is necessary
				 */
// 				std::cout << "check if uzawa loop is necessary: " << std::endl;
				this->mpNewtonRaphsonStrategy->Solve();
				Update( originalPosition );
				if( IsConverged( 0,  originalPosition) )//neu
				{
					uzawaConverged = true;
					mpNewtonRaphsonStrategy->SetReformDofSetAtEachStepFlag(true);
					Clean( originalPosition );
					return 0.0;
				}
				/**
				 **********************************************************************
				 * Beginning of UZAWA loop
				 **********************************************************************
				 */
                
				mpNewtonRaphsonStrategy->SetReformDofSetAtEachStepFlag(false);
                for( int uzawaStep=1; 
                     uzawaStep < BaseType::GetModelPart().GetProperties(1)[MAX_UZAWA_ITERATIONS]; 
                     uzawaStep++ )
                {
					
					
                    std::cout << "I am inside the uzawa loop, iteration no. " << uzawaStep << std::endl;
                    
                    /**
                     * Solving the standard newton-raphson iteration 
                     */
                    this->mpNewtonRaphsonStrategy->Solve();
                    
                    /**
                     * - updating the lagrange multipliers
                     */
                    Update( originalPosition );
                    
                    /**
                     * checking convergence!
                     */
//                     if( IsConverged( uzawaStep,  originalPosition) )
                    if( IsConverged( uzawaStep,  originalPosition) )//neu
                    {
                       uzawaConverged = true;
                       break;
                    }

                }
                if( ! uzawaConverged )
                {
                    std::cout << "uzawa algorithm failes to converge within maximum number of iterations" << std::endl;
                }
//                 if( IsConverged( uzawaStep ) )
//                 {
//                     std::cout << "uzawa loop has been terminated successfully!" << std::endl;
//                 }
//                 else
//                 {
//                     std::cout << "uzawa algorithm failes to converge within maximum number of iterations" << std::endl;
//                 }
                /**
                 **********************************************************************
                 * END of UZAWA loop
                 **********************************************************************
                 */ 
                /**
                 * - cleaning up the conditions
                 */
                Clean( originalPosition );
				mpNewtonRaphsonStrategy->SetReformDofSetAtEachStepFlag(true);
                
                KRATOS_CATCH("")
						
				return 0.0;
            }
            
            /** 
             * this should be considered as a "post solution" convergence 
             * check which is useful for coupled analysis
             * - the convergence criteria used is the one used inside the "solve" step
             */
            //**********************************************************************
            bool IsConverged( int step , int lastRealCondition)
            {
                bool Converged = false;
                bool friction = false;
                double ratio = 0;
                double absolute = 0;
                double energy_contact = 0;
                double ratio_friction=0;
                double absolute_friction=0;
                double energy_friction = 0;
//                 double gap = 0;
                //Index is at least 1 to avoid fail in case of no contact conditions
                int Index = 1;
                int Index2 = 1;

                ConditionsArrayType& ConditionsArray = BaseType::GetModelPart().Conditions();
                
                for( typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin();
                     it!=ConditionsArray.ptr_end(); ++it)
                {
                    if( (*it)->GetValue( IS_CONTACT_SLAVE ) )
                    {
                        for( IndexType i=0; i<(*it)->GetGeometry().IntegrationPoints().size();
                             i++ )
                        {
//                             KRATOS_WATCH((*it)->GetValue( LAMBDAS )[i]);
//                             KRATOS_WATCH((*it)->GetValue( LAMBDAS_T )(i,0));
//                             KRATOS_WATCH((*it)->GetValue( LAMBDAS_T )(i,1));
                            if((*it)->GetValue( GAPS )[i]>0)
                            {
                                energy_contact += ((*it)->GetValue( GAPS )[i])
                                        * ((*it)->GetValue( GAPS )[i])
                                        *0.5*((*it)->GetValue(PENALTY)[i]);
                                absolute += (*it)->GetValue( LAMBDAS )[i];
                                ratio += (*it)->GetValue( DELTA_LAMBDAS )[i];
                                Index++;
                            }
                        }
                    }
                }
                
                if( BaseType::GetModelPart().GetProperties(1)[FRICTION_COEFFICIENT] > 0.0 )
                {
                    friction=true;
                    
                    for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin() + lastRealCondition; 
                         it!=ConditionsArray.ptr_end(); ++it)
                    {
                        int i=(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX);
                        if((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(GAPS)[i]  > 0)
                        {
                            Matrix m = (*it)->GetValue(CONTACT_LINK_M);
                            Matrix mInv = ZeroMatrix(2);
                            double det=0.0;
                            MathUtils<double>::InvertMatrix2( m, mInv, det );
                            
                            absolute_friction +=
                                sqrt((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)
                                *mInv(0,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)
                                +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)
                                *mInv(0,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)
                                +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)
                                *mInv(1,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)
                                +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)
                                *mInv(1,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1));
                            
                            ratio_friction += sqrt((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)
                                *mInv(0,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                *mInv(0,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                    +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                *mInv(1,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                    +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                *mInv(1,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(DELTA_LAMBDAS_T )(i,1));       
                           
                            energy_friction += 0.5*
                                    ((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(PENALTY_T)[i])*
                                    ((m(0,0) 
                                * ((*it)->GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT)[0] 
                                - (*it)->GetValue( MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT)[0])
                                    + m(0,1)
                                * ((*it)->GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT)[1]
                                - (*it)->GetValue( MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT)[1]))*
                                    (m(0,0) 
                                * ((*it)->GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT)[0] 
                                - (*it)->GetValue( MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT)[0])
                                    + m(0,1)
                                * ((*it)->GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT)[1]
                                - (*it)->GetValue( MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT)[1]))
                                +
                                    (m(1,0) 
                                * ((*it)->GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT)[0] 
                                - (*it)->GetValue( MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT)[0])
                                    + m(1,1)
                                * ((*it)->GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT)[1]
                                - (*it)->GetValue( MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT)[1]))*
                                    (m(1,0) 
                                * ((*it)->GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT)[0] 
                                - (*it)->GetValue( MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT)[0])
                                    + m(1,1)
                                * ((*it)->GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT)[1]
                                - (*it)->GetValue( MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT)[1])));
                            Index2++;
                        }
                        (*it)->GetValue( MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT)
                                = (*it)->GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT);
                    }
                }
                if( this->mEchoLevel > 2 )
                {
                    std::cout << "absolute Lambda: " << absolute/Index << std::endl;
                    std::cout << "relative Lambda: " << ratio/absolute << std::endl;
                    std::cout << "energy criterion: " << energy_contact/Index << std::endl;
                    std::cout << "K_Contact: " << (this)->GetModelPart().GetProperties(1)[K_CONTACT] << std::endl;
                    if( friction )
					{
						std::cout << "absolute Lambda friction: " << absolute_friction/Index << std::endl;
						std::cout << "relative Lambda friction: " << ratio_friction/absolute_friction << std::endl;
						std::cout << "energy criterion friction: " << energy_friction/Index2 << std::endl;
						std::cout << "K_Contact_Friction: " << (this)->GetModelPart().GetProperties(1)[K_CONTACT_T] << std::endl;
					}
                }
				else if( this->mEchoLevel > 0 )
				{
					std::cout << "relative Lambda: " << ratio/absolute << std::endl;
					if( friction )
					{
						std::cout << "relative Lambda friction: " << ratio_friction/absolute_friction << std::endl;
					}
				}
                
                if( BaseType::GetModelPart().GetProperties(1)[FRICTION_COEFFICIENT] > 0.0 )
                {
                    if( (fabs(absolute/Index) < 1e-3 || fabs(ratio/absolute) < 1e-3)
                         && (fabs(absolute_friction/Index) < 1e-3 || fabs(ratio_friction/absolute_friction) < 1e-2))
                    {
                        if( this->mEchoLevel > 0 )
                            std::cout << "Contact condition converged after " << step+1 << " steps" << std::endl;
                            Converged = true;
                    }
                    else
                    {
                        if( this->mEchoLevel > 0 )
                            std::cout << "Contact has been detected, next solution step required (UZAWA STEP: "<< step+1 << ")" << std::endl;
                    }
                }
                else
                {
                    if( fabs(absolute/Index) < 1e-3 || fabs(ratio/absolute) < 1e-3 )
                     {
                        if( this->mEchoLevel > 0 )
                            std::cout << "Contact condition converged after " << step+1 << " steps" << std::endl;
                            Converged = true;
                     }
                    else
                    {
                        if( this->mEchoLevel > 0 )
                            std::cout << "Contact has been detected, next solution step required (UZAWA STEP: "<< step+1 << ")" << std::endl;
                    }
                }
                return( Converged );
            }
    
            //*********************************************************************************
            /** 
             * this operations should be called before printing the results 
             * when non trivial results (e.g. stresses)
             * need to be calculated given the solution of the step
             * 
             * This operations should be called only when needed, 
             * before printing as it can involve a non negligible cost
             */
            void CalculateOutputData()
            {
                /*
                ProcessInfo& pCurrentProcessInfo = GetModelPart().GetProcessInfo();
        
                DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();
                GetScheme()->CalculateOutputData(GetModelPart(),rDofSet,mA,mDx,mb);		
                */
            }
    
            //*********************************************************************************
    
        protected:
            
            /**
             * there are no protected class members
             */
    
        private:
            
            /**
             * Member Variables 
             */
            
//             typename TSchemeType::Pointer mpScheme;
        
//             typename TLinearSolver::Pointer mpLinearSolver;
        
//             typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;
        
//             typename TConvergenceCriteriaType::Pointer mpConvergenceCriteria;
            
            /**
             * standard Newton-Raphson strategy which is called
             * for the inner iteration loop.
             */
            typename ResidualBasedNewtonRaphsonStrategyType::Pointer mpNewtonRaphsonStrategy;
        
            TSystemVectorType mDx;
            TSystemVectorType mb;
            TSystemMatrixType mA;
        
            TSystemVectorType mOldResidualVector;
            
            
            
            /** 
             * Flag telling if it is needed to reform the DofSet at each
             * solution step or if it is possible to form it just once
             * - true  => reforme at each time step
             * - false => form just one (more efficient)
             * 
             * Default = false
             */
            bool mReformDofSetAtEachStep;
        
            /** 
             * Flag telling if it is needed or not to compute the reactions
             * 
             * default = true
             */
            bool mCalculateReactionsFlag;
        
            bool mSolutionStepIsInitialized;
        
            /**
             * maximum number of newton-raphson iterations
             * default = 30
             */
            unsigned int mMaxIterationNumber;
        
            bool mInitializeWasPerformed;
            
            //**********************************************************************
            //**********************************************************************
            /**
             * Setting up the contact conditions.
             * This function performs the contact search and creates
             * temporary, virtual linking conditions for each integration point
             * on the slave surfaces.
             */
            int SetUpContactConditions()
            {
                BaseType::GetModelPart().GetProperties(1)[K_CONTACT]= INT_MAX;
                BaseType::GetModelPart().GetProperties(1)[K_CONTACT_T]= INT_MAX;
                //getting the array of the conditions
                ConditionsArrayType& ConditionsArray = BaseType::GetModelPart().Conditions();
                //setting up candidate master conditions
                ConditionsArrayType MasterConditionsArray;
                //setting up an array of linking conditions
                ConditionsArrayType LinkingConditions;
                
                int lastRealCondition = BaseType::GetModelPart().Conditions().size();
//                 KRATOS_WATCH( lastRealCondition );
                
                
                if(lastRealCondition != 0 )
                {
                
                    for( typename ConditionsArrayType::ptr_iterator it = 
                         ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end();
                         ++it )
                    {
						if( (*it)->GetValue( ACTIVATION_LEVEL ) == 0 )
						{
							if( (*it)->GetValue( IS_CONTACT_MASTER ) )
							{
// 								std::cout << "contact master detected..." << std::endl;
								MasterConditionsArray.push_back( *it );
							}
						}
                    }
                    
                    GeometryType::Pointer tempGeometry =  GeometryType::Pointer( new Geometry<Node<3> >() );
                
                    PropertiesType::Pointer tempProperties = ConditionsArray.begin()->pGetProperties();
                
                    typename ConditionsArrayType::ptr_iterator it_begin=ConditionsArray.ptr_begin();
                    typename ConditionsArrayType::ptr_iterator it_end=ConditionsArray.ptr_end();
                    for (typename ConditionsArrayType::ptr_iterator it = it_begin; it!=it_end; ++it)
                    {
						if( (*it)->GetValue( ACTIVATION_LEVEL ) == 0 )
						{
							if( (*it)->GetValue( IS_CONTACT_SLAVE ) )
							{
								for( IndexType i = 0; i < (*it)->GetGeometry().IntegrationPoints().size();
															 i++ )
								{
									Point<3> MasterContactLocalPoint;
									Point<3> SlaveContactLocalPoint;
									(*it)->GetValue(PENALTY)[i]=BaseType::GetModelPart().GetProperties(1)[INITIAL_PENALTY];
									(*it)->GetValue(PENALTY_T)[i]=BaseType::GetModelPart().GetProperties(1)[INITIAL_PENALTY_T];
									Condition::Pointer CurrentMaster ;
									if( SearchPartner( 
																   (**it),
									MasterConditionsArray, i,
									MasterContactLocalPoint,
									SlaveContactLocalPoint,
									CurrentMaster
													 )
									  )
									{
										IndexType newId = (BaseType::GetModelPart().Conditions().end()-1)->Id()+LinkingConditions.size()+1;
										//creating contact link element
										Condition::Pointer newLink = Condition::Pointer( new 
												ContactLink3D(
												newId,
										tempGeometry,
										tempProperties,
										CurrentMaster, 
										*it,
										MasterContactLocalPoint,
										SlaveContactLocalPoint, 
										i
															 ) );
										LinkingConditions.push_back( newLink );
									}
								}
							}
						}
                    }
                    
                    if( BaseType::GetModelPart().GetProperties(1)[CONTACT_DOUBLE_CHECK] )
                    {
                        if( this->mEchoLevel > 1 )
                        {
                            std::cout << "double-check of contact problem activated: additional contact links are set up..." << std::endl;
                        }
                        MasterConditionsArray.clear();
                        for( typename ConditionsArrayType::ptr_iterator it = 
                             ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end();
                             ++it )
                        {
							if( (*it)->GetValue( ACTIVATION_LEVEL ) == 0 )
							{
								if( (*it)->GetValue( IS_CONTACT_SLAVE ) )
								{
									MasterConditionsArray.push_back( *it );
								}
							}
                        }
                    
                        GeometryType::Pointer tempGeometry =  GeometryType::Pointer( new Geometry<Node<3> >() );
                
                        PropertiesType::Pointer tempProperties = ConditionsArray.begin()->pGetProperties();
                
                        typename ConditionsArrayType::ptr_iterator it_begin=ConditionsArray.ptr_begin();
                        typename ConditionsArrayType::ptr_iterator it_end=ConditionsArray.ptr_end();
                        for (typename ConditionsArrayType::ptr_iterator it = it_begin; it!=it_end; ++it)
                        {
							if( (*it)->GetValue( ACTIVATION_LEVEL ) == 0 )
							{
								if( (*it)->GetValue( IS_CONTACT_MASTER ) )
								{
									for( IndexType i = 0; i < (*it)->GetGeometry().IntegrationPoints().size();
																	i++ )
									{
										Point<3> MasterContactLocalPoint;
										Point<3> SlaveContactLocalPoint;
										Condition::Pointer CurrentMaster ;
										if( SearchPartner( 
																		  (**it),
										MasterConditionsArray, i,
										MasterContactLocalPoint,
										SlaveContactLocalPoint,
										CurrentMaster
														 )
										  )
										{
											IndexType newId = (BaseType::GetModelPart().Conditions().end()-1)->Id()+LinkingConditions.size()+1;
											Condition::Pointer newLink = Condition::Pointer( new 
													ContactLink3D(
													newId,
											tempGeometry,
											tempProperties,
											CurrentMaster, 
											*it,
											MasterContactLocalPoint,
											SlaveContactLocalPoint, 
											i
																 ) );
											LinkingConditions.push_back( newLink );
										}
									}
								}
							}
                        }
                    }
                    
                    //adding linking to model_part
                    for( typename ConditionsArrayType::ptr_iterator it=LinkingConditions.ptr_begin();
                         it != LinkingConditions.ptr_end(); ++it )
                    {
                        BaseType::GetModelPart().Conditions().push_back( *it );
                    }
                    LinkingConditions.clear();
                }
                return lastRealCondition;
            }//SetUpContactConditions
            
            /**
             * searches a contact partner for a given slave condition
             */
            bool SearchPartner(
                    Condition& Slave,
                    ConditionsArrayType& AllMasterElements, 
                    const IndexType& IntegrationPointIndex,
                    Point<3>& MasterContactLocalPoint,
                    Point<3>& SlaveContactLocalPoint,
                    Condition::Pointer& CurrentMaster
                              )
            {
               
                KRATOS_TRY
                bool PartnerExists = false;
                
                //checking for existent master surfaces
                if( AllMasterElements.size() > 0 )
                {
                    PartnerExists = true;
                    
                    SlaveContactLocalPoint 
                            = Slave.GetGeometry().IntegrationPoints()[IntegrationPointIndex];
//                    KRATOS_WATCH(SlaveContactLocalPoint);
                    //calculating global coordinates of current integration point 
                    Point<3> SlaveContactGlobalPoint;
                    SlaveContactGlobalPoint = Slave.GetGeometry().GlobalCoordinates(
                            SlaveContactGlobalPoint, SlaveContactLocalPoint );
//                    KRATOS_WATCH(SlaveContactGlobalPoint);
				   
                    Point<3> GlobalCandidate;
//                     KRATOS_WATCH( GlobalCandidate );
                    
                    //defining set of possible master surface elements
                    ConditionsArrayType::Pointer MasterSet( 
                            new ConditionsArrayType() );
                    double minDist = INT_MAX;
                    //loop over all master surfaces (global search)
                    for( ConditionsArrayType::ptr_iterator it =
                         AllMasterElements.ptr_begin(); 
                         it != AllMasterElements.ptr_end();
                         ++it )
                    {
                        //loop over all nodes in current master surface
                        for( unsigned int n=0; n<(*it)->GetGeometry().PointsNumber(); n++ )
                        {
                            double dist = (((*it)->GetGeometry().GetPoint(n).X()-SlaveContactGlobalPoint[0])
                                        * ((*it)->GetGeometry().GetPoint(n).X()-SlaveContactGlobalPoint[0])
                                        + ((*it)->GetGeometry().GetPoint(n).Y()-SlaveContactGlobalPoint[1])
                                        * ((*it)->GetGeometry().GetPoint(n).Y()-SlaveContactGlobalPoint[1])
                                        + ((*it)->GetGeometry().GetPoint(n).Z()-SlaveContactGlobalPoint[2])
                                        * ((*it)->GetGeometry().GetPoint(n).Z()-SlaveContactGlobalPoint[2])
                                          );
                            if( fabs(dist-minDist) < 1e-4 )
                            {
                                MasterSet->push_back(*it);
//                                KRATOS_WATCH(minDist);
//                                 std::cout << "candidate master element added..." << std::endl;
                            } 
                            else if( dist < minDist )
                            {
                                MasterSet->clear();
                                
//                                 std::cout << "MasterSet is cleared: lower distance found" << std::endl;
                                
                                GlobalCandidate = (*it)->GetGeometry().GetPoint(n);
//                                 KRATOS_WATCH( GlobalCandidate );
                                minDist = dist;
                                MasterSet->push_back(*it);
                            }
                        }
//                         std::cout << "minimum distance: " << minDist << std::endl;
                    }
//                     std::cout << "number of checked master elements: " << MasterSet->size() << std::endl;
                    //searching contact partner (local search)
                    Point<3> MasterContactGlobalPoint;
                    bool LocalPartnerExists = false;
//                    KRATOS_WATCH( GlobalCandidate );
                    for( ConditionsArrayType::ptr_iterator it = MasterSet->ptr_begin(); it != MasterSet->ptr_end(); ++it )
                    {
                        if( ClosestPoint( *it, MasterContactGlobalPoint, 
                            MasterContactLocalPoint, SlaveContactGlobalPoint, GlobalCandidate ) )
                        {
// 							std::cout << "Found contact partner: " << MasterContactGlobalPoint << std::endl;
// 							std::cout << "for current Point: " << SlaveContactGlobalPoint << std::endl;
                            CurrentMaster = *it;
                            LocalPartnerExists = true;
                            break;
                        }
                
                    }
                    PartnerExists = ( PartnerExists && LocalPartnerExists ); 
                }
        
                return PartnerExists;
                KRATOS_CATCH("")
            }//SearchPartner
            
            /**
             * This function calculates the local coordinates of an orthogonal 
             * projection on an arbitrary surface.
             */
            bool ClosestPoint( Condition::Pointer& Surface,
                               GeometryType::CoordinatesArrayType& rResultGlobal, 
                               GeometryType::CoordinatesArrayType& rResultLocal,
                               const GeometryType::CoordinatesArrayType& rSlaveContactGlobalPoint,
                               const GeometryType::CoordinatesArrayType& rCandidateGlobal
                             )
            {
                double Xi1 = 0.0;
                double Xi2 = 0.0;
                double deltaXi1 = 0.0;
                double deltaXi2 = 0.0;
                Matrix localCoords;
//                 KRATOS_WATCH( Surface->GetValue( IS_CONTACT_MASTER ) );
//                 KRATOS_WATCH( Surface->GetGeometry() );
                localCoords = (Surface->GetGeometry()).PointsLocalCoordinates( localCoords );
                //determining local coordinates for rResult
                for( unsigned int n=0; n<Surface->GetGeometry().PointsNumber(); n++ )
                {
                    if(    fabs(rCandidateGlobal[0]-Surface->GetGeometry().GetPoint(n).X()) < 1e-7
                           && fabs(rCandidateGlobal[1]-Surface->GetGeometry().GetPoint(n).Y()) < 1e-7
                           && fabs(rCandidateGlobal[2]-Surface->GetGeometry().GetPoint(n).Z()) < 1e-7
                      )
                    {
                        Xi1 = localCoords(n,0);
                        Xi2 = localCoords(n,1);
                        break;
                    }
                }
                //setting up LocalPoint
                rResultLocal[0] = Xi1;
                rResultLocal[1] = Xi2;
                rResultLocal[2] = 0.0;
                //setting up rResult
                rResultGlobal = rCandidateGlobal;
//                 KRATOS_WATCH( rCandidateGlobal );
                //searching for orthogonal projection
                for( int k=0; k<1000; k++ )
				{
                    //setting up tangential vectors
                    Vector t1 = ZeroVector(3);//first tangential vector
                    Vector t2 = ZeroVector(3);//second tangential vector
                    //derivatives of tangential vectors
                    Vector dt11 = ZeroVector(3);
                    Vector dt12 = ZeroVector(3);
                    Vector dt21 = ZeroVector(3);
                    Vector dt22 = ZeroVector(3);
            
                    //retrieving first order derivatives in current solution point 
                    Matrix DN = ZeroMatrix(Surface->GetGeometry().PointsNumber(),2);
                    Surface->GetGeometry().ShapeFunctionsLocalGradients( DN, rResultLocal );
                    //retrieving second order derivatives in current solution point 
                    GeometryType::ShapeFunctionsSecondDerivativesType D2N;
                    Surface->GetGeometry().ShapeFunctionsSecondDerivatives( D2N, rResultLocal );
                    for( unsigned int n=0; n<Surface->GetGeometry().PointsNumber(); n++ )
                    {
                        //contribution to tangential vectors
                        t1[0] += Surface->GetGeometry().GetPoint(n).X()*DN(n,0);
                        t1[1] += Surface->GetGeometry().GetPoint(n).Y()*DN(n,0);
                        t1[2] += Surface->GetGeometry().GetPoint(n).Z()*DN(n,0);
                        t2[0] += Surface->GetGeometry().GetPoint(n).X()*DN(n,1);
                        t2[1] += Surface->GetGeometry().GetPoint(n).Y()*DN(n,1);
                        t2[2] += Surface->GetGeometry().GetPoint(n).Z()*DN(n,1);
                        //contribution to derivatives of tangential vectors
                        dt11[0] += Surface->GetGeometry().GetPoint(n).X()*D2N[n](0,0);
                        dt11[1] += Surface->GetGeometry().GetPoint(n).Y()*D2N[n](0,0);
                        dt11[2] += Surface->GetGeometry().GetPoint(n).Z()*D2N[n](0,0);
                        dt12[0] += Surface->GetGeometry().GetPoint(n).X()*D2N[n](0,1);
                        dt12[1] += Surface->GetGeometry().GetPoint(n).Y()*D2N[n](0,1);
                        dt12[2] += Surface->GetGeometry().GetPoint(n).Z()*D2N[n](0,1);
                        dt21[0] += Surface->GetGeometry().GetPoint(n).X()*D2N[n](1,0);
                        dt21[1] += Surface->GetGeometry().GetPoint(n).Y()*D2N[n](1,0);
                        dt21[2] += Surface->GetGeometry().GetPoint(n).Z()*D2N[n](1,0);
                        dt22[0] += Surface->GetGeometry().GetPoint(n).X()*D2N[n](1,1);
                        dt22[1] += Surface->GetGeometry().GetPoint(n).Y()*D2N[n](1,1);
                        dt22[2] += Surface->GetGeometry().GetPoint(n).Z()*D2N[n](1,1);
                    }
                    //defining auxiliary terms
                    double A1 = ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*t1[0])
                                +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*t1[1])
                                +((rSlaveContactGlobalPoint[2]-rResultGlobal[2])*t1[2]);
                    double A2 = ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*t2[0])
                                +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*t2[1])
                                +((rSlaveContactGlobalPoint[2]-rResultGlobal[2])*t2[2]);
                    double B11 = (-t1[0]*t1[0]-t1[1]*t1[1]-t1[2]*t1[2])
                                + ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*dt11[0])
                                +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt11[1])
                                +((rSlaveContactGlobalPoint[2]-rResultGlobal[2])*dt11[2]);
                    double B12 = (-t2[0]*t1[0]-t2[1]*t1[1]-t2[2]*t1[2])
                                + ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*dt12[0])
                                +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt12[1])
                                +((rSlaveContactGlobalPoint[2]-rResultGlobal[2])*dt12[2]);
                    double B21 = (-t1[0]*t2[0]-t1[1]*t2[1]-t1[2]*t2[2])
                                + ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*dt21[0])
                                +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt21[1])
                                +((rSlaveContactGlobalPoint[2]-rResultGlobal[2])*dt21[2]);
                    double B22 = (-t2[0]*t2[0]-t2[1]*t2[1]-t2[2]*t2[2])
                                + ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*dt22[0])
                                +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt22[1])
                                +((rSlaveContactGlobalPoint[2]-rResultGlobal[2])*dt22[2]);
                    //calculating update for Xi
                    deltaXi1 = -A1*B22/(B11*B22-B12*B21)+A2*B12/(B11*B22-B12*B21);
                    deltaXi2 =  A2*B21/(B11*B22-B12*B21)-A2*B11/(B11*B22-B12*B21);
                    //updating Xi
                    Xi1 += deltaXi1;
                    Xi2 += deltaXi2;
                    //updating LocalPoint
                    rResultLocal[0] = Xi1;
                    rResultLocal[1] = Xi2;
                    //updating rResult
                    rResultGlobal = ZeroVector( 3 );
                    rResultGlobal = Surface->GetGeometry().GlobalCoordinates( rResultGlobal, rResultLocal );
					if( fabs(deltaXi1) < 1e-7 && fabs(deltaXi2) < 1e-7 )
                    {
                        //check whether contact point lies within elementary boundaries
                        if( fabs(Xi1) <= 1.0 && fabs(Xi2) <= 1.0 )
                        {
//                     std::cout << "found matching point after " << k << " iteratons" << std::endl;
                            return true;
                        }
                        else
                        {
//                     std::cout << "no matching point inside this element" << std::endl;
                            return false;
                        }
                    }
                }
                return false;
            }
            
            
            /**
             * This function cleans up the conditions list after uzawa iteration
             */
            void Clean( int lastRealCondition )
            {
                /**
                 * cleaning up model part 
                 */
                BaseType::GetModelPart().Conditions().erase(
                        BaseType::GetModelPart().Conditions().begin()+lastRealCondition,
                BaseType::GetModelPart().Conditions().end() );
				/**
				 * clear, if necessary, the strategies
				 */
				if(mReformDofSetAtEachStep == true )
				{
// 					std::cout << "Clearing System" << std::endl;
					TSparseSpace::Clear(mA);
					TSparseSpace::Clear(mDx);
					TSparseSpace::Clear(mb);
				}
                
            }
            
            /**
             * This function updates the lagrangian multipliers.
             */
            void Update( int lastRealCondition )
            {
                
                ConditionsArrayType& ConditionsArray = BaseType::GetModelPart().Conditions();
                
                /**
                 * Update of normal stress lagrange multipliers
                 */
                for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
                {
                    if( (*it)->GetValue( IS_CONTACT_SLAVE ) )
                    {
                        for( IndexType IntPoint = 0; 
                             IntPoint < (*it)->GetGeometry().IntegrationPoints().size();
                             ++IntPoint )
                        {

                            double newLambda = (*it)->GetValue( LAMBDAS )[IntPoint]
                                        +(*it)->GetValue( GAPS )[IntPoint]
                                        *(*it)->GetValue( PENALTY )[IntPoint];
                            if(newLambda < 0.0)
                            {
                                newLambda=0.0;
                            }
//                             std::cout << "before update: ";
//                             KRATOS_WATCH( (*it)->GetValue( LAMBDAS )[IntPoint] );
  
                            (*it)->GetValue( DELTA_LAMBDAS )[IntPoint] = newLambda
                                    - (*it)->GetValue( LAMBDAS )[IntPoint];
                            (*it)->GetValue( LAMBDAS )[IntPoint] = newLambda;
//                             std::cout << "after update: ";
//                             KRATOS_WATCH( (*it)->GetValue( LAMBDAS )[IntPoint] );
                        }
                    }
                }
                
                 /**
                 * Update of frictional stress lagrange multipliers
                 * 
                 */


//                 KRATOS_WATCH(BaseType::GetModelPart().GetMesh().GetProperties(2)[FRICTION_COEFFICIENT] );
                if( BaseType::GetModelPart().GetProperties(1)[FRICTION_COEFFICIENT] > 0.0 )
                {
                    for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin()+lastRealCondition; 
                         it!=ConditionsArray.ptr_end(); ++it)
                    {
                       Vector newLambda_T_trial=ZeroVector(2);

                        Matrix m = (*it)->GetValue(CONTACT_LINK_M);
                        
                        newLambda_T_trial[0] =
                                (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(LAMBDAS_T)( (*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX),0) 
                                + (*it)->GetValue( CONTACT_LINK_SLAVE)->GetValue( PENALTY_T )[(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX)] 
                                * (m(0,0) 
                                * ((*it)->GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT)[0] 
                                - (*it)->GetValue( MASTER_CONTACT_LOCAL_POINT)[0])
                                + m(0,1)
                                * ((*it)->GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT)[1]
                                - (*it)->GetValue( MASTER_CONTACT_LOCAL_POINT)[1]));
                             
                        newLambda_T_trial[1] =
                                (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(LAMBDAS_T)( (*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX),1)
                                + (*it)->GetValue( CONTACT_LINK_SLAVE)->GetValue( PENALTY_T )[(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX)]
                                * (m(1,0) 
                                * ((*it)->GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT)[0] 
                                - (*it)->GetValue( MASTER_CONTACT_LOCAL_POINT)[0])
                                + m(1,1)
                                * ((*it)->GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT)[1]
                                - (*it)->GetValue( MASTER_CONTACT_LOCAL_POINT)[1]));
                        
                        Matrix mInv = ZeroMatrix(2);
                        double det=0.0;
                        MathUtils<double>::InvertMatrix2( m, mInv, det );

                        double NormLambda_T = sqrt( newLambda_T_trial[0]*mInv(0,0)*newLambda_T_trial[0]
                                    + newLambda_T_trial[0]*mInv(0,1)*newLambda_T_trial[1]
                                    + newLambda_T_trial[1]*mInv(1,0)*newLambda_T_trial[0]
                                    + newLambda_T_trial[1]*mInv(1,1)*newLambda_T_trial[1]);//;

                        
                        if( ( NormLambda_T <= (((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS )[(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX)]) )*(*it)->GetProperties()[FRICTION_COEFFICIENT]))
                        {
                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( DELTA_LAMBDAS_T )( (*it)->GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX),0) = newLambda_T_trial[0]-(*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX),0);
                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( DELTA_LAMBDAS_T )( (*it)->GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX),1) = newLambda_T_trial[1]-(*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX),1);
                            
                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX),0) = newLambda_T_trial[0];
                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX),1) = newLambda_T_trial[1];
                        }
                        else
                        {
                            newLambda_T_trial[0] = (*it)->GetProperties()[FRICTION_COEFFICIENT]
                                    * (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS )[(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX)]
                                    *newLambda_T_trial[0]/(NormLambda_T);

                            
                            newLambda_T_trial[1] = (*it)->GetProperties()[FRICTION_COEFFICIENT]
                                    * (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS )[(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX)]
                                    *newLambda_T_trial[1]/(NormLambda_T);
                            
                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( DELTA_LAMBDAS_T )( (*it)->GetValue(         
                                    CONTACT_SLAVE_INTEGRATION_POINT_INDEX),0) = newLambda_T_trial[0]-(*it)->GetValue( 
                                    CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue( 
                                    CONTACT_SLAVE_INTEGRATION_POINT_INDEX),0);
                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( DELTA_LAMBDAS_T )( (*it)->GetValue( 
                                    CONTACT_SLAVE_INTEGRATION_POINT_INDEX),1) = newLambda_T_trial[1]-(*it)->GetValue( 
                                    CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue( 
                                    CONTACT_SLAVE_INTEGRATION_POINT_INDEX),1);

                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue(
                                    CONTACT_SLAVE_INTEGRATION_POINT_INDEX),0) = newLambda_T_trial[0];
                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue(
                                    CONTACT_SLAVE_INTEGRATION_POINT_INDEX),1) = newLambda_T_trial[1];
                        }
                    }
                }
                
                if(BaseType::GetModelPart().GetProperties(1)[CONTACT_RAMP])
                {
					std::cout << "ramping penalties" << std::endl;
                   
                double alpha=BaseType::GetModelPart().GetProperties(1)[RAMP_CRITERION];
                double alpha_T=BaseType::GetModelPart().GetProperties(1)[RAMP_CRITERION_T];
                double beta=BaseType::GetModelPart().GetProperties(1)[RAMP_FACTOR];
                double beta_T=BaseType::GetModelPart().GetProperties(1)[RAMP_FACTOR_T];
                double Kmax=0.0;
                double Kmax_T= 0.0;
                double max_penalty=BaseType::GetModelPart().GetProperties(1)[MAXIMUM_PENALTY];
                double max_penalty_T=BaseType::GetModelPart().GetProperties(1)[MAXIMUM_PENALTY];
                bool check_penalty = false;
                bool check_penalty_T = false;
                    for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
                    {
                        if( (*it)->GetValue( IS_CONTACT_SLAVE ) )
                        {
                            for( IndexType IntPoint = 0; 
                                IntPoint < (*it)->GetGeometry().IntegrationPoints().size();
                                ++IntPoint )
                            {
                                if(((*it)->GetValue( DELTA_LAMBDAS )[IntPoint]) > Kmax)
                                {
                                  Kmax = (*it)->GetValue( DELTA_LAMBDAS )[IntPoint];
                                }
                            }
                        }
                    }
                    if( BaseType::GetModelPart().GetProperties(1)[FRICTION_COEFFICIENT] > 0.0 )
                    {
                        for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin()+lastRealCondition; 
                         it!=ConditionsArray.ptr_end(); ++it)
                        {
                            int i=(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX);
                            Matrix m = (*it)->GetValue(CONTACT_LINK_M);
                            Matrix mInv = ZeroMatrix(2);
                            double det=0.0;
                            MathUtils<double>::InvertMatrix2( m, mInv, det );
                            double K_trial = sqrt((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                        *mInv(0,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                        +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                        *mInv(0,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                        +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                        *mInv(1,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                        +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                        *mInv(1,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(DELTA_LAMBDAS_T )(i,1));
                            if(K_trial > Kmax_T)
                            {
                                Kmax_T=K_trial;
                            }
                        }
                    }

                    if( Kmax >= (this)-> GetModelPart().GetProperties(1)[K_CONTACT])
                        {
                            for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
                            {
                                if( (*it)->GetValue( IS_CONTACT_SLAVE ) )
                                {
                                    for( IndexType IntPoint = 0; 
                                     IntPoint < (*it)->GetGeometry().IntegrationPoints().size();
                                     ++IntPoint )
                                    { 
                                        if( ((*it)->GetValue( DELTA_LAMBDAS )[IntPoint]) > 
                                               ((this)->GetModelPart().GetProperties(1)[K_CONTACT]/alpha))
                                        {
                                            (*it)->GetValue( PENALTY)[IntPoint] = beta*(*it)->GetValue( PENALTY )[IntPoint] ;
											std::cout << "updated penalty to " << (*it)->GetValue( PENALTY)[IntPoint] << std::endl;
                                            
                                            if((*it)->GetValue( PENALTY)[IntPoint] > max_penalty)
                                            {
                                                (*it)->GetValue( LAMBDAS )[IntPoint]=((*it)->GetValue( LAMBDAS )[IntPoint]
                                                        -(*it)->GetValue( DELTA_LAMBDAS )[IntPoint])
                                                         +(*it)->GetValue( DELTA_LAMBDAS )[IntPoint]
                                                        /(max_penalty/
                                                        ((*it)->GetValue( PENALTY)[IntPoint])*beta);
                                                (*it)->GetValue( PENALTY)[IntPoint]=max_penalty;
                                            }
                                            else
                                            {
                                            (*it)->GetValue( LAMBDAS )[IntPoint]=((*it)->GetValue( LAMBDAS )[IntPoint]
                                                    -(*it)->GetValue( DELTA_LAMBDAS )[IntPoint])
                                                     +(*it)->GetValue( DELTA_LAMBDAS )[IntPoint]
                                            /beta;
                                            }
                                        }
                                    }
                                } 
                         }
                         check_penalty =true;
                    }
                    
                    if( BaseType::GetModelPart().GetProperties(1)[FRICTION_COEFFICIENT] > 0.0 )  
                    {
                        if(Kmax_T >= (this)-> GetModelPart().GetProperties(1)[K_CONTACT_T])
                        {
                            for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin()+lastRealCondition; 
                                 it!=ConditionsArray.ptr_end(); ++it)
                            {
                                int i=(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX);
                                Matrix m = (*it)->GetValue(CONTACT_LINK_M);
                                Matrix mInv = ZeroMatrix(2);
                                double det=0.0;
                                MathUtils<double>::InvertMatrix2( m, mInv, det );
                                double K = sqrt((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            *mInv(0,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            *mInv(0,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                            +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                            *mInv(1,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                            *mInv(1,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(DELTA_LAMBDAS_T )(i,1));
                                
                                    
                                 if( K > ((this)-> GetModelPart().GetProperties(1)[K_CONTACT_T]/alpha_T))
                                        {
                                            (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i]=beta_T* (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T )[i] ;

                                            if((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i]>max_penalty_T)
                                            {
                                                (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)
                                                        -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0))
                                                        +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                                        /(max_penalty_T/((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i])*beta_T);
                                                (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)
                                                        -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1))
                                                        +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                                        /(max_penalty_T/((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i])*beta_T);
                                                (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i]=max_penalty_T;
                                            }
                                            else
                                            {
                                                (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)
                                                        -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0))
                                                        +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                                /beta_T;
                                                (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)
                                                        -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1))
                                                        +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                                /beta_T;
                                            }
                                        }
                            }
                            check_penalty_T =true;
                        }
                    }

                    
                    if(Kmax <= (this)-> GetModelPart().GetProperties(1)[K_CONTACT]/alpha || !check_penalty)
                    {
                        (this)->GetModelPart().GetProperties(1)[K_CONTACT] = Kmax;
                        check_penalty =true;
                    }
                    
                    if( BaseType::GetModelPart().GetProperties(1)[FRICTION_COEFFICIENT] > 0.0 )  
                    {
                        if(Kmax_T <= (this)-> GetModelPart().GetProperties(1)[K_CONTACT_T]/alpha_T || !check_penalty_T)
                        {
                            (this)->GetModelPart().GetProperties(1)[K_CONTACT_T] = Kmax_T;
                            check_penalty_T =true;
                        }
                    }
                    
                    if(!check_penalty)
                    {
                        for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
                        {
                            if( (*it)->GetValue( IS_CONTACT_SLAVE ) )
                            {
                                for( IndexType IntPoint = 0; 
                                     IntPoint < (*it)->GetGeometry().IntegrationPoints().size();
                                     ++IntPoint )
                                { 
                                    if(((*it)->GetValue( DELTA_LAMBDAS )[IntPoint]) > ((this)->GetModelPart().GetProperties(1)[K_CONTACT]/alpha))
                                    {
                                        (*it)->GetValue( PENALTY )[IntPoint]=beta* ((*it)->GetValue( PENALTY )[IntPoint]);
                                       
                                        if((*it)->GetValue( PENALTY)[IntPoint]>max_penalty)
                                        {
                                            (*it)->GetValue( LAMBDAS )[IntPoint]=((*it)->GetValue( LAMBDAS )[IntPoint])
                                                     -(*it)->GetValue( DELTA_LAMBDAS )[IntPoint]
                                                     +(*it)->GetValue( DELTA_LAMBDAS )[IntPoint]
                                                    /(max_penalty/((*it)->GetValue( PENALTY)[IntPoint])*beta);
                                            (*it)->GetValue( PENALTY)[IntPoint]=max_penalty;
                                        }
                                        else
                                        {
                                            (*it)->GetValue( LAMBDAS )[IntPoint]=(*it)->GetValue( LAMBDAS )[IntPoint]
                                                     -(*it)->GetValue( DELTA_LAMBDAS )[IntPoint]
                                                     +(*it)->GetValue( DELTA_LAMBDAS )[IntPoint]
                                                    /beta;
                                        }
                                    } 
                                }
                            } 
                        }
                        (this)->GetModelPart().GetProperties(1)[K_CONTACT] = Kmax;
                        check_penalty =true;
                    }
                    if( BaseType::GetModelPart().GetProperties(1)[FRICTION_COEFFICIENT] > 0.0 )  
                    {
                        if(!check_penalty_T)
                        {
                            for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin()+lastRealCondition; 
                                 it!=ConditionsArray.ptr_end(); ++it)
                            {
                                int i=(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX);
                                Matrix m = (*it)->GetValue(CONTACT_LINK_M);
                                Matrix mInv = ZeroMatrix(2);
                                double det=0.0;
                                MathUtils<double>::InvertMatrix2( m, mInv, det );
                                double K = sqrt((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            *mInv(0,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            *mInv(0,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                            +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                            *mInv(1,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                            *mInv(1,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(DELTA_LAMBDAS_T )(i,1));
                                
                                        if(K > ((this)->GetModelPart().GetProperties(1)[K_CONTACT_T]/alpha_T) )
                                        {
                                            (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T )[i] = beta_T* ((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T )[i]);
                                            if((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i] > max_penalty_T)
                                            {
                                                (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0))
                                                        -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                                        +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                                        /(max_penalty_T/((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i])*beta_T);
                                                (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1))
                                                        -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                                        +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                                    /(max_penalty_T/((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i])*beta_T);
                                                (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i] = max_penalty_T;
                                            }
                                            else
                                            {
                                                (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0))
                                                        -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                                        +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                                   /beta_T;
                                                (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1))
                                                        -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                                        +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                                   /beta_T;
                                            }
                                        } 
							
                            }
                            (this)->GetModelPart().GetProperties(1)[K_CONTACT_T] = Kmax_T;
                            check_penalty =true;
                        }
                    }
                }//if(BaseType::GetModelPart().GetMesh().GetProperties(1)[CONTACT_RAMP])
                    return;        
            }//UpdateAndClean
            
            //**********************************************************************
            //**********************************************************************
            void Initialize()
            {
//                 KRATOS_TRY
                //pointers needed in the solution
//                 typename TSchemeType::Pointer pScheme = GetScheme();
//                 typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;
//                 ProcessInfo& pCurrentProcessInfo = GetModelPart().GetProcessInfo();
            
                //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
//                 if (pScheme->SchemeIsInitialized() == false) 
//                     pScheme->Initialize(GetModelPart());
            
                //Initialize The Elements - OPERATIONS TO BE DONE ONCE
//                 if (pScheme->ElementsAreInitialized() == false) 
//                     pScheme->InitializeElements(GetModelPart());
            
                //initialisation of the convergence criteria
//                 if (mpConvergenceCriteria->mConvergenceCriteriaIsInitialized == false)
//                     mpConvergenceCriteria->Initialize(GetModelPart());
            
//                 KRATOS_CATCH("")
            }
        
            //**********************************************************************
            //**********************************************************************
            void InitializeSolutionStep()
            {
//                 KRATOS_TRY
            
//                         typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
//                 typename TSchemeType::Pointer pScheme = GetScheme();
            
//                 ProcessInfo& pCurrentProcessInfo = GetModelPart().GetProcessInfo();
            
                //setting up the Vectors involved to the correct size 
//                 pBuilderAndSolver->ResizeAndInitializeVectors( mA, mDx, mb, 
//                         GetModelPart().Elements(), GetModelPart().Conditions(),
//                         GetModelPart().GetProcessInfo());
            
                //initial operations ... things that are constant over the Solution Step
//                 pBuilderAndSolver->InitializeSolutionStep(GetModelPart(),mA,mDx,mb);
//                 if(GetOldResidual().size() == 0)
//                 {
//                     GetOldResidual().resize(mb.size());
//                     TSparseSpace::SetToZero(GetOldResidual());
//                 }
            
                //initial operations ... things that are constant over the Solution Step
//                 pScheme->InitializeSolutionStep(GetModelPart(),mA,mDx,mb);
            
//                 KRATOS_CATCH("")
            }
        
            //**********************************************************************
            //**********************************************************************
            void MaxIterationsExceeded()
            {
                std::cout << "***************************************************" << std::endl;
                std::cout << "****** ATTENTION: max uzawa steps exceeded ********" << std::endl;
                std::cout << "***************************************************" << std::endl;
            }
        
            //**********************************************************************
            //**********************************************************************
            TSystemVectorType&  GetOldResidual()
            {
                return mOldResidualVector;
            }
        
            //**************************************************************************
            //**************************************************************************
            //to save the final RHS...
            void SetOldResidual(TSystemVectorType& rResidual)
            {
            //mResidual = b;
                mOldResidualVector = rResidual;
            //TSparseSpace::SetToZero(mOldResidualVector);
            //mOldResidualVector.resize(rResidual.size());
            //std::copy(rResidual.begin(), rResidual.end(), mOldResidualVector.begin());
            }
            
            /**
             * Copy constructor. 
             */
            ResidualBasedUzawaNewtonRaphsonStrategy( 
                    const ResidualBasedUzawaNewtonRaphsonStrategy& Other);
    };/* Class ResidualBasedUzawaNewtonRaphsonStrategy */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_UZAWA_NEWTON_RAPHSON_STRATEGY  defined */

