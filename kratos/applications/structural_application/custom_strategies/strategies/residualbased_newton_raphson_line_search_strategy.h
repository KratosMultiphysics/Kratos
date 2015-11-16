// Kratos Multi-Physics
// 
// Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement: 
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// 	
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY 
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/* *********************************************************
 *
 *   Last Modified by:    $Author: Nelson Lafontaine //inglafontaine@gmail.com    $
 *   Date:                $Date: 2009-28-04 $
 *   Revision:            $Revision: 1.20   $
 *
 * ***********************************************************/


#if !defined(KRATOS_RESIDUALBASED_NEWTON_RAPHSON_LINE_SEARCHES_STRATEGY)
#define  KRATOS_RESIDUALBASED_NEWTON_RAPHSON_LINE_SEARCHES_STRATEGY


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
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "custom_utilities/line_searches_utility.h"

#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>



namespace Kratos
{

template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver // = LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedNewtonRaphsonLineSearchesStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>,
  public LineSearchesUtility<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedNewtonRaphsonLineSearchesStrategy);

    /**@name Type Definitions */
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    /** Constructor.
     */
    ResidualBasedNewtonRaphsonLineSearchesStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        unsigned int MaxNewtonRapshonIterations,
        unsigned int MaxLineSearchIterations,
        double tolls, // energy tolerance factor on LineSearch (0.8 is ok)
        double amp, // maximum amplification factor
        double etmxa, // maximum allowed step length
        double etmna, // minimum allowed step length
        bool CalculateReactions,
        bool ReformDofSetAtEachStep,
        bool MoveMeshFlag,
        bool ApplyLineSearches
    )
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part,
                pScheme,
                pNewLinearSolver,
                pNewConvergenceCriteria,
                MaxNewtonRapshonIterations,
                CalculateReactions,
                ReformDofSetAtEachStep,
                MoveMeshFlag)
        ,
        LineSearchesUtility<TSparseSpace, TDenseSpace, TLinearSolver>(
            MaxLineSearchIterations,
            tolls,
            amp,
            etmxa,
            etmna,
            MoveMeshFlag,
            ApplyLineSearches)
    {
        KRATOS_TRY
        mKeepSystemConstantDuringIterations = false;

        //this->SetMaxIterationNumber(MaxNewtonRapshonIterations);
        //this->SetParametersLineSearches(MaxLineSearchIterations,tolls,amp,etmxa,etmna,MoveMeshFlag,ApplyLineSearches);
        KRATOS_CATCH("")
    }

    /** Destructor.
     */
    virtual ~ResidualBasedNewtonRaphsonLineSearchesStrategy()
    {
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
        typename TSchemeType::Pointer pScheme = ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::GetBuilderAndSolver();

        //int solstep = pCurrentProcessInfo.GetCurrentSolutionStep();
        DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();

        //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        //if the operations needed were already performed this does nothing

        if (this->mInitializeWasPerformed == false)
            this->Initialize();

        //set up the system, operation performed just once unless it is required
        //to reform the dof set at each iteration
        if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false ||
                this->mReformDofSetAtEachStep == true)
        {
            //setting up the list of the DOFs to be solved
            pBuilderAndSolver->SetUpDofSet(pScheme, BaseType::GetModelPart());

            //shaping correctly the system
            pBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());
        }

        //prints informations about the current time
        if (this->GetEchoLevel() != 0)
        {
            std::cout << " " << std::endl;
            std::cout << "CurrentTime = " << BaseType::GetModelPart().GetProcessInfo()[TIME] << std::endl;
        }

        //updates the database with a prediction of the solution
        ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Predict();

        //initialize solution step
        if (this->mSolutionStepIsInitialized == false)
        {
            this->InitializeSolutionStep();
        }

        TSystemMatrixType& mA = *(this->mpA);
        TSystemVectorType& mDx = *(this->mpDx);
        TSystemVectorType& mb = *(this->mpb);


        //initializing the parameters of the Newton-Raphson cicle
        unsigned int iteration_number = 0;
        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool is_converged = false;
        bool ResidualIsUpdated = false;
        bool Satisfactory_Line_Search = false;


        //so = TSparseSpace::Dot(mDx,mb);
        TSystemVectorPointerType pX_old; // iteracion anterior
        TSystemVectorPointerType pDelta_p;

        // asignando tamaño
        pBuilderAndSolver->ResizeAndInitializeVectors(this->mpA, pX_old, pDelta_p, BaseType::GetModelPart().Elements(), BaseType::GetModelPart().Conditions(), BaseType::GetModelPart().GetProcessInfo());

        TSystemVectorType& X_old = *pX_old;
        TSystemVectorType& Delta_p = *pDelta_p;

        // setting to zero elements
        TSparseSpace::SetToZero(X_old);
        TSparseSpace::SetToZero(Delta_p);

        this->BackupDatabase(rDofSet, X_old);
        //true  = 1
        //false = 0

        while ((is_converged == false) && (iteration_number++ < (this->mMaxIterationNumber)))
        {

            //setting the number of iteration
            std::cout << "Newton_Rapshon_Iteration_Number:" << iteration_number << std::endl;
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
            is_converged = this->mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

            //call the linear system solver to find the correction mDx for the
            //it is not called if there is no system to solve
            if (SparseSpaceType::Size(mDx) != 0)
            {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
                {
                    if (iteration_number == 1 || GetKeepSystemConstantDuringIterations() == false)
                    {
                        //mA = 0.00;
                        TSparseSpace::SetToZero(mA);
                        TSparseSpace::SetToZero(mDx);
                        TSparseSpace::SetToZero(mb);

                        pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
                        TSparseSpace::Copy(mDx, Delta_p);
                    }
                    else
                    {
                        TSparseSpace::SetToZero(mDx);
                        TSparseSpace::SetToZero(mb);

                        pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
                        TSparseSpace::Copy(mDx, Delta_p);
                    }
                }
                else
                {
                    TSparseSpace::SetToZero(mDx);
                    TSparseSpace::SetToZero(mb);

                    pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
                    TSparseSpace::Copy(mDx, Delta_p);
                }
            }
            else
            {
                std::cout << "ATTENTION: No Free DOFs!! " << std::endl;
            }

            // Using LineSearch
            // booleano
            // KRATOS_WATCH(mb)
            //KRATOS_WATCH(mDx)
            if (this->mApplyLineSearches == true)
            {

               Satisfactory_Line_Search = this->LineSearches(BaseType::GetModelPart(),
                                          pScheme,
                                          pBuilderAndSolver,
                                          rDofSet,
                                          X_old, Delta_p, mDx, mb, mA);
               
//                 if ( Satisfactory_Line_Search== true)
//                       {
//                             std::cout<<"***************************************************"<<std::endl;
//                             std::cout<<"******Line Searches Has Succesfully Finished********"<<std::endl;
//                             std::cout<<"***************************************************"<<std::endl;
//                       }

                 
                std::cout << "Line-Search Step Factor:" << this->meta << std::endl;
                rDofSet = pBuilderAndSolver->GetDofSet();
                this->SetDatabaseToValue(rDofSet, X_old);
                TSparseSpace::Assign(mDx, this->meta, mDx);
                pScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

            }
            else
            {
                //Updating the results stored in the database
                rDofSet = pBuilderAndSolver->GetDofSet();
                pScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
            }

            //move the mesh if needed

            if (BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
            TSparseSpace::SetToZero(X_old);
            this->BackupDatabase(rDofSet, X_old);
            //KRATOS_WATCH(X_old);

            pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

            ResidualIsUpdated = false;

            if (is_converged == true)
            {

                if (this->mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {
                    TSparseSpace::SetToZero(mb);

                    pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb);
                    ResidualIsUpdated = true;
                    //std::cout << "mb is calculated" << std::endl;
                }

                is_converged = this->mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
            }


        }


        //plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= (this->mMaxIterationNumber))
            ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::MaxIterationsExceeded();

        //recalculate residual if needed
        // (note that some convergence criteria need it to be recalculated)
        if (ResidualIsUpdated == false)
        {
            TSparseSpace::SetToZero(mb);

            pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb);

            //std::cout << "mb is calculated" << std::endl;
        }

        //calculate reactions if required
        if (this->mCalculateReactionsFlag == true)
        {
            pBuilderAndSolver->CalculateReactions(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }
        //Finalisation of the solution step,
        //operations to be done after achieving convergence, for example the
        //Final Residual Vector (mb) has to be saved in there
        //to avoid error accumulation
        pScheme->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);
        pBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

        //Cleaning memory after the solution
        pScheme->Clean();

        //reset flags for next step
        this->mSolutionStepIsInitialized = false;

        //deallocate the systemvectors
        if (this->mReformDofSetAtEachStep == true)
        {

            SparseSpaceType::Clear(this->mpA);
            SparseSpaceType::Clear(this->mpDx);
            SparseSpaceType::Clear(this->mpb);
            //TSparseSpace::ClearData(mA);
            //TSparseSpace::ClearData(mDx);
            //TSparseSpace::ClearData(mb);
        }


        return 0.00;

        KRATOS_CATCH("")
    }

    void SetKeepSystemConstantDuringIterations(bool value)
    {
        mKeepSystemConstantDuringIterations = value;
    }

    bool GetKeepSystemConstantDuringIterations()
    {
        return mKeepSystemConstantDuringIterations;
    }

private:
    //flag to allow keeping system matrix constant during iterations
    bool mKeepSystemConstantDuringIterations;
};
} /* namespace Kratos.*/

#endif


