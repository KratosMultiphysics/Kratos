// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Nelson Lafontaine
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_BASED_RESIDUAL_ARC_LENGHT_STRATEGY )
#define  KRATOS_BASED_RESIDUAL_ARC_LENGHT_STRATEGY

/* System includes */
#include <limits>
#include<iostream>
#include<iomanip>

/* External includes */

/* Project includes */
// #include "structural_mechanics_application.h"
#include "includes/define.h"
#include "includes/kernel.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

// Default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "custom_utilities/line_searches_utility.hpp"

#include <cmath>

namespace Kratos
{
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class ResidualBasedArcLengthStrategy
    : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef ConvergenceCriteria<TSparseSpace,TDenseSpace> TConvergenceCriteriaType;

    typedef LineSearchesUtility<TSparseSpace, TDenseSpace, TLinearSolver> TlineSearchesType;

    /** Counted pointer of ClassName */
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

    typedef long double RealType;

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /************************************* CONSTRUCTOR *********************************/
    /***********************************************************************************/

    ResidualBasedArcLengthStrategy(
            ModelPart& model_part,
            typename TSchemeType::Pointer pScheme,
            typename TLinearSolver::Pointer pNewLinearSolver,
            typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
            unsigned int Ide,
            unsigned int MaxIterations,
            unsigned int MaxRecursive,
            RealType factor_delta_lmax,
            bool CalculateReactions     = true,
            bool ReformDofSetAtEachStep = true,
            bool MoveMeshFlag           = true
            )
        : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag),
         mAuxElementModelPart(model_part.GetOwnerModel().CreateModelPart("ResidualBasedArcLengthStrategy_AuxElementModelPart")),
         mAuxConditionModelPart(model_part.GetOwnerModel().CreateModelPart("ResidualBasedArcLengthStrategy_AuxConditionModelPart"))
    {
        KRATOS_TRY;
        
        

        // Set flags to default values
        SetMaxIterationNumber(MaxIterations);
        mCalculateReactionsFlag = CalculateReactions;
        InitializeAuxiliaryModelParts(model_part);

        mfactor_delta_lmax      = factor_delta_lmax;
        mIde	                = Ide;
        mMaxRecursive           = MaxRecursive;

        mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        // Saving the convergence criteria to be used
        mpConvergenceCriteria = pNewConvergenceCriteria;

        // Saving the scheme
        mpScheme = pScheme;

        // Saving the linear solver
        mpLinearSolver = pNewLinearSolver;

        // Setting up the default builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer
                             (
                                 new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(mpLinearSolver)
                             );

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed    = false;
        mInit                      = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // Be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        KRATOS_CATCH("");
    }

    /************************************* DESTRUCTOR **********************************/
    /***********************************************************************************/
    
    ~ResidualBasedArcLengthStrategy() override 
    {
        Model& current_model = BaseType::GetModelPart().GetOwnerModel();
        current_model.DeleteModelPart("ResidualBasedArcLengthStrategy_AuxElementModelPart");
        current_model.DeleteModelPart("ResidualBasedArcLengthStrategy_AuxConditionModelPart");
    }

    /************************************* OPERATIONS **********************************/
    /***********************************************************************************/

    //Set and Get Scheme ... containing Builder, Update and other
    void SetScheme(typename TSchemeType::Pointer pScheme )
    {
        mpScheme = pScheme;
    };

    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    // Set and Get the BuilderAndSolver
    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver )
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
    {
        mCalculateReactionsFlag = CalculateReactionsFlag;
    }
    bool GetCalculateReactionsFlag()
    {
        return mCalculateReactionsFlag;
    }

    void SetReformDofSetAtEachStepFlag(bool flag)
    {
        mReformDofSetAtEachStep = flag;
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
    }
    bool GetReformDofSetAtEachStepFlag()
    {
        return mReformDofSetAtEachStep;
    }

    void SetMaxIterationNumber(unsigned int  MaxIterationNumber)
    {
        mMaxIterationNumber = MaxIterationNumber;
    }

    unsigned int GetMaxIterationNumber()
    {
        return mMaxIterationNumber;
    }

    // Level of echo for the solving strategy
    // 0 -> Mute... no echo at all
    // 1 -> Printing time and basic informations
    // 2 -> Printing linear solver data
    // 3 -> Print of debug informations:
    // Echo of stiffness matrix, Dx, b...
    void SetEchoLevel(int Level) override
    {
        BaseType::mEchoLevel = Level;
        GetBuilderAndSolver()->SetEchoLevel(Level);
    }

    void VariablesArcLength()
    {
        KRATOS_TRY;

        mlambda_old        = 0.00;
        mlambda            = 1.00;
        mdelta_lambda      = 1.00;
        meta               = 1.00;

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Operation to predict the solution ... if it is not called a trivial predictor is used in which the
    * values of the solution step of interest are assumed equal to the old values
    */

    void Predict() override
    {
        KRATOS_TRY;

        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();

        TSystemMatrixType& mA  = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb  = *mpb;

        GetScheme()->Predict(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

        // Move the mesh if needed
        if(this->MoveMeshFlag() == true)
        {
            BaseType::MoveMesh();
        }

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It initialises the vector of auxiliar values
    * @return pAux: Vector with auxiliar values
    */

    void InitializeAuxVectors(TSystemVectorPointerType& pAux)
    {
      if (pAux == NULL) // If the pointer is not initialized initialize it to an empty matrix
      {
          TSystemVectorPointerType pNewAux = TSystemVectorPointerType(new TSystemVectorType(0));
          pAux.swap(pNewAux);
      }

      TSystemVectorType& Aux = *pAux;
      if(Aux.size() !=  GetBuilderAndSolver()->GetEquationSystemSize())
      {
          Aux.resize(GetBuilderAndSolver()->GetEquationSystemSize(), false);
      }
   }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It solves the problem
    */

    double Solve() override
    {
        KRATOS_TRY;

        //std::cout<<std::fixed<<std::setw(15)<<std::scientific<<std::setprecision(9);
        if (this->GetEchoLevel() > 0)
        {
            std::cout<<"************************************************************************"<<std::endl;
            std::cout<<"Begininning Arc Lenght Method.A Pseudo-Line Searches Included. Please Wait...."<<std::endl;
            std::cout<<"************************************************************************"<<std::endl;
        }

        TSystemVectorPointerType  pq;               //  Fext
        TSystemVectorPointerType  pSigma_q;         //  Displacement conditions
        TSystemVectorPointerType  pSigma_h;         //  Displacemenet produced due to the imbalance
        TSystemVectorPointerType  ph;	            //  Ortogonal component of h
        TSystemVectorPointerType  pe;               //  Out of balance load  lambda*Fext - Fint
        TSystemVectorPointerType  pE;               //  Lamda_old + Delta_lambda) * Fext
        TSystemVectorPointerType  pAux_q;
        TSystemVectorPointerType  pAux_h;
        TSystemVectorPointerType  pq_Inc_Aux;

        // Initialize member variables
        mIterationNumber              = 0;
        mReduceArcLenght              = true;

        // Initialize other variables
        RealType Ao                   = 0.00;
        RealType A                    = 1.00;
        //RealType aux                  = 0.00;
        RealType miu                  = 0.00;
        RealType g                    = 0.00;
        RealType toler                = 0.1;   // Residual tolerance
        RealType toler_l              = 0.001; // Lambda tolerance
        RealType res                  = 1.00;
        RealType num                  = 1.00;
        RealType den                  = 1.00;
        RealType lambda_error         = 1.00;
        RealType lambda_old_iter      = 0.00;
        RealType fact                 = 0.00;

        bool local_converged          = false;
        bool local_converged_e        = false; // Residual
        bool local_converged_h        = false; // Orthogonal residual
        bool local_converged_l        = false; // Lambda
        bool is_converged             = false;
        unsigned int recursive        = 0;

        //vector<RealType> Parameters;
//         RealType old_residual = 0.00;
//         RealType new_residual = 0.00;

        //unsigned int MaxLineSearchIter = 50;
        //RealType tolls                 = 0.80;
        //RealType amp                   = 1.618;
        //RealType etmxa                 = 3.0;
        //RealType etmna                 = 0.10;
        //bool MoveMeshFlag              = BaseType::MoveMeshFlag();
        //bool ApplyLineSearches         = true;

        // Line searches utility
        //TLineSearchesType Searches;
        //TLineSearchesType LineSearches(MaxLineSearchIter, tolls, amp, etmxa, etmna, MoveMeshFlag, ApplyLineSearches);

        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        //ModelPart& r_model_part = BaseType::GetModelPart();

        DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();

        // Creating models part for analysis
        InitializeAuxiliaryModelParts(BaseType::GetModelPart());
        mstep = BaseType::GetModelPart().GetProcessInfo()[STEP];

        if (this->GetEchoLevel() > 0)
        {
            std::cout<<" STEP NUMBER                   = " << mstep << std::endl;
            std::cout<<" DESIRED ITERATIONS            = " << mIde  << std::endl;
            std::cout<<" MAX. ITERATIONS               = " << mMaxIterationNumber << std::endl;
            std::cout<<" MAX. RECURSIVE ITERATIONS     = " << mMaxRecursive << std::endl;
            std::cout<<" CURRENT TIME                  = " << BaseType::GetModelPart().GetProcessInfo()[TIME] << std::endl;
        }

        // Initialisation of the convergence criteria and variables of arc lenght
        if(mInitializeWasPerformed == false)
        {
            Initialize();
        }

        // Set up the system, operation performed just once unless it is required to reform the dof set at each iteration
        if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false || mReformDofSetAtEachStep == true )
        {
            // Setting up the list of the DOFs to be solved
            pBuilderAndSolver->SetUpDofSet(pScheme,BaseType::GetModelPart());

            // Shaping correctly the system
            pBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());
        }

        // Updates the database with a prediction of the solution
        Predict();

        // Initialize solution step
        if (mSolutionStepIsInitialized == false)
        {
            InitializeSolutionStep();
        }

        // Initializing the local variables
        InitializeAuxVectors(pSigma_q);
        InitializeAuxVectors(pSigma_h);
        InitializeAuxVectors(ph);
        InitializeAuxVectors(pe);
        InitializeAuxVectors(pE);
        InitializeAuxVectors(pq);

        InitializeAuxVectors(pAux_q);
        InitializeAuxVectors(pAux_h);
        InitializeAuxVectors(pq_Inc_Aux);

        // Main data
        TSystemVectorType& mDelta_p      = *mpDelta_p;    /// P  current change
        TSystemVectorType& mDelta_pold   = *mpDelta_pold; /// P  =  u_(step+1)-u_(step)
        TSystemVectorType& mX_old        = *mpX_old;      /// old = positions X+u
        TSystemMatrixType& mA            = *mpA;
        TSystemVectorType& mDx           = *mpDx;
        TSystemVectorType& mb            = *mpb;

        /// Local axiliareis cvector
        TSystemVectorType& Sigma_q       = *pSigma_q;
        TSystemVectorType& Sigma_h       = *pSigma_h;
        TSystemVectorType& h             = *ph;
        TSystemVectorType& e             = *pe;
        TSystemVectorType& E             = *pE;
        TSystemVectorType& q             = *pq;
        TSystemVectorType& Aux_q         = *pAux_q;
        TSystemVectorType& Aux_h         = *pAux_h;
        TSystemVectorType& q_Inc_Aux     = *pq_Inc_Aux;

        //// Do nothing. It is called in order to have an order sequence
        //pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
        //is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        // Function to perform the building and the solving phase.
        if(BaseType::mRebuildLevel >1 || BaseType::mStiffnessMatrixIsBuilt == false)
        {
            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mDx);
            TSparseSpace::SetToZero(mb);
            TSparseSpace::SetToZero(Sigma_q);
            TSparseSpace::SetToZero(Sigma_h);
            TSparseSpace::SetToZero(h);
            TSparseSpace::SetToZero(E);
            TSparseSpace::SetToZero(e);
            TSparseSpace::SetToZero(q);
            TSparseSpace::SetToZero(Aux_q);
            TSparseSpace::SetToZero(Aux_h);
            TSparseSpace::SetToZero(q_Inc_Aux);

            pBuilderAndSolver->Build(pScheme,mAuxElementModelPart,mA,mb);
            pBuilderAndSolver->BuildRHS(pScheme,mAuxConditionModelPart, q);
            TSparseSpace::Copy(q ,Aux_q); // Aux = q;
        }
        else
        {
            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mDx);
            TSparseSpace::SetToZero(mb);
            TSparseSpace::SetToZero(Sigma_q);
            TSparseSpace::SetToZero(Sigma_h);
            TSparseSpace::SetToZero(h);
            TSparseSpace::SetToZero(E);
            TSparseSpace::SetToZero(e);
            TSparseSpace::SetToZero(q);
            TSparseSpace::SetToZero(Aux_q);
            TSparseSpace::SetToZero(Aux_h);
            TSparseSpace::SetToZero(q_Inc_Aux);

            pBuilderAndSolver->Build(pScheme,mAuxElementModelPart,mA,mb);
            pBuilderAndSolver->BuildRHS(pScheme,mAuxConditionModelPart, q);
            TSparseSpace::Copy(q ,Aux_q); // Aux = q;
        }

        // Out of balance force for the first iteration in each step.
        // mb is not zero in the second step.
        noalias(q_Inc_Aux) = q; //+ mb;
        //KRATOS_WATCH(q);

        //WARNING: For some reason q is set
        //KRATOS_WATCH(mA);
        TSparseSpace::Copy(q_Inc_Aux, Aux_q);
        pBuilderAndSolver->SystemSolve(mA, Sigma_q, q_Inc_Aux);
        TSparseSpace::Copy(Aux_q, q_Inc_Aux);
        //noalias(Sigma_q) += mDelta_pold; /// should be the total acumulated

        //Iteration Cicle... performed only for NonLinearProblems
        do
        {
            mIterationNumber = 0;
            if(recursive++ >= mMaxRecursive)
            {
                break;
            }
            while(  is_converged == false && mIterationNumber++ < mMaxIterationNumber)
            {
                // Setting the number of iteration

                if (this->GetEchoLevel() > 0)
                {
                    std::cout<<"\n STEP NUMBER       = " << mstep <<"  ITERATIONS NUMBER = " << mIterationNumber << "  RECURSIVE NUMBER = " << recursive << std::endl;
                }
                BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = mIterationNumber;

                // Setting variables in the begining of the iteraction
                pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
                meta = 1.00;

                local_converged = false;
                if(mIterationNumber == 1 && mInit == false)// mstep == 1)
                {
                    mInit            = true;
                    mlambda_old      = 0.00;
                    Ao               = TSparseSpace::Dot(Sigma_q, Sigma_q);        // Ao = inner_prod (Sigma_q, Sigma_q);
                    mdelta_l         = sqrt(2.00 * Ao * mdelta_lambda * mdelta_lambda);
                    mdelta_lold      = mdelta_l;
                    mdelta_lmax      = mdelta_l*mfactor_delta_lmax;                // Maximum size of the arc-length
                    TSparseSpace::Assign(mDelta_p, mdelta_lambda,Sigma_q);         // mDelta_p = mdelta_lambda*Sigma_q;
                    mdelta_lambda_old = mdelta_lambda;
                    TSparseSpace::Copy(mDelta_p ,mDelta_pold);                     // WARNING = only for the fisrt step and the first iteraction mDelta_pold      = mDelta_p;
                    TSparseSpace::Copy(mDelta_p, mDx);

                    //TSparseSpace::Copy(mRHS_cond,h);

                    if (this->GetEchoLevel() > 1)
                    {
                        std::cout << " Solution Formulation at Origin " << std::endl;
                        std::cout << "   Arc length      = " << mdelta_l << std::endl;
                        std::cout << "   A_o             = " << Ao << std::endl;
                        std::cout << "   Delta Lamba Old = " << mdelta_lambda_old     << std::endl;
                        std::cout << "   Delta Lamba     = " << mdelta_lambda << std::endl;
                        std::cout << "   Eta Factor      = " << meta << std::endl;
                        std::cout << "   Delta p      = " << mDelta_p << std::endl;
                     }

                    is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(),rDofSet,mA, mDx, q_Inc_Aux);
                    //mpConvergenceCriteria->GetConvergenceData(Parameters);
                    //old_residual = Parameters[0];
                }
                else if(mIterationNumber == 1 && mInit == true)// mstep != 1)
                {
                    //RealType aux1  =  TSparseSpace::Dot(mX_old, mX_old);
                    RealType aux2  =  TSparseSpace::Dot(mDelta_pold, mDelta_pold);   //inner_prod(mDelta_pold,mDelta_pold);
                    Ao           =  aux2/(mlambda_old * mlambda_old);
                    miu          =  mdelta_l/std::sqrt(aux2 + Ao*mdelta_lambda_old*mdelta_lambda_old);

                    //KRATOS_WATCH(miu);
                    //KRATOS_WATCH(Ao);
                    miu = (miu >= 1.00) ? 1.00: miu;

                    TSparseSpace::Assign(mDelta_p, miu, mDelta_pold);            // mDelta_p     = miu*mDelta_pold;
                    mdelta_lambda = miu * mdelta_lambda_old;
                    // TSparseSpace::Assign(mDelta_p, 1.00, Sigma_q);

                    if (this->GetEchoLevel() > 1)
                    {
                        std::cout << " Solution Iteration from  Converged Point "    << std::endl;
                        std::cout << "   Arc length      = " << mdelta_l             << std::endl;
                        std::cout << "   Lamda Old       = " << mlambda_old          << std::endl;
                        std::cout << "   A_o             = " << Ao                   << std::endl;
                        std::cout << "   Delta Lamba Old = " << mdelta_lambda_old    << std::endl;
                        std::cout << "   Delta Lamba     = " << mdelta_lambda        << std::endl;
                        std::cout << "   Eta Factor      = " << meta                 << std::endl;
                        std::cout << "   Miu Factor      = " << miu                  << std::endl;
                    }

                    TSparseSpace::Copy(mDelta_p, mDx);

                    is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(),rDofSet,mA, mDx, q_Inc_Aux);
                    //mpConvergenceCriteria->GetConvergenceData(Parameters);
                    //old_residual = Parameters[0];
                }
                else
                {
                    TSparseSpace::SetToZero(mDx);

                    // Compute Sigma_h
                    TSparseSpace::Copy(h, Aux_h); /// Aux = h;
                    pBuilderAndSolver->SystemSolve(mA, Sigma_h, h);
                    TSparseSpace::Copy(Aux_h ,h); /// h = Aux

                    // Recursive function
                    //KRATOS_WATCH(g)
                    //KRATOS_WATCH(Ao)
                    //KRATOS_WATCH(Sigma_h)
                    //KRATOS_WATCH(mA)
                    Recursive_Function_Arc_Length(pq, pSigma_q, pSigma_h, mpDx, g, Ao);
                }

                if (mIterationNumber != 1)
                {
                    lambda_old_iter = mlambda;
                }
                mlambda =  mlambda_old + mdelta_lambda;

                //KRATOS_WATCH(mlambda);
                //KRATOS_WATCH(mlambda_old);
                //KRATOS_WATCH(mdelta_lambda);

                //TSparseSpace::Assign(E, mlambda, q);
                pScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

                if(BaseType::MoveMeshFlag() == true)
                {
                    BaseType::MoveMesh();
                }

                TSparseSpace::SetToZero(mb);
                pBuilderAndSolver->BuildRHS(pScheme,mAuxElementModelPart, mb);
                TSparseSpace::ScaleAndAdd(mlambda, q, A, mb, e); // Residual convergence Fint-lambda*Fext
                //noalias(e) = E + mb; /// WARNING = in Kratos Fint is compted like -mb

                // Finalize the iteration
                pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);

                // Convergence for e
                den        = TSparseSpace::Dot(q, q);
                num        = std::sqrt(TSparseSpace::Dot(e,e));
                //KRATOS_WATCH(num)
                //KRATOS_WATCH(den)
                //KRATOS_WATCH(mlambda)

                res        = 100.00 * std::abs(num/(mlambda* std::sqrt(den)));
                if (this->GetEchoLevel() > 1)
                {
                    std::cout << " Convergence reached for e = " << res << "  Required Convergence = " << toler << std::endl;
                }

                if(res < toler)
                {
                    local_converged_e = true;
                }

                // Convergence for h
                // Computing g and h
                num        = TSparseSpace::Dot(e, q);
                g          = num/den;
                //noalias(h) = e - g*q;
                TSparseSpace::ScaleAndAdd(A, e,-g, q, h);

                if(std::sqrt(TSparseSpace::Dot(h,h)) < 1.0E-10)
                {
                    TSparseSpace::Copy(e,h);
                    g = 0.0;
                }

                num        = std::sqrt(TSparseSpace::Dot(h,h));
                res        = 100.00 * std::abs(num/(mlambda * std::sqrt(den)));
                if (this->GetEchoLevel() > 1)
                {
                    std::cout << " Convergence reached for h = " << res << "  Required Convergence = " << toler << std::endl;
                }

                if(res < toler)
                {
                    local_converged_h = true;
                }

                if(std::abs(lambda_old_iter) < 1E-10)
                {
                    fact = 1.00;
                }
                else
                {
                    fact = lambda_old_iter;
                }
                lambda_error    = std::abs(100.00 * (mlambda - lambda_old_iter)/(fact + 1E-20));
                local_converged_l = (lambda_error < toler_l) ?  true : false;
                if (this->GetEchoLevel() > 1)
                {
                    std::cout << " Convergence reached for l = " << lambda_error << "  Required Convergence = " << toler_l << std::endl;
                }

                // Compute the truth criteria of convergence
                local_converged = bool(local_converged_e || local_converged_h || local_converged_l);
                //local_converged = true;
                is_converged    = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(),rDofSet,mA, mDx, e);
                is_converged    = bool((is_converged == true && local_converged == true) && mIterationNumber > 2);

                //mpConvergenceCriteria->GetConvergenceData(Parameters);
                //new_residual = Parameters[0];
                //KRATOS_WATCH(new_residual)
                //KRATOS_WATCH(old_residual)

                //if(new_residual>old_residual)
                //{
                //    std::cout << "Calling Line Searches " << std::endl; // In a line search, mDx is the total acumulated in that step
                //    Searches.mso = TSparseSpace::Dot(mDelta_p, mb);
                //    KRATOS_WATCH(Searches.mso)
                //    Searches.LineSearches(BaseType::GetModelPart(), pScheme, pBuilderAndSolver, rDofSet, mX_old, mDelta_p, mDx, mb, mA);
                //    mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(),rDofSet,mA, mDx, mb);
                //    mpConvergenceCriteria->GetConvergenceData(Parameters);
                //    new_residual = Parameters[0];
                //}

//                 old_residual = new_residual;

                if(is_converged == true)
                {
                    mReduceArcLenght = false;
                }

                // Optional recomputing K(mA y sigma_q )
                if(is_converged == false)
                {
                    TSparseSpace::SetToZero(mDx);
                    TSparseSpace::SetToZero(mb);
                    TSparseSpace::SetToZero(mA);
                    pBuilderAndSolver->Build(pScheme,mAuxElementModelPart,mA,mb);
                    pBuilderAndSolver->SystemSolve(mA, Sigma_q, q_Inc_Aux);
                    //noalias(Sigma_q) += mDelta_pold; /// Should be the total acumulated
                    TSparseSpace::Copy(Aux_q, q_Inc_Aux); /// q = Aux_q
                }

                // Optional recomputing K(mA y sigma_q )
                if(is_converged==false)
                {
                    TSparseSpace::SetToZero(mDx);
                    TSparseSpace::SetToZero(mb);
                    TSparseSpace::SetToZero(mA);
                    if (this->GetEchoLevel() > 0)
                    {
                        std::cout << " Rebuilding Tangent Matrix " << std::endl;
                    }
                    pBuilderAndSolver->Build(pScheme,mAuxElementModelPart,mA,mb);
                    TSparseSpace::Copy(Aux_q, q_Inc_Aux); /// q = Aux_q
                    pBuilderAndSolver->SystemSolve(mA, Sigma_q, q_Inc_Aux);
                    //noalias(Sigma_q) += mDelta_pold; /// Should be the total acumulated
                    TSparseSpace::Copy(Aux_q, q_Inc_Aux); /// q = Aux_q
                }

                if (is_converged == false && mIterationNumber >= mMaxIterationNumber)
                {
                    mReduceArcLenght = true;
                    mdelta_l          = std::sqrt(RealType(mIde)/RealType(mIterationNumber)) * mdelta_l;
                    //mdelta_l          = (recursive + 1) * mdelta_lold; // Increasing the arc-length
                    mdelta_lambda     = mdelta_lambda_old;
                    meta              = 1.00;
                    mlambda           = mlambda_old;

                    std::cout << "***************************************************" << std::endl;
                    std::cout << "******* ATTENTION: Max Iterations Exceeded ********" << std::endl;
                    std::cout << "***************************************************" << std::endl;

                    TSparseSpace::SetToZero(mDx);
                    TSparseSpace::SetToZero(mb);
                    TSparseSpace::SetToZero(mA);
                    SetDatabaseToValue(rDofSet, mX_old);
                    pBuilderAndSolver->Build(pScheme,mAuxElementModelPart,mA,mb);

                    // Solve again
                    TSparseSpace::Copy(Aux_q, q_Inc_Aux); /// q = Aux_q
                    pBuilderAndSolver->SystemSolve(mA, Sigma_q, q_Inc_Aux);
                    TSparseSpace::Copy(Aux_q, q_Inc_Aux); /// q = Aux_q

                    if (this->GetEchoLevel() > 0)
                    {
                        std::cout << "Arc lenght modified  = " << mdelta_l  <<std::endl;
                    }
                }
            }
    } // end while

    while(mReduceArcLenght == true);

    if(is_converged == true)
    {
        // Calculate reactions if required
        if (mCalculateReactionsFlag ==true)
        {
            pBuilderAndSolver->CalculateReactions(pScheme,BaseType::GetModelPart(),mA,mDx,mb);
        }

        // Finalisation of the solution step, operations to be done after achieving convergence, for example the
        // Final Residual Vector (mb) has to be saved in there
        FinalizeSolutionStep();

        // Cleaning memory after the solution
        pScheme->Clean();
        mSolutionStepIsInitialized = false;

        // Deallocate the systemvectors
        if (mReformDofSetAtEachStep == true)
        {
            Clear();
        }
        else
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "*********** WARNING: Step no converged ************" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }
    }

    return 0.00;

    KRATOS_CATCH("");

    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    *  This should be considered as a "post solution" convergence check which is useful for coupled analysis
    *  the convergence criteria used is the one used inside the "solve" step
    */

    bool IsConverged() override
    {
        KRATOS_TRY;

        TSystemMatrixType& mA  = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb  = *mpb;

        if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
        {
            GetBuilderAndSolver()->BuildRHS(GetScheme(),BaseType::GetModelPart(),mb);
        }

        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();

        return mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

        KRATOS_CATCH("");

    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * This operations should be called before printing the results when non trivial results (e.g. stresses)
    * need to be calculated given the solution of the step
    * This operations should be called only when needed, before printing as it can involve a non negligible cost
    */

    void CalculateOutputData() override
    {
        TSystemMatrixType& mA  = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb  = *mpb;

        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();
        GetScheme()->CalculateOutputData(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It clears the variables of the arc length
    */

    void Clear() override
    {
        KRATOS_TRY;
        if (this->GetEchoLevel() > 0)
        {
            std::cout << "Arc Length Strategy  Clear function used" << std::endl;
        }

        TSystemMatrixType& mA          = *mpA;
        TSystemVectorType& mDx         = *mpDx;
        TSystemVectorType& mb          = *mpb;
        TSystemVectorType& mDelta_p    = *mpDelta_p;
        TSystemVectorType& mDelta_pold = *mpDelta_pold;

        SparseSpaceType::Clear(mpA);
        SparseSpaceType::Resize(mA, 0, 0);

        SparseSpaceType::Clear(mpDx);
        SparseSpaceType::Resize(mDx, 0);

        SparseSpaceType::Clear(mpb);
        SparseSpaceType::Resize(mb, 0);

        SparseSpaceType::Clear(mpDelta_p);
        SparseSpaceType::Resize(mDelta_p, 0);

        SparseSpaceType::Clear(mpDelta_pold);
        SparseSpaceType::Resize(mDelta_pold, 0);

        // Setting to zero the internal flag to ensure that the dof sets are recalculated
        GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        GetBuilderAndSolver()->Clear();

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Recursive function inside the arc length that computes the new delta_lambda
    * mDx are the increments, dx = dx_old + mDx()
    * mDelta_p is the accumulated
    * @param pq:
    * @param pSigma_q:
    * @param pSigma_h:
    * @param pdx_aux:
    * @param g:
    * @param Ao:
    * @return mdelta_lambda: The increment in the load factor
    */

    void Recursive_Function_Arc_Length(
            TSystemVectorPointerType& pq,
            TSystemVectorPointerType& pSigma_q,
            TSystemVectorPointerType& pSigma_h,
            TSystemVectorPointerType& pdx_aux,
            RealType& g,
            RealType& Ao
            )
    {
        KRATOS_TRY;

        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        DofsArrayType& rDofSet                                    = pBuilderAndSolver->GetDofSet();
        typename TSchemeType::Pointer pScheme                     = GetScheme();

        // Coeficients to solve the  2nd equation.
        RealType a    = 0.00;
        RealType b    = 0.00;
        RealType c    = 0.00;
        RealType disc = 0.00;           // Discriminant of the quadratic equation
        RealType x    = 0.00;
        std::vector<RealType> x_sol(2);              // Solution of the second order equation
        x_sol.resize(2,false);

        RealType lambda_cr       = 0.00;
        RealType delta_lambda_cr = 0.00;
        RealType delta_lcr       = 0.00;
        RealType miu             = 0.00;
        bool  imag = false;

        // Aux Variables
        TSystemVectorPointerType pAux_Vector;
        TSystemVectorPointerType pDelta_p;
        TSystemVectorPointerType pDelta_p1; // For first roots
        TSystemVectorPointerType pDelta_p2; // For second roots

        InitializeAuxVectors(pAux_Vector);
        InitializeAuxVectors(pDelta_p);
        InitializeAuxVectors(pDelta_p1);
        InitializeAuxVectors(pDelta_p2);

        TSystemVectorType& Aux_Vector  = *pAux_Vector;
        TSystemVectorType& Delta_p     = *pDelta_p;
        TSystemVectorType& Delta_p1    = *pDelta_p1;
        TSystemVectorType& Delta_p2    = *pDelta_p2;
        TSystemVectorType& dx_aux      = *pdx_aux;

        TSparseSpace::SetToZero(Aux_Vector);
        TSparseSpace::SetToZero(Delta_p);
        TSparseSpace::SetToZero(Delta_p1);
        TSparseSpace::SetToZero(Delta_p2);

        // Variables vectoriales y matriciles
        TSystemMatrixType& mA            = *mpA;
        //TSystemVectorType& mDx           = *mpDx;
        TSystemVectorType& mb            = *mpb;
        TSystemVectorType& mDelta_p      = *mpDelta_p;
        TSystemVectorType& mDelta_pold   = *mpDelta_pold;
        TSystemVectorType& Sigma_q       = *pSigma_q;
        TSystemVectorType& Sigma_h       = *pSigma_h;
        TSystemVectorType& mX_old        = *mpX_old;
        TSystemVectorType& q             = *pq;

        // Constants needed for the Ublas operations
        RealType A = 1.00;
        //RealType B = 1.00;

        // Calculate_Current_Delta(rDofSet, Delta_p);

        a = Ao + TSparseSpace::Dot(Sigma_q, Sigma_q);
        TSparseSpace::ScaleAndAdd(A, mDelta_p, meta, Sigma_h, Aux_Vector); // Aux_Vector = A * mDelta_p + meta*Sigma_h
        b = 2.00 * (Ao * (mdelta_lambda-g) + TSparseSpace::Dot(Sigma_q, Aux_Vector));
        c = Ao * (mdelta_lambda-g) * (mdelta_lambda-g) - mdelta_l * mdelta_l + TSparseSpace::Dot(Aux_Vector, Aux_Vector);

        //KRATOS_WATCH(Ao);
        //KRATOS_WATCH(meta);
        //KRATOS_WATCH(Sigma_q);
        //KRATOS_WATCH(mDelta_p);
        //KRATOS_WATCH(Sigma_h);
        //KRATOS_WATCH(g);
        //KRATOS_WATCH(mdelta_lambda);
        //KRATOS_WATCH(Aux_Vector);
        //KRATOS_WATCH(mdelta_l);
        //KRATOS_WATCH(a);
        //KRATOS_WATCH(b);
        //KRATOS_WATCH(c);

        disc = b * b - 4.00 * a * c;
        if (disc >= 0.00)
        {
            StructuralMechanicsMathUtilities::SolveSecondOrderEquation(a,b,c,x_sol);

            TSparseSpace::ScaleAndAdd(x_sol[0],Sigma_q,meta,Sigma_h,Delta_p1); //Delta_p1 = x_sol(0)*Sigma_q + meta*Sigma_h
            TSparseSpace::ScaleAndAdd(x_sol[1],Sigma_q,meta,Sigma_h,Delta_p2); //Delta_p2 = x_sol(1)*Sigma_q + meta*Sigma_h

            if (this->GetEchoLevel() > 1)
            {
                std::cout<<" Real roots found " << std::endl;
                std::cout<<" First Solution  = " << x_sol[0] <<  std::endl;
                std::cout<<" Second Solution = " << x_sol[1] <<  std::endl;
            }

            // Choose the x value: the larges dot product
            // WARNING: The old code use the current incremental displacement
            // First roots
            noalias(Delta_p) = mDelta_p + Delta_p1;
            RealType a1        = TSparseSpace::Dot(Delta_p, mDelta_pold);

            //KRATOS_WATCH(Delta_p1[0]);
            //KRATOS_WATCH(Delta_p2[0]);
            //KRATOS_WATCH(mDelta_pold);
            //KRATOS_WATCH(a1);

            // Second roots
            TSparseSpace::SetToZero(Delta_p);
            noalias(Delta_p) = mDelta_p + Delta_p2;
            RealType a2        = TSparseSpace::Dot(Delta_p,mDelta_pold);
            //KRATOS_WATCH(a2);

            if(a1 > a2)
            {
                x = x_sol[0];
                noalias(mDelta_p)+= Delta_p1;
                TSparseSpace::Copy(Delta_p1, dx_aux);
            }
            else
            {
              x = x_sol[1];
              noalias(mDelta_p)+= Delta_p2;
              TSparseSpace::Copy(Delta_p2, dx_aux);
            }

            mdelta_lambda += - g + x;

            if (this->GetEchoLevel() > 1)
            {
                std::cout << " Solution Chosen = " << x <<  std::endl;
                std::cout << " New DeltaLamda  = " << mdelta_lambda <<  std::endl;
            }
        }
        else
        {
            std::cout << "WARNING: No real roots were found " << std::endl;
            std::cout << "Introducting a pseudo-line search to avoid complex roots " << std::endl;
            if (this->GetEchoLevel() > 1)
            {
                std::cout << "Calculating eta" << std::endl;
                Calculate_eta(Ao, pSigma_q, pSigma_h, g, imag);
            }

            if (imag == false)
            {
                if (this->GetEchoLevel() > 1)
                {
                    std::cout << " eta was found with value = " << meta << std::endl;
                    std::cout << " Calling again the Recursive Function " << std::endl;
                }
                Recursive_Function_Arc_Length(pq, pSigma_q, pSigma_h, pdx_aux, g, Ao);
            }
            else
            {
                std::cout << "WARNING: No real roots were found. Avoid " << std::endl;
                TSystemVectorPointerType pDelta_pcr;
                InitializeAuxVectors(pDelta_pcr);
                TSystemVectorType& Delta_pcr = *pDelta_pcr;
                TSparseSpace::SetToZero(Delta_pcr);
                noalias(Delta_pcr) = mDelta_p + Sigma_h;

                //this->BackupDatabase(rDofSet,mX_old);

                TSparseSpace::Copy(Sigma_h, dx_aux);

                // Update results
                pScheme->Update(BaseType::GetModelPart(),rDofSet,mA, dx_aux ,mb);
                if(this->MoveMeshFlag() == true)
                {
                    BaseType::MoveMesh();
                }

                TSparseSpace::SetToZero(mb);

                pBuilderAndSolver->BuildRHS(pScheme,mAuxElementModelPart,mb);
                lambda_cr          = -TSparseSpace::Dot(mb,q)/TSparseSpace::Dot(q,q);
                delta_lambda_cr    = lambda_cr - mlambda_old;
                delta_lcr          = std::sqrt(TSparseSpace::Dot(Delta_pcr,Delta_pcr) + Ao*(delta_lambda_cr)*(delta_lambda_cr));
                miu                = mdelta_l/delta_lcr;

                this->BackupDatabase(rDofSet,mX_old);
                TSparseSpace::Copy(mDelta_p, dx_aux);

                noalias(mDelta_p) = miu * Delta_pcr;
                mdelta_lambda     = miu * delta_lambda_cr;
                mdelta_l          = delta_lcr;

                if (this->GetEchoLevel() > 1)
                {
                    std::cout << "   Arc Length      = " << mdelta_l             << std::endl;
                    std::cout << "   Lamba Old       = " << mlambda_old          << std::endl;
                    std::cout << "   Delta Lamba     = " << mdelta_lambda        << std::endl;
                    std::cout << "   Lamda_cr        = " << lambda_cr            << std::endl;
                    std::cout << "   Delta_Lamda_cr  = " << delta_lambda_cr      << std::endl;
                    std::cout << "   Miu Factor      = " << miu                  << std::endl;
                    std::cout << "   Eta             = " << meta                 << std::endl;
                }

                //TSparseSpace::SetToZero(mDx);
                //TSparseSpace::Copy(mDelta_p, mDx);
            }
        }

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Returns the LHS of the problem
    */

    TSystemMatrixType& GetSystemMatrix()
    {
        TSystemMatrixType& mA = *mpA;

        return mA;
    }

    /***********************************************************************************/
    /***********************************************************************************/

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
    TSystemVectorPointerType mpX_old;
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
    bool mInit;

    bool mSolutionStepIsInitialized;
    unsigned int mMaxIterationNumber;
    unsigned int mstep;
    unsigned int mMaxRecursive;
    unsigned int mIterationNumber;
    bool mReduceArcLenght;

    RealType mdelta_l;         // Arc length
    RealType mdelta_lold;      // Arc length from the previous increment
    RealType mdelta_lmax;      // Maximum arc length allowed
    RealType mfactor_delta_lmax;
    RealType meta;
    int mIde;
    RealType mlambda;
    RealType mlambda_old;
    RealType mdelta_lambda;
    RealType mdelta_lambda_old;
    ModelPart& mAuxElementModelPart;
    ModelPart& mAuxConditionModelPart;

    /*@} */
    /**@name Private Operators*/
    /*@{ */

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Initilise the variables, schemes and convergence criterias
    */

    void Initialize() override
    {
        KRATOS_TRY;

        // Pointers needed in the solution
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;

        // Initialize The Scheme - OPERATIONS TO BE DONE ONCE
        if (pScheme->SchemeIsInitialized() == false)
        {
            pScheme->Initialize(BaseType::GetModelPart());
        }

        // Initialize The Elements - OPERATIONS TO BE DONE ONCE
        if (pScheme->ElementsAreInitialized() == false)
        {
            pScheme->InitializeElements(BaseType::GetModelPart());
        }

        // Initialize Conditions
        if (pScheme->ConditionsAreInitialized() == false)
        {
            pScheme->InitializeConditions(mAuxConditionModelPart);
        }

        // Initialisation of the convergence criteria
        if (mpConvergenceCriteria->IsInitialized() == false)
        {
            pConvergenceCriteria->Initialize(BaseType::GetModelPart());
        }

        mInitializeWasPerformed = true;

        VariablesArcLength(); // Initialising the variables of the arc length

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It initialises the solution step
    */

    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        if (this->GetEchoLevel() > 0)
        {
            std::cout<< "Initializing Solution Step " << std::endl;
        }
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;
        DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();

        // Setting up the Vectors involved to the correct size with value cero
        pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, mpA,mpDx,mpb,BaseType::GetModelPart());
        InitializeAuxVectors(mpDelta_p);
        InitializeAuxVectors(mpDelta_pold);
        InitializeAuxVectors(mpX_old);

        TSystemMatrixType& mA            = *mpA;          // Stiffness Matrix
        TSystemVectorType& mDx           = *mpDx;         // External Force
        TSystemVectorType& mb            = *mpb;          // Internal Force
        TSystemVectorType& mDelta_p      = *mpDelta_p;    // P  current change
        TSystemVectorType& mDelta_pold   = *mpDelta_pold; // P  =  u_(step+1)-u_(step)
        TSystemVectorType& mX_old        = *mpX_old;      // old = positions X+u

        TSparseSpace::SetToZero(mDelta_p);
        TSparseSpace::SetToZero(mDelta_pold);
        TSparseSpace::SetToZero(mX_old);

        Calculate_Previous_Delta(rDofSet, mDelta_pold); // Store the last converged delta P
        BackupDatabase(rDofSet, mX_old);                // Store the actual point x = X + u_n = X_0 + Deltap_1 + Deltap_2....+ Deltap_n
        meta = 1.00;                                    // Reseting meta = 1.00; always we begining with 1.00

        // Initial operations
        pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
        pScheme->InitializeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
        pConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(),rDofSet, mA, mDx, mb);

        mSolutionStepIsInitialized = true;

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It finalises the arc length for the currrent step
    * @param mIterationNumber: The iteration number in the non-linear step
    * @param mReduceArcLenght: Boolean that tells if the arc length has been computed with the reduced method
    */

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        typename TBuilderAndSolverType::Pointer pBuilderAndSolver       = GetBuilderAndSolver();
        typename TSchemeType::Pointer pScheme                           = GetScheme();
        typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;
        DofsArrayType& rDofSet                                          = GetBuilderAndSolver()->GetDofSet();

        TSystemMatrixType& mA            = *mpA;
        TSystemVectorType& mDx           = *mpDx;
        TSystemVectorType& mb            = *mpb;

        RealType factor     = 1.00;
        mdelta_lambda_old =  mdelta_lambda;
        mlambda_old       =  mlambda;

        // KRATOS_WATCH(mlambda_old)

        factor           = std::sqrt(RealType(mIde)/RealType(mIterationNumber));

        // Controling the size of the arc
        if (factor > 1.0)
        {
            factor = 1.00;
        }
        if (factor < 0.75)
        {
            factor = 0.75;
        }

        mdelta_l = factor * mdelta_lold;
        if (mdelta_lold > mdelta_lmax)
        {
            mdelta_lold = mdelta_lmax;
            mdelta_l    = mdelta_lmax;
        }

        mReduceArcLenght = false;

        pScheme->FinalizeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
        pBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
        pConvergenceCriteria->FinalizeSolutionStep(BaseType::GetModelPart(),rDofSet, mA, mDx, mb);

        BaseType::GetModelPart().GetProcessInfo()[LAMBDA] = mlambda;

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It initialises an auxiliar model part
    * @param ThisModelPart: Model part
    */

    void InitializeAuxiliaryModelParts(ModelPart& ThisModelPart)
    {
        mAuxElementModelPart.SetBufferSize(ThisModelPart.GetBufferSize());
        mAuxConditionModelPart.SetBufferSize(ThisModelPart.GetBufferSize());

        mAuxElementModelPart.Nodes() =              ThisModelPart.Nodes();
        mAuxConditionModelPart.Nodes()    =        ThisModelPart.Nodes();
        mAuxElementModelPart.PropertiesArray()   = ThisModelPart.PropertiesArray();
        mAuxConditionModelPart.PropertiesArray() = ThisModelPart.PropertiesArray();
        mAuxElementModelPart.Elements() = ThisModelPart.Elements();
        mAuxConditionModelPart.Conditions() = ThisModelPart.Conditions();
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It calculates the parameter eta
    * @param Ao:
    * @param pSigma_q:
    * @param pSigma_h:
    * @param g:
    * @param imag:
    * @return meta: The eta parameter
    */

    void Calculate_eta(
            RealType& Ao,
            TSystemVectorPointerType& pSigma_q,
            TSystemVectorPointerType& pSigma_h,
            RealType& g,
            bool& imag
            )
    {
        RealType a_prima = 0.00;
        RealType b_prima = 0.00;
        RealType c_prima = 0.00;
        RealType disc    = 0.00;
        std::vector<RealType> solution;
        solution.resize(2, false);

        TSystemVectorType& Sigma_q     = *pSigma_q;
        TSystemVectorType& Sigma_h     = *pSigma_h;
        TSystemVectorType& mDelta_p    = *mpDelta_p;
        //KRATOS_WATCH(Sigma_q);
        //KRATOS_WATCH(Sigma_h);

        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
//        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();

        // Necessary to find the roots. mDelta_p becomes aDx from the previous iteration
        /* for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        if(i_dof->IsFree())
            mDelta_p[i_dof->EquationId()] = i_dof->GetSolutionStepValue(0)-i_dof->GetSolutionStepValue(1);
        */

        /// WARNING: Verify current mDelta_p
        //Calculate_Delta_pold(rDofSet,mDelta_p);

        RealType param_a = TSparseSpace::Dot(Sigma_q,  Sigma_q);
        RealType param_b = TSparseSpace::Dot(Sigma_h,  Sigma_h);
        RealType param_c = TSparseSpace::Dot(Sigma_q,  Sigma_h);
        RealType param_d = TSparseSpace::Dot(mDelta_p, Sigma_h);
        RealType param_e = TSparseSpace::Dot(mDelta_p, mDelta_p);
        RealType param_f = TSparseSpace::Dot(Sigma_q,  mDelta_p);

        a_prima = (Ao + param_a) * param_b - param_c * param_c;
        b_prima = 2.00 * ( (Ao + param_a)*param_d - ((Ao * (mdelta_lambda-g) + param_f)) * param_c);
        c_prima = (Ao + param_a) * (param_e - mdelta_l*mdelta_l) - (2.00 * Ao * (mdelta_lambda - g) + param_f) * param_f + param_a * Ao * (mdelta_lambda-g) * (mdelta_lambda-g);

        // Before
        //a_prima = (Ao + inner_prod(Sigma_q,Sigma_q))*(inner_prod(Sigma_h,Sigma_h))-(inner_prod(Sigma_q,Sigma_h))*(inner_prod(Sigma_q,Sigma_h));
        //b_prima = 2.00*((Ao + inner_prod(Sigma_q,Sigma_q))*(inner_prod(mDelta_p,Sigma_h))-((Ao*(mdelta_lambda-g)+ inner_prod(Sigma_q,mDelta_p)))*(inner_prod(Sigma_q,Sigma_h)));
        //c_prima = (Ao + inner_prod(Sigma_q,Sigma_q))*((inner_prod(mDelta_p,mDelta_p)-mdelta_l*mdelta_l))-(2.00*Ao*(mdelta_lambda-g)+inner_prod(Sigma_q,mDelta_p))*inner_prod(Sigma_q,mDelta_p) + inner_prod(Sigma_q,Sigma_q)*Ao*(mdelta_lambda-g)*(mdelta_lambda-g);

        disc    = b_prima * b_prima - 4.00 * a_prima * c_prima;
        //KRATOS_WATCH(a_prima);
        //KRATOS_WATCH(b_prima);
        //KRATOS_WATCH(c_prima);
        //KRATOS_WATCH(disc);

        if (disc >= 0.00)
        {
            imag = false;
            StructuralMechanicsMathUtilities::SolveSecondOrderEquation(a_prima, b_prima, c_prima, solution);

            if(solution[0] < 0.00)
            {
                solution[0] = 0.00;
            }

            if(solution[1] < 0.00)
            {
                solution[1] = 0.00;
            }

            //KRATOS_WATCH(solution);
            //KRATOS_WATCH(a_prima);
            //KRATOS_WATCH(b_prima);
            //KRATOS_WATCH(c_prima);

            if (solution[0] > solution[1])
            {
                RealType a =   solution[1];
                solution[1]     =   solution[0];
                solution[0]     =   a;
            }

            RealType shi  = 0.05 * std::abs((solution[1] - solution[0]));
            if (solution[1] < 1.00)
            {
                meta  = solution[1] - shi;
            }
            else if ((1.00 > (-b_prima/a_prima)) && (solution[1] > 1.00))
            {
                meta  = solution[1] + shi;
            }
            else if ((1.00 < (-b_prima/a_prima)) && (solution[0] < 1.00))
            {
                meta  = solution[0] -shi;
            }
            else
            {
                meta  = solution[0] + shi;
            }
        }
        else
        {
            std::cout << "eta is an imaginary number" << std::endl;;
            meta   = 1.00;
            imag  = true;
        }
    }

      /***********************************************************************************/
      /***********************************************************************************/

      /**
      * Computed the increment of displacements from the last converged point to the previous step  in the iteraction  i+1
      * @param rDofSet: Set of degrees of freedom
      * @return: Delta_pold: Increment of displacements
      */

      void Calculate_Previous_Delta(
              DofsArrayType const & rDofSet,
              TSystemVectorType& Delta_pold
              )
      {
        KRATOS_TRY;

        for(typename DofsArrayType::const_iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; i_dof ++)
        {
            if(i_dof->IsFree())
            {
                Delta_pold[i_dof->EquationId()] = i_dof->GetSolutionStepValue(1) - i_dof->GetSolutionStepValue(2);
            }
        }

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Computed the increment of displacements from the last converged point to the current step  in the iteraction  i+1
    * @param rDofSet: Set of degrees of freedom
    * @return: Delta_pold: Increment of displacements
    */

    void Calculate_Current_Delta(
        DofsArrayType const & rDofSet,
        TSystemVectorType& Delta_pold
    )
    {
        KRATOS_TRY;

        for(typename DofsArrayType::const_iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
	{
            if(i_dof->IsFree())
	    {
                Delta_pold[i_dof->EquationId()] = i_dof->GetSolutionStepValue(0) - i_dof->GetSolutionStepValue(1);
            }
	}

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Set the values of the database to the corresponding to X_old
    * @return rDofSet: Set of degrees of freedom
    * @param X_old: The old displacements
    */

    void SetDatabaseToValue(
      DofsArrayType& rDofSet,
      const TSystemVectorType& X_old
    )
    {
        KRATOS_TRY;

        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree())
            {
                i_dof->GetSolutionStepValue() = X_old[i_dof->EquationId()];
            }
        }

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It allows to write the old displacements in X_old
    * @param rDofSet: Set of degrees of freedom
    * @return X_old: The old displacements
    */
    void BackupDatabase(
      DofsArrayType const& rDofSet,
      TSystemVectorType& X_old
    )
    {
        KRATOS_TRY;

        for(typename DofsArrayType::const_iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree())
            {
                X_old[i_dof->EquationId()] = i_dof->GetSolutionStepValue();
            }
        }

        KRATOS_CATCH("");
    }

    ResidualBasedArcLengthStrategy(const ResidualBasedArcLengthStrategy& Other);

    /*@} */

}; /* Class ResidualBasedArcLengthStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_ARC_LENGTH_STRATEGY  defined */

