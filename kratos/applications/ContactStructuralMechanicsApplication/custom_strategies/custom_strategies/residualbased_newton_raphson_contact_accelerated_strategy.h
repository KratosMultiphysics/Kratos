// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: structural_mechanics_application/license.txt
//
//  Main authors:    Mohamed Khalil
//                   Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_RESIDUALBASED_NEWTON_RAPHSON_CONTACT_ACCELERATED_STRATEGY)
#define KRATOS_RESIDUALBASED_NEWTON_RAPHSON_CONTACT_ACCELERATED_STRATEGY

/* System Includes */

/* External Includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "custom_utilities/contact_utilities.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"

// Default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

// Convergence accelerators
#include "../FSIapplication/custom_utilities/convergence_accelerator.hpp"

// TODO: Try to derive it directly from the one in the core
// TODO: Extend the descriptions

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
    
/** \brief  Short class definition.
This class 
*/

template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
         
class ResidualBasedNewtonRaphsonContactAcceleratedStrategy :
    public SolvingStrategy< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;
    
    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedNewtonRaphsonContactAcceleratedStrategy );

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;
    
    typedef ConvergenceAccelerator<TDenseSpace> TConvergenceAcceleratorType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    
    /**
     * Default constructor 
     * @param rModelPart: The model part of the problem
     * @param pScheme: The integration scheme
     * @param pNewLinearSolver: The linear solver employed
     * @param pNewConvergenceCriteria: The convergence criteria employed
     * @param MaxIterationNumber: The maximum number of iterations
     * @param CalculateReactions: The flag for the reaction calculation
     * @param ReformDofSetAtEachStep: The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag: The flag that allows to move the mesh
     * @param ConvergenceAccelerator: The convergence accelerator to use
     */
    
    ResidualBasedNewtonRaphsonContactAcceleratedStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        unsigned int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false,
        typename TConvergenceAcceleratorType::Pointer pConvergenceAccelerator = nullptr,
        unsigned int MaxNumberConvergenceAccelerationIterations = 20,
        double ReductionCoefficient = 0.1,
        typename TConvergenceCriteriaType::Pointer pConvergenceCriteriaForConvergenceAccelerator = nullptr,
        bool UpdateSystem = false
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, MoveMeshFlag)
    {
        KRATOS_TRY

        mKeepSystemConstantDuringIterations = false;

        // Set flags to default values
        SetMaxIterationNumber(MaxIterations);
        mCalculateReactionsFlag = CalculateReactions;
        mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        //Saving the convergence criteria to be used
        mpConvergenceCriteria = pNewConvergenceCriteria;

        // Saving the scheme
        mpScheme = pScheme;
        
        // Saving convergence accelerator
        mpConvergenceAccelerator = pConvergenceAccelerator;
        
        // Setting convergence acceleration parameters 
        mMaxNumberConvergenceAccelerationIterations = MaxNumberConvergenceAccelerationIterations;
        mReductionCoefficient = ReductionCoefficient;
        mpConvergenceCriteriaForConvergenceAccelerator = pConvergenceCriteriaForConvergenceAccelerator;
        mUpdateSystem = UpdateSystem;

        // Saving the linear solver
        mpLinearSolver = pNewLinearSolver;

        // Setting up the default builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer
                             (
                                 new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > (mpLinearSolver)
                             );

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        //By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("");
    }

    /**
     * Default constructor 
     * @param rModelPart: The model part of the problem
     * @param pScheme: The integration scheme
     * @param pNewLinearSolver: The linear solver employed
     * @param pNewConvergenceCriteria: The convergence criteria employed
     * @param MaxIterationNumber: The maximum number of iterations
     * @param CalculateReactions: The flag for the reaction calculation
     * @param ReformDofSetAtEachStep: The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag: The flag that allows to move the mesh
     * @param ConvergenceAccelerator: The convergence accelerator to use
     */
    
    ResidualBasedNewtonRaphsonContactAcceleratedStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        unsigned int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false,
        typename TConvergenceAcceleratorType::Pointer pConvergenceAccelerator = nullptr,
        unsigned int MaxNumberConvergenceAccelerationIterations = 20,
        double ReductionCoefficient = 0.1,
        typename TConvergenceCriteriaType::Pointer pConvergenceCriteriaForConvergenceAccelerator = nullptr,
        bool UpdateSystem = false
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, MoveMeshFlag)
    {
        KRATOS_TRY;

        mKeepSystemConstantDuringIterations = false;

        // Set flags to default values
        SetMaxIterationNumber(MaxIterations);
        mCalculateReactionsFlag = CalculateReactions;
        mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        // Saving the convergence criteria to be used
        mpConvergenceCriteria = pNewConvergenceCriteria;

        // Saving the scheme
        mpScheme = pScheme;

        // Saving convergence accelerator
        mpConvergenceAccelerator = pConvergenceAccelerator;
        
        // Setting convergence acceleration parameters 
        mMaxNumberConvergenceAccelerationIterations = MaxNumberConvergenceAccelerationIterations;
        mReductionCoefficient = ReductionCoefficient;
        mpConvergenceCriteriaForConvergenceAccelerator = pConvergenceCriteriaForConvergenceAccelerator;
        mUpdateSystem = UpdateSystem;
        
        // Saving the linear solver
        mpLinearSolver = pNewLinearSolver;

        // Setting up the default builder and solver
        mpBuilderAndSolver = pNewBuilderAndSolver;

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("");
    }

    /** 
     * Destructor.
     */
    virtual ~ResidualBasedNewtonRaphsonContactAcceleratedStrategy()
    {
    }
    
    // Set and Get Scheme ... containing Builder, Update and other

    void SetScheme(typename TSchemeType::Pointer pScheme)
    {
        mpScheme = pScheme;
    };

    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };
    void SetConvergenceAccelerator(typename TConvergenceAcceleratorType::Pointer pConvergenceAccelerator)
    {
        mpConvergenceAccelerator = pConvergenceAccelerator;
    };

    typename TConvergenceAcceleratorType::Pointer GetConvergenceAccelerator()
    {
        return mpConvergenceAccelerator;
    };

    // Set and Get the BuilderAndSolver

    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    // Set and Get the flag for calculating the reactions 
    
    void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
    {
        mCalculateReactionsFlag = CalculateReactionsFlag;
    }

    bool GetCalculateReactionsFlag()
    {
        return mCalculateReactionsFlag;
    }

    // Set and get for the flag to reform the DOF at each step
    
    void SetReformDofSetAtEachStepFlag(bool flag)
    {
        mReformDofSetAtEachStep = flag;
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
    }

    bool GetReformDofSetAtEachStepFlag()
    {
        return mReformDofSetAtEachStep;
    }

    // Set and get of the maximum iteration number 
    
    void SetMaxIterationNumber(unsigned int MaxIterationNumber)
    {
        mMaxIterationNumber = MaxIterationNumber;
    }

    unsigned int GetMaxIterationNumber()
    {
        return mMaxIterationNumber;
    }
    
    /**
     * It levels the echo for the solving strategy
     * @param Level:  0 -> mute... no echo at all
     *                1 -> printing time and basic informations
     *                2 -> printing linear solver data
     *                3 -> Print of debug informations:
     *                     Echo of stiffness matrix, Dx, b...
     */
        
    void SetEchoLevel(int Level)
    {
        BaseType::mEchoLevel = Level;
        GetBuilderAndSolver()->SetEchoLevel(Level);
    }
    
    /******************** OPERATIONS ACCESSIBLE FROM THE INPUT: ************************/
    /***********************************************************************************/

    /**
     * Operation to predict the solution ... if it is not called a trivial predictor is used in which the
     * values of the solution step of interest are assumed equal to the old values
     */
    
    void Predict()
    {
        KRATOS_TRY;
        
        // OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        // If the operations needed were already performed this does nothing
        if (mInitializeWasPerformed == false)
        {
            Initialize();
        }

        // Initialize solution step
        if (mSolutionStepIsInitialized == false)
        {
            InitializeSolutionStep();
        }

        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        GetScheme()->Predict(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        // Move the mesh if needed
        if (this->MoveMeshFlag() == true) 
        {
            BaseType::MoveMesh();
        }
        
        KRATOS_CATCH("");
    }

    /**
     * Initialization of member variables and prior operations
     */
    void Initialize()
    {
        KRATOS_TRY;

        if (mInitializeWasPerformed == false)
        {
            // Pointers needed in the solution
            typename TSchemeType::Pointer pScheme = GetScheme();
            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;

            //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            if (pScheme->SchemeIsInitialized() == false)
            {
                pScheme->Initialize(BaseType::GetModelPart());
            }

            //Initialize The Elements - OPERATIONS TO BE DONE ONCE
            if (pScheme->ElementsAreInitialized() == false)
            {
                pScheme->InitializeElements(BaseType::GetModelPart());
            }

            //Initialize The Conditions - OPERATIONS TO BE DONE ONCE
            if (pScheme->ConditionsAreInitialized() == false)
            {
                pScheme->InitializeConditions(BaseType::GetModelPart());
            }

            // Initialisation of the convergence criteria
            if (mpConvergenceCriteria->mConvergenceCriteriaIsInitialized == false)
            {
                mpConvergenceCriteria->Initialize(BaseType::GetModelPart());
            }
            
            // Initialisation of the convergence criteria for the convergence accelerator
            if (mpConvergenceCriteriaForConvergenceAccelerator->mConvergenceCriteriaIsInitialized == false)
            {
                mpConvergenceCriteriaForConvergenceAccelerator->Initialize(BaseType::GetModelPart());
            }
            
            mpConvergenceAccelerator->Initialize(); // TODO: Check if correct

            mInitializeWasPerformed = true;
        }

        KRATOS_CATCH("");
    }

    /**
     * The problem of interest is solved.
     * This function calls sequentially: Initialize(), InitializeSolutionStep(), Predict(), SolveSolutionStep() and FinalizeSolutionStep().
     * All those functions can otherwise be called separately.
     */
    
    double Solve()
    {
        Initialize();
        InitializeSolutionStep();
        Predict();
        SolveSolutionStep();
        FinalizeSolutionStep();
        return 0.00;
    }

    /**
     * Clears the internal storage
     */
    
    void Clear()
    {
        KRATOS_TRY;

        SparseSpaceType::Clear(mpA);
        TSystemMatrixType& mA = *mpA;
        SparseSpaceType::Resize(mA, 0, 0);

        SparseSpaceType::Clear(mpDx);
        TSystemVectorType& mDx = *mpDx;
        SparseSpaceType::Resize(mDx, 0);

        SparseSpaceType::Clear(mpb);
        TSystemVectorType& mb = *mpb;
        SparseSpaceType::Resize(mb, 0);

        // Setting to zero the internal flag to ensure that the dof sets are recalculated
        GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        GetBuilderAndSolver()->Clear();

        GetScheme()->Clear();

        if (this->GetEchoLevel() > 0)
        {
            if(BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
            {
                std::cout << "Newton Raphson strategy Clear function used" << std::endl;
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * This should be considered as a "post solution" convergence check which is useful for coupled analysis
     * the convergence criteria used is the one used inside the "solve" step
     */
    
    bool IsConverged()
    {
        KRATOS_TRY;

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;


        if (mpConvergenceCriteria->mActualizeRHSIsNeeded == true)
        {
            GetBuilderAndSolver()->BuildRHS(GetScheme(), BaseType::GetModelPart(), mb);
        }

        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();

        return mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
        
        KRATOS_CATCH("");
    }

    /**
     * This operations should be called before printing the results when non trivial results (e.g. stresses)
     * need to be calculated given the solution of the step
     * This operations should be called only when needed, before printing as it can involve a non negligible cost
     */
    
    void CalculateOutputData()
    {
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();
        GetScheme()->CalculateOutputData(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
    }

    /**
     * Performs all the required operations that should be done (for each step) before solving the solution step.
     * A member variable should be used as a flag to make sure this function is called only once per step.
     */
    
    void InitializeSolutionStep()
    {
        KRATOS_TRY;

        // Pointers needed in the solution
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

        // Set up the system, operation performed just once unless it is required
        // to reform the dof set at each iteration
        if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false || mReformDofSetAtEachStep == true)
        {
            // Setting up the list of the DOFs to be solved
            double setup_dofs_begintime = OpenMPUtils::GetCurrentTime();
            pBuilderAndSolver->SetUpDofSet(pScheme, BaseType::GetModelPart());
            if (this->GetEchoLevel() > 0 && rank == 0)
            {
                double setup_dofs_endtime = OpenMPUtils::GetCurrentTime();
                std::cout << "setup_dofs_time : " << setup_dofs_endtime- setup_dofs_begintime << std::endl;
            }

            // Shaping correctly the system
            double setup_system_begin = OpenMPUtils::GetCurrentTime();
            pBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());
            if (this->GetEchoLevel() > 0 && rank == 0)
            {
                double setup_system_end = OpenMPUtils::GetCurrentTime();
                std::cout << rank << ": setup_system_time : " << setup_system_end- setup_system_begin << std::endl;
            }
        }

        // Prints informations about the current time
        if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0 )
        {
            std::cout << " " << std::endl;
            std::cout << "CurrentTime = " << BaseType::GetModelPart().GetProcessInfo()[TIME] << std::endl;
        }

        if (mSolutionStepIsInitialized == false)
        {
            // typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
            // typename TSchemeType::Pointer pScheme = GetScheme();
            // int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

            // Setting up the Vectors involved to the correct size
            double system_matrix_resize_begin = OpenMPUtils::GetCurrentTime();
            pBuilderAndSolver->ResizeAndInitializeVectors(mpA, mpDx, mpb, BaseType::GetModelPart().Elements(), BaseType::GetModelPart().Conditions(), BaseType::GetModelPart().GetProcessInfo());
            if (this->GetEchoLevel() > 0 && rank == 0)
            {
                double system_matrix_resize_end = OpenMPUtils::GetCurrentTime();
                std::cout << rank << ": system_matrix_resize_time : " << system_matrix_resize_end- system_matrix_resize_begin << std::endl;
            }
            
            TSystemMatrixType& mA = *mpA;
            TSystemVectorType& mDx = *mpDx;
            TSystemVectorType& mb = *mpb;

            // Initial operations ... things that are constant over the Solution Step
            pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

            // Initial operations ... things that are constant over the Solution Step
            pScheme->InitializeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

            // TODO: Check if correct
            mpConvergenceAccelerator->InitializeSolutionStep();
            
            mSolutionStepIsInitialized = true;
        }

        KRATOS_CATCH("");
    }

    /**
     * Performs all the required operations that should be done (for each step) after solving the solution step.
     * A member variable should be used as a flag to make sure this function is called only once per step.
    */
    
    void FinalizeSolutionStep()
    {
        KRATOS_TRY;

        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        //Finalisation of the solution step,
        //operations to be done after achieving convergence, for example the
        //Final Residual Vector (mb) has to be saved in there
        //to avoid error accumulation

        pScheme->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);
        pBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

        //Cleaning memory after the solution
        pScheme->Clean();

        //reset flags for next step
        mSolutionStepIsInitialized = false;

        if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
        {
            SparseSpaceType::Clear(mpA);
            SparseSpaceType::Clear(mpDx);
            SparseSpaceType::Clear(mpb);

            this->Clear();
        }

        KRATOS_CATCH("");
    }
 
    /**
     * This is the initialization of the iteration cycle
     * @return is_converged: True when the problem has converged
     * @return ResidualIsUpdated: True when the residual has been updated
     * @return iteration_number: Number of iterations done 
     * @param pScheme: The integration scheme
     * @param pNewLinearSolver: The linear solver employed
     * @param pBuilderAndSolver: The builder and solver considered
     * @param rDofSet: The set of degrees of freedom
     * @param mA: The LHS of the problem
     * @return mDx: The solution to the problem
     * @param mb: The RHS of the problem
     */
    
    void InitiliazeCycle(
        bool& is_converged,
        bool & ResidualIsUpdated,
        unsigned int& iteration_number,
        typename TSchemeType::Pointer& pScheme,
        typename TBuilderAndSolverType::Pointer& pBuilderAndSolver,
        DofsArrayType& rDofSet,
        TSystemMatrixType& mA,
        TSystemVectorType& mDx,
        TSystemVectorType& mb
        )
    {
        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
      
        CoutSolvingProblem();
        
//         BaseType::GetModelPart().GetProcessInfo().SetNonLinearIterationNumber(iteration_number);

        pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
        is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        // Function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
        {
            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mDx);
            TSparseSpace::SetToZero(mb);

            pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }
        else
        {
            TSparseSpace::SetToZero(mDx); // mDx=0.00;
            TSparseSpace::SetToZero(mb);

            pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }

        if (this->GetEchoLevel() == 3) // If it is needed to print the debug info
        {
            std::cout << "SystemMatrix = " << mA << std::endl;
            std::cout << "solution obtained = " << mDx << std::endl;
            std::cout << "RHS  = " << mb << std::endl;
        }
        else if (this->GetEchoLevel() == 4) // Print to matrix market file
        {
            std::stringstream matrix_market_name;
            matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << iteration_number << ".mm";
            TSparseSpace::WriteMatrixMarketMatrix((char*) (matrix_market_name.str()).c_str(), mA, false);

            std::stringstream matrix_market_vectname;
            matrix_market_vectname << "b_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << iteration_number << ".mm.rhs";
            TSparseSpace::WriteMatrixMarketVector((char*) (matrix_market_vectname.str()).c_str(), mb);
        }

        // Update results
        rDofSet = pBuilderAndSolver->GetDofSet();
        pScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        // Move the mesh if needed
        if (BaseType::MoveMeshFlag() == true)
        {
            BaseType::MoveMesh();
        }

        pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

        if (is_converged == true)
        {
            // Initialisation of the convergence criteria
            rDofSet = pBuilderAndSolver->GetDofSet();
            mpConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

            if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
            {
                TSparseSpace::SetToZero(mb);

                pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb);
            }

            is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
        }
    }
 
    /**
     * This is the inner iteration cycle considered in the SolveSolutionStep
     * @return is_converged: True when the problem has converged
     * @return ResidualIsUpdated: True when the residual has been updated
     * @return iteration_number: Number of iterations done 
     * @param pScheme: The integration scheme
     * @param pNewLinearSolver: The linear solver employed
     * @param pBuilderAndSolver: The builder and solver considered
     * @param rDofSet: The set of degrees of freedom
     * @param mA: The LHS of the problem
     * @return mDx: The solution to the problem
     * @param mb: The RHS of the problem
     */
    
    void IterationCycle(
        bool& is_converged,
        bool & ResidualIsUpdated,
        unsigned int& iteration_number,
        typename TSchemeType::Pointer& pScheme,
        typename TBuilderAndSolverType::Pointer& pBuilderAndSolver,
        DofsArrayType& rDofSet,
        TSystemMatrixType& mA,
        TSystemVectorType& mDx,
        TSystemVectorType& mb
        )
    {
        // Iteration Cicle... performed only for NonLinearProblems
        while (is_converged == false && iteration_number++ < mMaxIterationNumber)
        {
            // Setting the number of iteration
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number; 
            
            CoutSolvingProblem();
           
            pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

            is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

            // Call the linear system solver to find the correction mDx for the
            // It is not called if there is no system to solve
            if (SparseSpaceType::Size(mDx) != 0)
            {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false )
                {
                    if( GetKeepSystemConstantDuringIterations() == false)
                    {
                        //mA = 0.00;
                        TSparseSpace::SetToZero(mA);
                        TSparseSpace::SetToZero(mDx);
                        TSparseSpace::SetToZero(mb);

                        pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
                    }
                    else
                    {
                        TSparseSpace::SetToZero(mDx);
                        TSparseSpace::SetToZero(mb);

                        pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
                    }
                }
                else
                {
                    TSparseSpace::SetToZero(mDx);
                    TSparseSpace::SetToZero(mb);

                    pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
                }
            }
            else
            {
                std::cout << "ATTENTION: no free DOFs!! " << std::endl;
            }
            
            // Updating the results stored in the database
            UpdateDatabase(pScheme, pBuilderAndSolver, rDofSet, mA, mDx, mb, BaseType::MoveMeshFlag());
            
            pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
            
            ResidualIsUpdated = false;

            if (is_converged == true)
            {
                if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {
                    TSparseSpace::SetToZero(mb);

                    pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb);
                    ResidualIsUpdated = true;
                    //std::cout << "mb is calculated" << std::endl;
                }
                is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
            }

            // The convergence accelerator is applied
            mpConvergenceCriteriaForConvergenceAccelerator->InitializeSolutionStep(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
            is_converged = ConvergenceAcceleratorStep(pScheme, pBuilderAndSolver, rDofSet, mA, mDx, mb); // TODO: Check if it is correct
        }
    }
    
    /**
     * Here the database is updated
     * @param pScheme: The integration scheme
     * @param pNewLinearSolver: The linear solver employed
     * @param pBuilderAndSolver: The builder and solver considered
     * @param rDofSet: The set of degrees of freedom
     * @param mA: The LHS of the problem
     * @return mDx: The solution to the problem
     */
    void UpdateDatabase(
        typename TSchemeType::Pointer& pScheme,
        typename TBuilderAndSolverType::Pointer& pBuilderAndSolver,
        DofsArrayType& rDofSet,
        TSystemMatrixType& mA,
        TSystemVectorType& mDx,
        TSystemVectorType& mb,
        const bool MoveMesh
    )
    {
        rDofSet = pBuilderAndSolver->GetDofSet();
        pScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        // Move the mesh if needed
        if (MoveMesh == true) 
        {
            BaseType::MoveMesh();
        }
    }
        
    /**
     * Here the accelerator is computed to get the new residual
     * @param pScheme: The integration scheme
     * @param pNewLinearSolver: The linear solver employed
     * @param pBuilderAndSolver: The builder and solver considered
     * @param rDofSet: The set of degrees of freedom
     * @param mA: The LHS of the problem
     * @return mDx: The solution to the problem
     */
    bool ConvergenceAcceleratorStep(
        typename TSchemeType::Pointer& pScheme,
        typename TBuilderAndSolverType::Pointer& pBuilderAndSolver,
        DofsArrayType& rDofSet,
        TSystemMatrixType& mA,
        TSystemVectorType& mDx,
        TSystemVectorType& mb
    )
    {
        // Here we start the non-linear iteration of the convergence accelerator
     
        bool is_converged = false;
        
        std::string CCForCAType = typeid(*mpConvergenceCriteriaForConvergenceAccelerator).name();
        
        std::cout << "\t--------------------------------" << std::endl;
        std::cout << "\tCONVERGENCE ACCELERATION STARTED" << std::endl;
        
        TSystemVectorType auxDx = mReductionCoefficient * mDx;
        
        unsigned int nl_it = 0; 
        while (is_converged == false && nl_it++ < mMaxNumberConvergenceAccelerationIterations)
        {            
            mpConvergenceAccelerator->InitializeNonLinearIteration();
                        
            is_converged = mpConvergenceCriteriaForConvergenceAccelerator->PreCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

            if (this->GetEchoLevel() != 0)
            {
                std::cout << "\tPerforming convergence acceleration non-linear iteration " << nl_it << std::endl;
            }
                
            if (CCForCAType.find("ResidualCriteria") != std::string::npos) // TODO: Check the name
            {
                // Calculate the new displacement
                mpConvergenceAccelerator->UpdateSolution(mb, auxDx);
                
                // Update the database
                UpdateDatabase(pScheme, pBuilderAndSolver, rDofSet, mA, auxDx, mb, BaseType::MoveMeshFlag());
                
                // Update residual variables
                if (mUpdateSystem == true)
                {
    //                 pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb);
                    pBuilderAndSolver->Build(pScheme, BaseType::GetModelPart(),mA, mb);
                }
                else
                {
                    mb -= prod(mA, auxDx);
                }
            }
            else if (CCForCAType.find("DisplacementConvergenceCriterion") != std::string::npos)
            {
                // Calculate the new displacement
                mpConvergenceAccelerator->UpdateSolution(mDx, auxDx);
                
                // Update the database
                UpdateDatabase(pScheme, pBuilderAndSolver, rDofSet, mA, auxDx, mb, BaseType::MoveMeshFlag());
                
                // Update residual variables
                mDx -= auxDx;
            }
            else if (CCForCAType.find("AndConvergenceCriterion") != std::string::npos) // TODO: Check the name
            {
                // Calculate the new displacement
                mpConvergenceAccelerator->UpdateSolution(mb, auxDx);
                mpConvergenceAccelerator->UpdateSolution(mDx, auxDx);
                
                // Update the database
                UpdateDatabase(pScheme, pBuilderAndSolver, rDofSet, mA, auxDx, mb, BaseType::MoveMeshFlag());
                
                // Update residual variables
                mDx -= auxDx;
                if (mUpdateSystem == true)
                {
    //                 pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb);
                    pBuilderAndSolver->Build(pScheme, BaseType::GetModelPart(),mA, mb);
                }
                else
                {
                    mb -= prod(mA, auxDx);
                }
            }
           
            mpConvergenceAccelerator->FinalizeNonLinearIteration(); 

            is_converged = mpConvergenceCriteriaForConvergenceAccelerator->PostCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);            
        }
        
        if (is_converged == true)
        {
            std::cout << "\tNON-LINEAR ITERATION FOR THE CONVERGENCE ACCELERATION ACHIEVED" << std::endl;
            std::cout << "\tTotal non-linear iterations "<< nl_it << std::endl;
            std::cout << "\t--------------------------------" << std::endl;
        }
        else
        {
            std::cout << "\tCONVERGENCE ACCELERATION FINISHED" << std::endl;
            std::cout << "\t--------------------------------" << std::endl;
        }
        
        mpConvergenceAccelerator->FinalizeSolutionStep();
        
        return is_converged;
    }
    
    /**
     * Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    
    bool SolveSolutionStep()
    {
        // Pointers needed in the solution
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();

        DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        //  Initializing the parameters of the Newton-Raphson cicle
        unsigned int iteration_number = 1;
        bool is_converged = false;
        bool ResidualIsUpdated = false;
        
        InitiliazeCycle(is_converged, ResidualIsUpdated, iteration_number, pScheme, pBuilderAndSolver, rDofSet, mA, mDx, mb);
        
        IterationCycle(is_converged, ResidualIsUpdated, iteration_number, pScheme, pBuilderAndSolver, rDofSet, mA, mDx, mb);

        // Plots a warning if the maximum number of iterations 
        if (iteration_number >= mMaxIterationNumber && is_converged == false  && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
        {            
            MaxIterationsExceeded();
        }

        // Recalculate residual if needed (note that some convergence criteria need it to be recalculated)
        if (ResidualIsUpdated == false)
        {
            // NOTE: The following part will be commented because it is time consuming
      
            //    TSparseSpace::SetToZero(mb);
            //    pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb);
        }

        // Calculate reactions if required
        if (mCalculateReactionsFlag == true)
        {
            pBuilderAndSolver->CalculateReactions(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }
        
        return is_converged;
    }

    ///@}
    ///@name Access
    ///@{

    TSystemMatrixType& GetSystemMatrix()
    {
        TSystemMatrixType& mA = *mpA;

        return mA;
    }

    void SetKeepSystemConstantDuringIterations(bool value)
    {
        mKeepSystemConstantDuringIterations = value;
    }

    bool GetKeepSystemConstantDuringIterations()
    {
        return mKeepSystemConstantDuringIterations;
    }
    
    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    
    typename TSchemeType::Pointer mpScheme;
    
    typename TConvergenceAcceleratorType::Pointer mpConvergenceAccelerator;
    
    unsigned int mMaxNumberConvergenceAccelerationIterations;
    
    typename TConvergenceCriteriaType::Pointer mpConvergenceCriteriaForConvergenceAccelerator;
    
    bool mUpdateSystem;
    
    double mReductionCoefficient;

    typename TLinearSolver::Pointer mpLinearSolver;

    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;

    typename TConvergenceCriteriaType::Pointer mpConvergenceCriteria;

    TSystemVectorPointerType mpDx;
    TSystemVectorPointerType mpb;
    TSystemMatrixPointerType mpA;

    /**
     * Flag telling if it is needed to reform the DofSet at each
     * solution step or if it is possible to form it just once
     * - true  => reforme at each time step
     * - false => form just one (more efficient)
     * Default = false
     */
    bool mReformDofSetAtEachStep;

    /**
    * Flag telling if it is needed or not to compute the reactions
    * default = true
    */
    bool mCalculateReactionsFlag;

    /**
    * Flag telling if the solution step is initialized
    */
    bool mSolutionStepIsInitialized;

    /**
    * Maximum number of iteratiosn, default 30
    */
    unsigned int mMaxIterationNumber;

    /**
    * Flag telling if the itialize was performed
    */
    bool mInitializeWasPerformed;

    /**
    * Flag to allow keeping system matrix constant during iterations
    */
    bool mKeepSystemConstantDuringIterations;

    ///@}
    ///@name Protected Operators
    ///@{

    void CoutSolvingProblem()
    {
        if (this->GetEchoLevel() != 0)
        {
            std::cout << "STEP: " << BaseType::GetModelPart().GetProcessInfo()[TIME_STEPS] << "\t NON LINEAR ITERATION: " << BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] << "\t TIME: " << BaseType::GetModelPart().GetProcessInfo()[TIME] << "\t DELTA TIME: " << BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME]  << std::endl;
        }
    }
    
    void MaxIterationsExceeded()
    {
        if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0 )
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "******* ATTENTION: max iterations exceeded ********" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }
    }

    /**
     * Function to perform expensive checks.
     * It is designed to be called ONCE to verify that the input is correct.
     */
    
    int Check()
    {
        KRATOS_TRY;

        BaseType::Check();

        GetBuilderAndSolver()->Check(BaseType::GetModelPart());

        GetScheme()->Check(BaseType::GetModelPart());

        mpConvergenceCriteria->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@{

    /** 
     * Copy constructor.
     */
    
    ResidualBasedNewtonRaphsonContactAcceleratedStrategy(const ResidualBasedNewtonRaphsonContactAcceleratedStrategy& Other)
    {
    };

private:

    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; /* Class ResidualBasedNewtonRaphsonContactAcceleratedStrategy */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  // namespace Kratos

#endif /* KRATOS_RESIDUALBASED_NEWTON_RAPHSON_CONTACT_ACCELERATED_STRATEGY */
