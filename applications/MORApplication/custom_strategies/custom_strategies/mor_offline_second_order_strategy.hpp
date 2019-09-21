//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(MOR_OFFLINE_SECOND_ORDER_STRATEGY)
#define MOR_OFFLINE_SECOND_ORDER_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/builtin_timer.h"
#include "utilities/qr_utility.h"

//default builder and solver
#include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"

namespace Kratos
{

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

/**
 * @class MorOfflineSecondOrderStrategy
 * @ingroup KratosCore
 * @brief This is the linear MOR matrix output strategy
 * @details This strategy builds the K and M matrices and outputs them
 * @author Aditya Ghantasala
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class MorOfflineSecondOrderStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@{
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MorOfflineSecondOrderStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef SystemMatrixBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef TDenseSpace DenseSpaceType;

    typedef typename TDenseSpace::MatrixType TDenseMatrixType;

    typedef typename TDenseSpace::MatrixPointerType TDenseMatrixPointerType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef TLinearSolver TLinearSolverType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    ///@}
    ///@name Life Cycle

    ///@{

    /**
     * Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration schemed
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    MorOfflineSecondOrderStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        bool MoveMeshFlag = false)
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, MoveMeshFlag)
    {
        KRATOS_TRY;
        
        // Saving the scheme
        mpScheme = pScheme;


        // Setting up the default builder and solver
         mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer(
            new TBuilderAndSolverType(pNewLinearSolver)); 

        // Saving the linear solver
        mpLinearSolver = pNewLinearSolver;            

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(false);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        // GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(0);

        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~MorOfflineSecondOrderStrategy() override
    {
        Clear();
    }

    /**
     * @brief Set method for the time scheme
     * @param pScheme The pointer to the time scheme considered
     */
    void SetScheme(typename TSchemeType::Pointer pScheme)
    {
        mpScheme = pScheme;
    };

    /**
     * @brief Get method for the time scheme
     * @return mpScheme: The pointer to the time scheme considered
     */
    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    /**
     * @brief Set method for the builder and solver
     * @param pNewBuilderAndSolver The pointer to the builder and solver considered
     */
    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    /**
     * @brief Get method for the builder and solver
     * @return mpBuilderAndSolver: The pointer to the builder and solver considered
     */
    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    /**
     * @brief Set method for the builder and solver
     * @param pNewBuilderAndSolver The pointer to the builder and solver considered
     */
    void SetLinearSolver(typename TLinearSolverType::Pointer pNewLinearSolver)
    {
        mpLinearSolver = pNewLinearSolver;
    };

    /**
     * @brief Get method for the builder and solver
     * @return mpBuilderAndSolver: The pointer to the builder and solver considered
     */
    typename TLinearSolverType::Pointer GetLinearSolver()
    {
        return mpLinearSolver;
    };

    /**
     * @brief This method sets the flag mInitializeWasPerformed
     * @param InitializePerformedFlag The flag that tells if the initialize has been computed
     */
    void SetInitializePerformedFlag(bool InitializePerformedFlag = true)
    {
        mInitializeWasPerformed = InitializePerformedFlag;
    }

    /**
     * @brief This method gets the flag mInitializeWasPerformed
     * @return mInitializeWasPerformed: The flag that tells if the initialize has been computed
     */
    bool GetInitializePerformedFlag()
    {
        return mInitializeWasPerformed;
    }

    /**
     * @brief This method sets the flag mReformDofSetAtEachStep
     * @param Flag The flag that tells if each time step the system is rebuilt
     */
    void SetReformDofSetAtEachStepFlag(bool Flag)
    {
        mReformDofSetAtEachStep = Flag;
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
    }

    /**
     * @brief This method returns the flag mReformDofSetAtEachStep
     * @return The flag that tells if each time step the system is rebuilt
     */
    bool GetReformDofSetAtEachStepFlag()
    {
        return mReformDofSetAtEachStep;
    }

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic informations
     * - 2: Printing linear solver data
     * - 3: Print of debug informations: Echo of stiffness matrix, Dx, b...
     */

    void SetEchoLevel(int Level) override
    {
        BaseType::mEchoLevel = Level;
        GetBuilderAndSolver()->SetEchoLevel(Level);
    }

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/

    /**
     * @brief Initialization of member variables and prior operations
     */
    virtual void Initialize() override
    {
        KRATOS_TRY;

        if (mInitializeWasPerformed == false)
        {
            mpA = TSparseSpace::CreateEmptyMatrixPointer();
            mpS = TSparseSpace::CreateEmptyMatrixPointer();
            mpM = TSparseSpace::CreateEmptyMatrixPointer();
            mpRHS = TSparseSpace::CreateEmptyVectorPointer();

            mpAr = TSparseSpace::CreateEmptyMatrixPointer();
            mpMr = TSparseSpace::CreateEmptyMatrixPointer();
            mpRHSr = TSparseSpace::CreateEmptyVectorPointer();
            mpBasis = TSparseSpace::CreateEmptyMatrixPointer();
            mpSr = TSparseSpace::CreateEmptyMatrixPointer();

            //pointers needed in the solution
            typename TSchemeType::Pointer p_scheme = GetScheme();

            //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            if (p_scheme->SchemeIsInitialized() == false)
                p_scheme->Initialize(BaseType::GetModelPart());

            //Initialize The Elements - OPERATIONS TO BE DONE ONCE
            if (p_scheme->ElementsAreInitialized() == false)
                p_scheme->InitializeElements(BaseType::GetModelPart());

            //Initialize The Conditions - OPERATIONS TO BE DONE ONCE
            if (p_scheme->ConditionsAreInitialized() == false)
                p_scheme->InitializeConditions(BaseType::GetModelPart());

            mInitializeWasPerformed = true;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Clears the internal storage
     */
    virtual void Clear() override
    {
        KRATOS_TRY;

        // if the preconditioner is saved between solves, it
        // should be cleared here.
        GetBuilderAndSolver()->GetLinearSystemSolver()->Clear();

        if (mpA != nullptr)
            SparseSpaceType::Clear(mpA);
        if (mpM != nullptr)
            SparseSpaceType::Clear(mpM);
        if (mpRHS != nullptr)
            SparseSpaceType::Clear(mpRHS);
        if (mpS != nullptr)
            SparseSpaceType::Clear(mpS);

        //setting to zero the internal flag to ensure that the dof sets are recalculated
        GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        GetBuilderAndSolver()->Clear();
        GetScheme()->Clear();

        mInitializeWasPerformed = false;

        KRATOS_CATCH("");
    }


    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    virtual void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        if (mSolutionStepIsInitialized == false)
        {
            //pointers needed in the solution
            typename TSchemeType::Pointer p_scheme = GetScheme();
            typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();

            const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

            //set up the system, operation performed just once unless it is required
            //to reform the dof set at each iteration
            BuiltinTimer system_construction_time;
            if (p_builder_and_solver->GetDofSetIsInitializedFlag() == false ||
                mReformDofSetAtEachStep == true)
            {
                //setting up the list of the DOFs to be solved
                BuiltinTimer setup_dofs_time;
                p_builder_and_solver->SetUpDofSet(p_scheme, BaseType::GetModelPart());
                KRATOS_INFO_IF("Setup Dofs Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << setup_dofs_time.ElapsedSeconds() << std::endl;

                //shaping correctly the system
                BuiltinTimer setup_system_time;
                p_builder_and_solver->SetUpSystem(BaseType::GetModelPart());
                KRATOS_INFO_IF("Setup System Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << setup_system_time.ElapsedSeconds() << std::endl;

                //setting up the Vectors involved to the correct size
                BuiltinTimer system_matrix_resize_time;
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpA, mpRHS, mpRHS,
                                                                 BaseType::GetModelPart());
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpM, mpRHS, mpRHS,
                                                                 BaseType::GetModelPart()); 
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpS, mpRHS, mpRHS,
                                                                 BaseType::GetModelPart()); 

                KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << system_matrix_resize_time.ElapsedSeconds() << std::endl;
            }

            KRATOS_INFO_IF("System Construction Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << system_construction_time.ElapsedSeconds() << std::endl;

            TSystemMatrixType& rA  = *mpA;
            TSystemMatrixType& rM = *mpM;
            TSystemVectorType& rRHS  = *mpRHS;

            //initial operations ... things that are constant over the Solution Step
            p_builder_and_solver->InitializeSolutionStep(BaseType::GetModelPart(), rA, rRHS, rRHS);
            p_builder_and_solver->InitializeSolutionStep(BaseType::GetModelPart(), rM, rRHS, rRHS);

            //initial operations ... things that are constant over the Solution Step
            p_scheme->InitializeSolutionStep(BaseType::GetModelPart(), rA, rRHS, rRHS);
            p_scheme->InitializeSolutionStep(BaseType::GetModelPart(), rM, rRHS, rRHS);

            mSolutionStepIsInitialized = true;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();

        TSystemMatrixType& rA  = *mpA;
        TSystemMatrixType& rM = *mpM;
        TSystemVectorType& rRHS  = *mpRHS;
        TSystemMatrixType& rS  = *mpS;

        //Finalisation of the solution step,
        //operations to be done after achieving convergence, for example the
        //Final Residual Vector (mb) has to be saved in there
        //to avoid error accumulation

        p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), rA, rRHS, rRHS);
        p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), rM, rRHS, rRHS);
        p_builder_and_solver->FinalizeSolutionStep(BaseType::GetModelPart(), rA, rRHS, rRHS);
        p_builder_and_solver->FinalizeSolutionStep(BaseType::GetModelPart(), rM, rRHS, rRHS);

        p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), rS, rRHS, rRHS);
        p_builder_and_solver->FinalizeSolutionStep(BaseType::GetModelPart(), rS, rRHS, rRHS);

        //Cleaning memory after the solution
        p_scheme->Clean();

        //reset flags for next step
        mSolutionStepIsInitialized = false;

        if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
        {
            SparseSpaceType::Clear(mpA);
            SparseSpaceType::Clear(mpRHS);
            SparseSpaceType::Clear(mpM);
            SparseSpaceType::Clear(mpS);

            this->Clear();
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_ERROR << "ERROR: Calling MOR offline strategy base class" << std::endl;
		return false;
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();

        GetBuilderAndSolver()->Check(BaseType::GetModelPart());

        GetScheme()->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Operators

    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access

    ///@{

    /**
     * @brief This method returns the LHS matrix
     * @return The LHS matrix
     */
    virtual TSystemMatrixType &GetSystemMatrix()
    {
        TSystemMatrixType &mA = *mpA;

        return mA;
    }

    virtual TSystemMatrixType &GetMassMatrix()
    {
        TSystemMatrixType &mM = *mpM;

        return mM;
    }

    virtual TSystemMatrixType &GetDampingMatrix()
    {
        TSystemMatrixType &mS = *mpS;

        return mS;
    }

    /**
     * @brief This method returns the RHS vector
     * @return The RHS vector
     */
    virtual TSystemVectorType& GetSystemVector()
    {
        TSystemVectorType& mb = *mpRHS;

        return mb;
    }

    virtual TSystemMatrixType& GetMr()
    {
        TSystemMatrixType &mMr = *mpMr;

        return mMr;
    }

    virtual TSystemMatrixPointerType& pGetMr()
    {
        return mpMr;
    }

    virtual TSystemMatrixType& GetKr()
    {
        TSystemMatrixType &mAr = *mpAr;

        return mAr;
    }

    virtual TSystemMatrixPointerType& pGetKr()
    {
        return mpAr;
    }

    virtual TSystemMatrixType& GetDr()
    {
        TSystemMatrixType &mSr = *mpSr;

        return mSr;
    };

    virtual TSystemMatrixPointerType& pGetDr()
    {
        return mpSr;
    }

    virtual TSystemVectorType& GetRHSr()
    {
        
        TSystemVectorType& mb = *mpRHSr;

        return mb;
    };

    virtual TSystemVectorPointerType& pGetRHSr()
    {
        return mpRHSr;
    }

    virtual TSystemMatrixType& GetBasis()
    {
        TSystemMatrixType &mBasis = *mpBasis;

        return mBasis;
    }

    virtual TSystemMatrixPointerType& pGetBasis()
    {
        return mpBasis;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

  private:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

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

    ///@}

  protected:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    typename TLinearSolver::Pointer mpLinearSolver; /// The pointer to the linear solver considered
    typename TSchemeType::Pointer mpScheme; /// The pointer to the time scheme employed
    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver; /// The pointer to the builder and solver employe

    TSystemVectorPointerType mpRHS; /// The RHS vector of the system of equations
    TSystemMatrixPointerType mpA; /// The Stiffness matrix of the system of equations
    TSystemMatrixPointerType mpM; /// The Mass matrix of the system of equations
    TSystemMatrixPointerType mpS; /// The Damping matrix of the system of equations


    // reduced matrices
    TSystemVectorPointerType mpRHSr; //reduced RHS
    TSystemMatrixPointerType mpAr;
    TSystemMatrixPointerType mpMr;
    // TDenseMatrixPointerType mpBasis;
    TSystemMatrixPointerType mpBasis;
    TSystemMatrixPointerType mpSr;

    int myTestInteger = 42;

    vector< double > mSamplingPoints;
    QR<double, row_major> mQR_decomposition;

    /**
     * @brief Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    * @details Default = false
        - true  : Reforme at each time step
        - false : Form just one (more efficient)
     */
    bool mReformDofSetAtEachStep;

    bool mSolutionStepIsInitialized; /// Flag to set as initialized the solution step

    bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief This method returns the components of the system of equations depending of the echo level
     * @param IterationNumber The non linear iteration in the solution loop
     */
    virtual void EchoInfo(const unsigned int IterationNumber)
    {
        TSystemMatrixType& rA  = *mpA;
        TSystemMatrixType& rM = *mpM;
        TSystemVectorType& rRHS  = *mpRHS;
        TSystemMatrixType& rS  = *mpS;

        if (this->GetEchoLevel() == 2) //if it is needed to print the debug info
        {
            KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        }
        else if (this->GetEchoLevel() == 3) //if it is needed to print the debug info
        {
            KRATOS_INFO("LHS") << "SystemMatrix = " << rA << std::endl;
            KRATOS_INFO("Mass")  << "Mass Matrix = " << rM << std::endl;
            KRATOS_INFO("Damping")  << "Damping Matrix = " << rS << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        }
        std::stringstream matrix_market_name;
        matrix_market_name << "K_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), rA, false);

        std::stringstream matrix_market_mass_name;
        matrix_market_mass_name << "M_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_mass_name.str()).c_str(), rM, false);   

        std::stringstream matrix_market_damping_name;
        matrix_market_damping_name << "D_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_damping_name.str()).c_str(), rS, false);        

        std::stringstream matrix_market_vectname;
        matrix_market_vectname << "RHS_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm.rhs";
        TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_vectname.str()).c_str(), rRHS);
    }

    /**
     * @brief This method prints information after reach the max number of iterations
     * @todo Replace by logger
     */

    virtual void MaxIterationsExceeded()
    {
        if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "******* ATTENTION: max iterations exceeded ********" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /**
     * Copy constructor.
     */

    MorOfflineSecondOrderStrategy(const MorOfflineSecondOrderStrategy &Other){};

    ///@}

}; /* Class MorOfflineSecondOrderStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* MOR_OFFLINE_SECOND_ORDER_STRATEGY  defined */