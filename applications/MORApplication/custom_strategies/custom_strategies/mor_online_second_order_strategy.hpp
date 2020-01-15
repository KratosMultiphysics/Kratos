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

#if !defined(MOR_ONLINE_SECOND_ORDER_STRATEGY)
#define MOR_ONLINE_SECOND_ORDER_STRATEGY

// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/solving_strategy.h"
#include "utilities/builtin_timer.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"

//default builder and solver
#include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"
#include "custom_strategies/custom_strategies/mor_offline_second_order_strategy.hpp"

// Application includes
#include "mor_application_variables.h"

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
 * @class MorOnlineSecondOrderStrategy
 * @ingroup KratosCore
 * @brief This is the MOR online strategy
 * @details The strategy calls an offline MOR strategy and solves the reduced system
 * @author Quirin Aumann
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver, //= LinearSolver<TSparseSpace,TDenseSpace>
          class OfflineStrategyType
          >
class MorOnlineSecondOrderStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@{
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MorOnlineSecondOrderStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef TUblasSparseSpace< std::complex< double > > ComplexSparseSpaceType;
    typedef TUblasDenseSpace< std::complex< double > > ComplexDenseSpaceType;

    typedef typename ComplexDenseSpaceType::MatrixType ComplexDenseMatrixType;
    typedef typename ComplexDenseSpaceType::VectorType ComplexDenseVectorType;
    typedef typename ComplexDenseSpaceType::MatrixPointerType ComplexDenseMatrixPointerType;
    typedef typename ComplexDenseSpaceType::VectorPointerType ComplexDenseVectorPointerType;

    // typedef MassAndStiffnessBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > TBuilderAndSolverType;

    // typedef LinearMorMatrixOutputStrategy< TSparseSpace, TDenseSpace, TLinearSolver > OfflineStrategyType;

    // typedef MorOfflineSecondOrderStrategy< TSparseSpace, TDenseSpace, TLinearSolver > OfflineStrategyType;

    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
    // typedef typename

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

    typedef typename TDenseSpace::MatrixPointerType LocalSystemMatrixPointerType;

    typedef typename TDenseSpace::VectorPointerType LocalSystemVectorPointerType;

    typedef LinearSolver<ComplexSparseSpaceType, ComplexDenseSpaceType> ComplexLinearSolverType;


    ///@}
    ///@name Life Cycle

    ///@{

    /**
     * Default constructor
     * TODO: rename to MorOnlineSecondOrderStrategy
     * @param rModelPart The model part of the problem
     * @param pLinearSolver The linear solver
     * @param pOfflineStrategy The offline MOR strategy
     */
    MorOnlineSecondOrderStrategy(
        ModelPart& rModelPart,
        typename TLinearSolver::Pointer pLinearSolver,
        typename OfflineStrategyType::Pointer pOfflineStrategy,
        bool expandSolution = true)
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, false)
    {
        KRATOS_TRY;

        // Setting up the default scheme
        // mpScheme = typename TSchemeType::Pointer(
        //     new TSchemeType());

        // Setting up the default builder and solver
        //TODO change to a builder and solver that can handle complex numbers
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer(
            new TBuilderAndSolverType(pLinearSolver));

        // Saving the linear solver
        mpLinearSolver = pLinearSolver;

        // Setting the offline strategy
        mpOfflineStrategy = pOfflineStrategy;

        // Set flags to correcty start the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;
        mReformDofSetAtEachStep = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(false);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // Do not rebuild the stiffness matrix
        this->SetRebuildLevel(0);

        // Set strategy flags
        mExpandSolution = expandSolution;

        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~MorOnlineSecondOrderStrategy() override
    {
        Clear();
    }

    /**
     * @brief Set method for the time scheme
     * @param pScheme The pointer to the time scheme considered
     */
    // void SetScheme(typename TSchemeType::Pointer pScheme)
    // {
    //     mpScheme = pScheme;
    // };

    /**
     * @brief Get method for the time scheme
     * @return mpScheme: The pointer to the time scheme considered
     */
    // typename TSchemeType::Pointer GetScheme()
    // {
    //     return mpScheme;
    // };

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

    void SetOfflineStrategy(typename OfflineStrategyType::Pointer pNewOfflineStrategy)
    {
        mpOfflineStrategy = pNewOfflineStrategy;
    };

    typename OfflineStrategyType::Pointer GetOfflineStrategy()
    {
        return mpOfflineStrategy;
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
        // GetBuilderAndSolver()->SetEchoLevel(Level);
    }

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/

    /**
     * @brief Initialization and solve of the offline strategy
     */
    void Initialize() override
    {
        KRATOS_TRY;

        if (mInitializeWasPerformed == false)
        {
            std::cout << "INITIALIZE online\n";
            typename OfflineStrategyType::Pointer p_offline_strategy = GetOfflineStrategy();
            p_offline_strategy->Initialize();
            p_offline_strategy->Solve();

            mpKr = p_offline_strategy->pGetKr();
            mpDr = p_offline_strategy->pGetDr();
            mpMr = p_offline_strategy->pGetMr();

            // cast rhs to complex
            LocalSystemVectorPointerType p_tmp_rhs = p_offline_strategy->pGetRHSr();
            KRATOS_WATCH(p_tmp_rhs)
            mpRHSr = ComplexDenseSpaceType::CreateEmptyVectorPointer();
            *mpRHSr = ComplexDenseVectorType( *p_tmp_rhs );

            mSystemSizeR = mpKr->size1();

            mpA = ComplexDenseSpaceType::CreateEmptyMatrixPointer();
            mpA->resize( mSystemSizeR, mSystemSizeR, false );
            mpDx = ComplexDenseSpaceType::CreateEmptyVectorPointer();
            mpDx->resize( mSystemSizeR, false );

            if (mExpandSolution)
            {
                mpBasis = p_offline_strategy->pGetBasis();
                mSystemSize = mpBasis->size1();
                mpResult = ComplexDenseSpaceType::CreateEmptyVectorPointer();
                mpResult->resize( mSystemSize, false );
            }

            mInitializeWasPerformed = true;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY;

        // if the preconditioner is saved between solves, it
        // should be cleared here.
        GetBuilderAndSolver()->GetLinearSystemSolver()->Clear();

        // if (mpKr != nullptr)
        //     SparseSpaceType::Clear(mpKr);
        // if (mpMr != nullptr)
        //     SparseSpaceType::Clear(mpMr);
        // if (mpDr != nullptr)
        //     SparseSpaceType::Clear(mpDr);
        // if (mpRHSr != nullptr)
        //     SparseSpaceType::Clear(mpRHSr);
        // if (mpA != nullptr)
        //     SparseSpaceType::Clear(mpA);
        // if (mpDx != nullptr)
        //     SparseSpaceType::Clear(mpDx);

        //setting to zero the internal flag to ensure that the dof sets are recalculated
        GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        GetBuilderAndSolver()->Clear();
        // GetScheme()->Clear();

        mInitializeWasPerformed = false;

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep() override
    {
        KRATOS_TRY;
        ModelPart& r_model_part = BaseType::GetModelPart();
        const int rank = r_model_part.GetCommunicator().MyPID();

        auto& r_process_info = r_model_part.GetProcessInfo();
        const double excitation_frequency = r_process_info[FREQUENCY];

        LocalSystemMatrixType& r_Kr = *mpKr;
        LocalSystemMatrixType& r_Dr = *mpDr;
        LocalSystemMatrixType& r_Mr = *mpMr;
        ComplexDenseVectorType& r_rhs = *mpRHSr;
        // ComplexDenseVectorType tmp = ComplexDenseVectorType( r_rhs );

        ComplexDenseMatrixType& r_A = *mpA;
        ComplexDenseVectorType& r_dx = *mpDx;

        BuiltinTimer system_construction_time;
        r_A = r_Dr;
        r_A *= std::complex<double>(0,excitation_frequency);
        r_A += r_Kr;
        r_A -= std::pow( excitation_frequency, 2.0 ) * r_Mr;

        KRATOS_INFO_IF("System Construction Time", BaseType::GetEchoLevel() > 0 && rank == 0)
            << system_construction_time.ElapsedSeconds() << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_TRY;

        // auto& r_process_info = BaseType::GetModelPart().GetProcessInfo();
        // const double excitation_frequency = r_process_info[FREQUENCY];

        // LocalSystemMatrixType& r_Kr = *mpKr;
        // LocalSystemMatrixType& r_Dr = *mpDr;
        // LocalSystemMatrixType& r_Mr = *mpMr;
        ComplexDenseVectorType& r_rhs = *mpRHSr;
        // ComplexDenseVectorType tmp = ComplexDenseVectorType( r_rhs );

        ComplexDenseMatrixType& r_A = *mpA;
        ComplexDenseVectorType& r_dx = *mpDx;

        // BuiltinTimer aa;
        // r_A = r_Dr;
        // r_A *= std::complex<double>(0,excitation_frequency);
        // r_A += r_Kr;
        // r_A -= std::pow( excitation_frequency, 2.0 ) * r_Mr;
        // KRATOS_INFO("solve system build") << aa.ElapsedSeconds() << "\n";
        // r_A = r_Kr - ( std::pow( excitation_frequency, 2.0 ) * r_Mr ); //TODO: damping is missing
        // r_A = r_Kr - ( std::pow( excitation_frequency, 2.0 ) * r_Mr ) + IMAG*excitation_frequency*r_Dr;

        // ComplexLinearSolverType 
        // TLinearSolver solver = GetBuilderAndSolver()->GetLinearSystemSolver();
        // solver.Solve( r_A, r_dx, r_rhs );
        BuiltinTimer ab;
        GetBuilderAndSolver()->GetLinearSystemSolver()->Solve( r_A, r_dx, r_rhs );
        KRATOS_INFO("solve system solve") << ab.ElapsedSeconds() << "\n";
        // KRATOS_WATCH(r_A)
        // KRATOS_WATCH(r_dx)
        // KRATOS_WATCH(r_rhs)

        std::cout << "ok hier\n";
        // GetBuilderAndSolver()->GetLinearSystemSolver()->Solve( r_A, r_dx, r_rhs ); //all have to be in the same space type!

        BuiltinTimer bb;
        if (mExpandSolution)
        {
            std::cout << "hello i am expanding\n";
            // TSystemMatrixType& r_basis = *mpBasis;
            ComplexDenseVectorType& r_result = *mpResult;
            r_result = prod( *mpBasis, *mpDx );
            // KRATOS_WATCH(r_result)
            AssignVariables(r_result);
        }
        else
        {
            /* code */
        }

        KRATOS_INFO("expand") << bb.ElapsedSeconds() << "\n";

        return true;
        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();
        KRATOS_WATCH(BaseType::GetModelPart())
        KRATOS_WATCH(GetBuilderAndSolver())
        GetBuilderAndSolver()->Check(BaseType::GetModelPart());

        // GetScheme()->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH("")
    }

    /**
     * @brief This method returns the components of the system of equations depending of the echo level
     */
    virtual void EchoInfo()
    {
        // TSystemMatrixType& rA  = *mpKr;
        // TSystemMatrixType& rM = *mpMr;
        // TSystemVectorType& rRHS  = *mpRHSr;
        // TSystemMatrixType& rS  = *mpDr;

        // if (this->GetEchoLevel() == 2) //if it is needed to print the debug info
        // {
        //     KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        // }
        // else if (this->GetEchoLevel() == 3) //if it is needed to print the debug info
        // {
        //     KRATOS_INFO("LHS") << "SystemMatrix = " << rA << std::endl;
        //     KRATOS_INFO("Dx")  << "Mass Matrix = " << mpMr << std::endl;
        //     KRATOS_INFO("Sx") << "Damping Matrix = " << rS << std::endl;
        //     KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        // }
        // std::stringstream matrix_market_name;
        // matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        // TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), rA, false);

        // std::stringstream matrix_market_mass_name;
        // matrix_market_mass_name << "M_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        // TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_mass_name.str()).c_str(), rM, false);

        // std::stringstream matrix_market_damping_name;
        // matrix_market_name << "S_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        // TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_damping_name.str()).c_str(), rS, false);

        // std::stringstream matrix_market_vectname;
        // matrix_market_vectname << "RHS_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << ".mm.rhs";
        // TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_vectname.str()).c_str(), *mpRHSr);
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
    // TSystemMatrixType &GetSystemMatrix()
    // {
    //     TSystemMatrixType &mA = *mpKr;

    //     return mA;
    // }

    /**
     * @brief This method returns the RHS vector
     * @return The RHS vector
     */
    // TSystemVectorType& GetSystemVector()
    // {
    //     TSystemVectorType& mb = *mpRHSr;

    //     return mb;
    // }

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
    // typename TSchemeType::Pointer mpScheme; /// The pointer to the time scheme employed
    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver; /// The pointer to the builder and solver employe
    typename OfflineStrategyType::Pointer mpOfflineStrategy;

    bool mExpandSolution;

    ComplexDenseMatrixPointerType mpA; /// The system matrix in reduced space
    ComplexDenseVectorPointerType mpDx; /// The result in reduced space
    ComplexDenseVectorPointerType mpRHSr; /// The RHS vector in reduced space
    LocalSystemMatrixPointerType mpKr; /// The Stiffness matrix in reduced space
    LocalSystemMatrixPointerType mpMr; /// The Mass matrix in reduced space
    LocalSystemMatrixPointerType mpDr; /// The Damping matrix in reduced space

    ComplexDenseVectorPointerType mpResult; /// The result expanded to original space
    // TSystemMatrixPointerType mpBasis;
    typename OfflineStrategyType::TReducedDenseMatrixPointerType mpBasis;

    unsigned int mSystemSize;
    unsigned int mSystemSizeR;

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

    void AssignVariables(ComplexDenseVectorType& rDisplacement, int step=0)
    {
        // KRATOS_WATCH(rDisplacement)
        auto& r_model_part = BaseType::GetModelPart();
        for( auto& node : r_model_part.Nodes() )
        {
            ModelPart::NodeType::DofsContainerType& rNodeDofs = node.GetDofs();

            for( auto it_dof = std::begin(rNodeDofs); it_dof != std::end(rNodeDofs); it_dof++ )
            {
                if( !(*it_dof)->IsFixed() )
                {
                    (*it_dof)->GetSolutionStepValue(step) = std::abs( rDisplacement((*it_dof)->EquationId()) );
                }
                else
                {
                    (*it_dof)->GetSolutionStepValue(step) = 0.0;
                }
            }
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

    MorOnlineSecondOrderStrategy(const MorOnlineSecondOrderStrategy &Other){};

    ///@}

}; /* Class MorOnlineSecondOrderStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* MOR_ONLINE_SECOND_ORDER_STRATEGY  defined */